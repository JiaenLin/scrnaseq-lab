// Reading .h5ad, in both layouts AnnData has used.
//
// Modern (AnnData ≥ 0.7): obs and var are groups, one dataset per column, with
// categoricals as a sub-group of categories + codes, and `encoding-type`
// attributes naming everything.
//
// Legacy (< 0.7): obs, var and obsm are each a single compound dataset, and a
// categorical column is bare integer codes whose levels live off in
// uns/<column>_categories. Files in this layout are still everywhere — the
// scanpy PBMC tutorial object is one — and an app that only reads the modern
// one fails on exactly the files people already have.

import {
  categorical, fromFactor, internGeneSet, maybeDiscrete, numeric, pickGeneAlias,
  type Column, type EmbeddingRef, type GeneAlias, type MatrixRef, type Scan,
} from './scan.ts'
import { classify, dropZeros, stridedSample, type CellMajor } from './matrix.ts'

/* eslint-disable @typescript-eslint/no-explicit-any */
type H5 = any

/** h5wasm returns attributes through getters; some builds wrap them. */
function attr(o: H5, name: string): unknown {
  const a = o?.attrs?.[name]
  if (a == null) return undefined
  return typeof a === 'object' && a !== null && 'value' in a ? (a as H5).value : a
}

const text = (v: unknown): string =>
  Array.isArray(v) ? String(v[0]) : v instanceof Uint8Array
    ? new TextDecoder().decode(v) : String(v)

/**
 * Which kind of HDF5 node this is, told by what it can do rather than by what
 * it is called.
 *
 * These used to compare `constructor.name` against 'Group' and 'Dataset',
 * which holds only for as long as nothing renames those classes. h5wasm is
 * bundled and minified into this app, so the names surviving is luck rather
 * than a guarantee — and the day it stopped being true the failure would be
 * silent: every group would read as "not a group", and obs would come back
 * empty with no error to explain it.
 */
const isGroup = (o: H5): boolean =>
  !!o && typeof o.keys === 'function' && typeof o.get === 'function'
const isDataset = (o: H5): boolean =>
  !!o && !isGroup(o) && (Array.isArray(o.shape) || typeof o.slice === 'function')

/** BigInt64Array and friends do not coerce, so every read goes through here. */
function toF32(v: unknown): Float32Array {
  if (v instanceof Float32Array) return v
  const a = v as ArrayLike<number | bigint>
  const out = new Float32Array(a.length)
  for (let i = 0; i < a.length; i++) out[i] = Number(a[i])
  return out
}

function toI32(v: unknown): Int32Array {
  if (v instanceof Int32Array) return v
  const a = v as ArrayLike<number | bigint>
  const out = new Int32Array(a.length)
  for (let i = 0; i < a.length; i++) out[i] = Number(a[i])
  return out
}

const toStrings = (v: unknown): string[] =>
  (v as unknown[]).map(x => (x instanceof Uint8Array ? new TextDecoder().decode(x) : String(x)))

export class H5adError extends Error {}

/**
 * The largest .h5ad the in-tab reader can accept.
 *
 * 2^31 - 1: h5wasm's emulated filesystem stores lengths as signed 32-bit.
 * Verified empirically against this build — a 2047 MB write succeeds and a
 * 2048 MB write throws.
 */
export const MAX_H5AD_BYTES = 2 ** 31 - 1

/** The compound-dataset layout: one array of rows, each row an array of fields. */
interface Compound { names: string[]; rows: unknown[][]; members: H5[] }

function compound(ds: H5): Compound | null {
  const members = ds?.metadata?.compound_type?.members
  if (!members) return null
  return { names: members.map((m: H5) => m.name as string), rows: ds.value as unknown[][], members }
}

const fieldOf = (c: Compound, name: string): unknown[] | null => {
  const i = c.names.indexOf(name)
  return i < 0 ? null : c.rows.map(r => r[i])
}

/**
 * Read obs into columns, whichever layout it is in.
 *
 * `uns` is passed because in the legacy layout that is the only place a
 * categorical column's level names exist — the column itself is bare integers,
 * and reading it without them turns "B cells" into 2.
 */
function readObs(obs: H5, uns: H5, nCells: number, note: (s: string) => void): Column[] {
  const out: Column[] = []

  const comp = compound(obs)
  if (comp) {
    for (let i = 0; i < comp.names.length; i++) {
      const name = comp.names[i]
      if (name === 'index' || name === '_index') continue
      const raw = comp.rows.map(r => r[i])
      const cats = uns && !(uns.get?.(`${name}_categories`) == null)
        ? toStrings(uns.get(`${name}_categories`).value) : null
      if (cats) {
        out.push(fromFactor(name, toI32(raw.map(Number)), cats))
      } else if (typeof raw[0] === 'string') {
        out.push(categorical(name, raw.map(String)))
      } else {
        out.push(maybeDiscrete(numeric(name, toF32(raw.map(Number))), nCells))
      }
    }
    return out
  }

  if (!isGroup(obs)) throw new H5adError('this file has no obs table, so there is nothing to group cells by')

  const indexName = text(attr(obs, '_index') ?? '_index')
  const order = attr(obs, 'column-order')
  const keys: string[] = order ? toStrings(order) : obs.keys().filter((k: string) => k !== indexName)

  for (const key of keys) {
    if (key === indexName) continue
    const o = obs.get(key)
    if (!o) continue
    try {
      if (isGroup(o)) {
        // A modern categorical: categories + codes, with -1 for NA.
        const cats = o.get('categories')
        const codes = o.get('codes')
        if (cats && codes) {
          out.push(fromFactor(key, toI32(codes.value), toStrings(cats.value)))
          continue
        }
        note(`obs/${key} is a group this reader does not understand; skipped`)
        continue
      }
      const v = o.value
      if (v?.length !== nCells) {
        note(`obs/${key} has ${v?.length ?? 0} values for ${nCells} cells; skipped`)
        continue
      }
      if (typeof v[0] === 'string' || v[0] instanceof Uint8Array) {
        out.push(categorical(key, toStrings(v)))
      } else {
        out.push(maybeDiscrete(numeric(key, toF32(v)), nCells))
      }
    } catch {
      note(`obs/${key} could not be read; skipped`)
    }
  }
  return out
}

/** Gene names from a var table, either layout. */
function readVar(v: H5, label: string): string[] {
  const comp = compound(v)
  if (comp) {
    const idx = fieldOf(comp, 'index') ?? fieldOf(comp, '_index')
    if (idx) return toStrings(idx)
  }
  if (isGroup(v)) {
    const name = text(attr(v, '_index') ?? '_index')
    const ds = v.get(name) ?? v.get('_index') ?? v.get('index')
    if (ds) return toStrings(ds.value)
  }
  throw new H5adError(`${label} has no gene names`)
}

/**
 * Every string-valued var column except the index.
 *
 * Read whole rather than sampled: a var table is one row per gene, so even an
 * atlas's is a few megabytes, and the alias is written into the bundle in full
 * anyway. Categoricals are expanded here — on this atlas the symbols are stored
 * as 31 017 categories over 31 053 codes, and reading only the categories would
 * misalign every gene after the first duplicate.
 */
function readVarStringColumns(v: H5, nGenes: number, note: (s: string) => void): {
  name: string; values: string[]
}[] {
  const out: { name: string; values: string[] }[] = []
  const comp = compound(v)
  if (comp) {
    for (const name of comp.names) {
      if (name === 'index' || name === '_index') continue
      const f = fieldOf(comp, name)
      if (!f || f.length !== nGenes) continue
      if (typeof f[0] !== 'string' && !(f[0] instanceof Uint8Array)) continue
      out.push({ name, values: toStrings(f) })
    }
    return out
  }
  if (!isGroup(v)) return out

  const indexName = text(attr(v, '_index') ?? '_index')
  for (const key of v.keys()) {
    if (key === indexName || key === '_index' || key === 'index') continue
    try {
      const o = v.get(key)
      if (!o) continue
      if (isGroup(o)) {
        const cats = o.get('categories')
        const codes = o.get('codes')
        if (!cats || !codes) continue
        const levels = toStrings(cats.value)
        const c = toI32(codes.value)
        if (c.length !== nGenes) continue
        out.push({ name: key, values: Array.from(c, i => (i >= 0 ? levels[i] ?? '' : '')) })
        continue
      }
      const val = o.value
      if (!val || val.length !== nGenes) continue
      if (typeof val[0] !== 'string' && !(val[0] instanceof Uint8Array)) continue
      out.push({ name: key, values: toStrings(val) })
    } catch {
      note(`var/${key} could not be read; skipped`)
    }
  }
  return out
}

/** Symbols beside accessions, or the reverse — whichever the var table adds. */
function readVarAlias(
  v: H5, genes: string[], note: (s: string) => void,
): GeneAlias | null {
  try {
    return pickGeneAlias(genes, readVarStringColumns(v, genes.length, note))
  } catch (e) {
    note(`the var table's other gene names could not be read (${(e as Error).message})`)
    return null
  }
}

/** A matrix, dense or sparse, in either the modern or the h5sparse layout. */
function readMatrix(m: H5, nCells: number, nGenes: number, key: string): CellMajor {
  if (isDataset(m)) {
    // Dense, cells × genes. Materialized whole, then thinned to nonzeros.
    const flat = toF32(m.value)
    const shape = m.shape as number[]
    if (shape[0] !== nCells || shape[1] !== nGenes) {
      throw new H5adError(`${key} is ${shape.join('×')}, expected ${nCells}×${nGenes}`)
    }
    let nnz = 0
    for (let i = 0; i < flat.length; i++) if (flat[i] !== 0) nnz++
    const indptr = new Int32Array(nCells + 1)
    const indices = new Int32Array(nnz)
    const data = new Float32Array(nnz)
    let at = 0
    for (let c = 0; c < nCells; c++) {
      indptr[c] = at
      const base = c * nGenes
      for (let g = 0; g < nGenes; g++) {
        const v = flat[base + g]
        if (v === 0) continue
        indices[at] = g
        data[at] = v
        at++
      }
    }
    indptr[nCells] = at
    return { nCells, nGenes, indptr, indices, data }
  }

  const enc = text(attr(m, 'encoding-type') ?? attr(m, 'h5sparse_format') ?? '')
  const data = toF32(m.get('data').value)
  const indices = toI32(m.get('indices').value)
  const indptr = toI32(m.get('indptr').value)

  if (enc.includes('csc')) {
    // Stored gene-major already; flip it to cell-major so one path follows.
    if (indptr.length !== nGenes + 1) {
      throw new H5adError(`${key} says CSC but its indptr is ${indptr.length}, not ${nGenes + 1}`)
    }
    const counts = new Int32Array(nCells + 1)
    for (let k = 0; k < indices.length; k++) counts[indices[k] + 1]++
    for (let c = 0; c < nCells; c++) counts[c + 1] += counts[c]
    const cursor = counts.slice(0, nCells)
    const ci = new Int32Array(indices.length)
    const cd = new Float32Array(data.length)
    for (let g = 0; g < nGenes; g++) {
      for (let k = indptr[g]; k < indptr[g + 1]; k++) {
        const at = cursor[indices[k]]++
        ci[at] = g
        cd[at] = data[k]
      }
    }
    return dropZeros({ nCells, nGenes, indptr: counts, indices: ci, data: cd })
  }

  if (indptr.length !== nCells + 1) {
    throw new H5adError(`${key} has an indptr of ${indptr.length} for ${nCells} cells — is it transposed?`)
  }
  return dropZeros({ nCells, nGenes, indptr, indices, data })
}


/** Row offsets of a CSR matrix, as counts per cell. */
function rowLengths(m: H5, nCells: number): Int32Array {
  const indptr = toI32(m.get('indptr').value)
  const out = new Int32Array(nCells)
  for (let c = 0; c < nCells; c++) out[c] = indptr[c + 1] - indptr[c]
  return out
}

/** Values read per slice. Large enough that chunk decompression amortizes. */
const SCAN_BLOCK = 8_000_000

/**
 * Read disjoint cell subsets from a CSR matrix in one forward pass.
 *
 * A window slides forward over the stored values, never backwards, so each
 * compressed chunk is decompressed at most once. Cells not in any group are
 * skipped outright — the window jumps to the next selected row — so a split
 * that touches a quarter of the cells reads far less than the whole file.
 */
function readGroupsCSR(
  m: H5, groups: Int32Array[], nCells: number, nGenes: number,
  onProgress?: (frac: number) => void,
): CellMajor[] {
  const indptr = toI32(m.get('indptr').value)
  const dataDs = m.get('data')
  const idxDs = m.get('indices')
  const total = indptr[nCells]

  // Which group each cell belongs to, and how big each group's arrays must be.
  const groupOf = new Int32Array(nCells).fill(-1)
  groups.forEach((cells, g) => { for (const c of cells) groupOf[c] = g })
  const out = groups.map(cells => {
    let nnz = 0
    for (const c of cells) nnz += indptr[c + 1] - indptr[c]
    return {
      nCells: cells.length, nGenes,
      indptr: new Int32Array(cells.length + 1),
      indices: new Int32Array(nnz),
      data: new Float32Array(nnz),
    } as CellMajor
  })
  const at = new Int32Array(groups.length)      // write cursor per group
  const row = new Int32Array(groups.length)     // next row index per group

  let winStart = 0, winEnd = 0
  let winData: Float32Array = new Float32Array(0)
  let winIdx: Int32Array = new Int32Array(0)

  for (let c = 0; c < nCells; c++) {
    const g = groupOf[c]
    if (g < 0) continue
    const p0 = indptr[c], p1 = indptr[c + 1]
    if (p0 < winStart || p1 > winEnd) {
      winStart = p0
      winEnd = Math.min(total, Math.max(p0 + SCAN_BLOCK, p1))
      winData = toF32(dataDs.slice([[winStart, winEnd]]))
      winIdx = toI32(idxDs.slice([[winStart, winEnd]]))
      onProgress?.(p0 / total)
    }
    const o = out[g]
    o.indptr[row[g]] = at[g]
    for (let k = p0; k < p1; k++) {
      const v = winData[k - winStart]
      if (v === 0) continue
      o.indices[at[g]] = winIdx[k - winStart]
      o.data[at[g]] = v
      at[g]++
    }
    row[g]++
  }
  // Explicit zeros were skipped, so the arrays may be shorter than reserved.
  return out.map((o, g) => {
    o.indptr[row[g]] = at[g]
    return at[g] === o.data.length ? o : {
      ...o, indices: o.indices.subarray(0, at[g]), data: o.data.subarray(0, at[g]),
    }
  })
}

/**
 * Read an embedding once and keep it.
 *
 * A 292 k-cell coordinate array is 2.3 MB — nothing beside the matrix — and a
 * split asks for it once per part per embedding. At 43 parts and two
 * embeddings that is 86 traversals of the same dataset through a block cache
 * that the matrix pass keeps evicting.
 */
function once(e: EmbeddingRef): EmbeddingRef {
  let held: Float32Array | null = null
  return { ...e, load: () => (held ??= e.load()) }
}

/** Enough values to classify a matrix, without reading all of it. */
function sampleMatrix(m: H5, want = 20000): Float64Array {
  if (isDataset(m)) {
    const shape = m.shape as number[]
    const rows = Math.min(shape[0], Math.max(1, Math.ceil(want / shape[1])))
    const step = Math.max(1, Math.floor(shape[0] / rows))
    const vals: number[] = []
    for (let r = 0; r < rows; r++) {
      const at = Math.min(shape[0] - 1, r * step)
      const row = m.slice([[at, at + 1], [0, shape[1]]]) as ArrayLike<number>
      for (let i = 0; i < row.length; i++) vals.push(Number(row[i]))
    }
    return Float64Array.from(vals)
  }
  const d = m.get('data')
  const n = (d.shape as number[])[0]
  // One strided read would mean n round trips into WASM; a few contiguous
  // blocks spread across the vector cost far less and sample just as fairly.
  const blocks = 8
  const per = Math.max(1, Math.floor(Math.min(want, n) / blocks))
  const vals: number[] = []
  for (let b = 0; b < blocks && vals.length < want; b++) {
    const at = Math.min(n - per, Math.floor((n / blocks) * b))
    if (at < 0) break
    const block = d.slice([[at, at + per]]) as ArrayLike<number>
    for (let i = 0; i < block.length; i++) vals.push(Number(block[i]))
  }
  return vals.length ? Float64Array.from(vals) : stridedSample(n, want, i => Number(d.slice([[i, i + 1]])[0]))
}

export interface H5Module {
  ready: Promise<{ FS: H5 }>
  File: new (name: string, mode: string) => H5
}

/** Read an .h5ad far enough to ask the user which columns to use. */
/**
 * Scan an already-open HDF5 handle.
 *
 * Split from the entry points below because how the file got opened — copied
 * into memory, or mounted lazily from a File inside a worker — changes nothing
 * about how it is read.
 */
export function scanOpenH5ad(f: H5, filename: string): Scan {
  const notes: string[] = []
  const note = (s: string) => notes.push(s)

  const obsRaw = f.get('obs')
  if (!obsRaw) throw new H5adError('this HDF5 file has no /obs — it is not an .h5ad')

  const X = f.get('X')
  if (!X) throw new H5adError('this file has no /X matrix')
  const shape = isDataset(X)
    ? (X.shape as number[])
    : (toI32(attr(X, 'shape') ?? attr(X, 'h5sparse_shape') ?? []) as unknown as number[])
  const nCells = Number(shape[0])
  const nGenes = Number(shape[1])
  if (!nCells || !nGenes) throw new H5adError('the /X matrix has no readable shape')

  const legacy = compound(obsRaw) != null
  note(legacy
    ? 'read in the legacy AnnData layout — obs and obsm are compound datasets and categorical levels come from uns'
    : 'read in the modern AnnData layout')

  const uns = f.get('uns')
  const columns = readObs(obsRaw, uns, nCells, note)

  // ---- gene lists -------------------------------------------------------
  const geneSets: Record<string, string[]> = {}
  const geneAliases: Record<string, GeneAlias> = {}
  const varDs = f.get('var')
  const varGenes = readVar(varDs, 'var')
  const varSet = internGeneSet(geneSets, varGenes, 'var')
  const varAlias = readVarAlias(varDs, varGenes, note)
  if (varAlias) geneAliases[varSet] = varAlias
  const rawVarDs = f.get('raw/var') ?? f.get('raw.var')
  let rawSet: string | null = null
  if (rawVarDs) {
    try {
      const rawGenes = readVar(rawVarDs, 'raw.var')
      rawSet = internGeneSet(geneSets, rawGenes, 'raw')
      // internGeneSet reuses an identical list, so raw and var can land on the
      // same key — do not overwrite an alias that is already there.
      if (!geneAliases[rawSet]) {
        const a = readVarAlias(rawVarDs, rawGenes, note)
        if (a) geneAliases[rawSet] = a
      }
    } catch { /* no raw genes */ }
  }

  // ---- candidate matrices ----------------------------------------------
  const matrices: MatrixRef[] = []
  const add = (key: string, m: H5, geneSet: string, cells: number, genes: number) => {
    if (!m) return
    try {
      const csr = !isDataset(m) &&
        !text(attr(m, 'encoding-type') ?? attr(m, 'h5sparse_format') ?? '').includes('csc')
      matrices.push({
        key, geneSet, nCells: cells, nGenes: genes,
        kind: classify(sampleMatrix(m)),
        load: () => readMatrix(m, cells, genes, key),
        ...(csr ? {
          nnzPerCell: () => rowLengths(m, cells),
          loadGroups: (groups, onProgress) =>
            readGroupsCSR(m, groups, cells, genes, onProgress),
        } : {}),
      })
    } catch (e) {
      note(`${key} could not be read (${(e as Error).message}); skipped`)
    }
  }

  add('X', X, varSet, nCells, nGenes)
  const layers = f.get('layers')
  if (isGroup(layers)) {
    for (const k of layers.keys()) add(`layers["${k}"]`, layers.get(k), varSet, nCells, nGenes)
  }
  const rawX = f.get('raw/X') ?? f.get('raw.X')
  if (rawX && rawSet) {
    add('raw.X', rawX, rawSet, nCells, geneSets[rawSet].length)
  }

  // ---- embeddings -------------------------------------------------------
  const embeddings: EmbeddingRef[] = []
  const obsm = f.get('obsm')
  const comp = obsm ? compound(obsm) : null
  if (comp) {
    comp.members.forEach((m: H5, i: number) => {
      const width = m.array_type?.shape?.[0] ?? 0
      if (width < 2) return
      embeddings.push(once({
        key: m.name, nDims: width,
        load: () => {
          const xy = new Float32Array(nCells * 2)
          for (let c = 0; c < nCells; c++) {
            const row = comp.rows[c][i] as ArrayLike<number>
            xy[2 * c] = Number(row[0])
            xy[2 * c + 1] = Number(row[1])
          }
          return xy
        },
      }))
    })
  } else if (isGroup(obsm)) {
    for (const k of obsm.keys()) {
      const ds = obsm.get(k)
      const s = ds?.shape as number[] | undefined
      if (!s || s.length < 2 || s[1] < 2 || s[0] !== nCells) continue
      embeddings.push(once({
        key: k, nDims: s[1],
        load: () => {
          const flat = toF32(ds.slice([[0, nCells], [0, 2]]))
          return flat.length === nCells * 2 ? flat : toF32(ds.value).filter((_v, i) => i % s[1] < 2)
        },
      }))
    }
  }

  // ---- provenance -------------------------------------------------------
  const unsKeys: string[] = isGroup(uns) ? uns.keys() : []
  const obsNames = columns.map(c => c.name)
  const provenance = {
    normalization: unsKeys.includes('log1p') ? 'log1p(CP10K)' : null,
    integration: embeddings.some(e => /harmony/i.test(e.key)) ? 'harmony' : null,
    doublets: obsNames.find(n => /doublet|scrublet/i.test(n)) ?? null,
    ambient: unsKeys.find(k => /soupx|cellbender|decontx/i.test(k)) ?? null,
  }

  return {
    format: 'h5ad', source: filename, nCells, columns, matrices, embeddings,
    geneSets, geneAliases, provenance, notes,
  }
}

// ---------------------------------------------------------------------------
// Two ways in.

/**
 * Copy the whole file into h5wasm's in-memory filesystem, then scan it.
 *
 * Simple and synchronous once the bytes are in hand, but it holds the file
 * twice and cannot exceed MAX_H5AD_BYTES. Used by the Node tests; the app now
 * takes the lazy path below.
 */
export async function scanH5ad(
  h5wasm: H5Module, bytes: Uint8Array, filename: string,
): Promise<Scan> {
  const { FS } = await h5wasm.ready

  // h5wasm is a 32-bit WebAssembly build of HDF5, and its emulated filesystem
  // holds the whole file. Lengths there are signed 32-bit, so 2 GB is a hard
  // ceiling — measured exactly: 2047 MB writes, 2048 MB does not.
  if (bytes.length >= MAX_H5AD_BYTES) {
    throw new H5adError(
      `this .h5ad is ${(bytes.length / 1e9).toFixed(1)} GB, which is past the 2 GB limit of the in-memory reader. Open it in the app, which mounts the file instead of copying it, or run tools/export_h5ad.py offline.`)
  }

  const path = `scan_${filename.replace(/[^\w.-]/g, '_')}`
  try { FS.unlink(path) } catch { /* first run */ }
  try {
    FS.writeFile(path, bytes)
  } catch {
    throw new H5adError(
      `this .h5ad (${(bytes.length / 1e9).toFixed(1)} GB) would not fit in the reader's memory.`)
  }
  return scanOpenH5ad(openOrFail(h5wasm, path), filename)
}

function openOrFail(h5wasm: H5Module, path: string): H5 {
  try {
    return new h5wasm.File(path, 'r')
  } catch {
    throw new H5adError('this file is not readable as HDF5. An .h5ad is an HDF5 file — if this came from Seurat, save it as .rds instead.')
  }
}

/**
 * Mount a File lazily and scan it, at any size.
 *
 * WORKERFS is Emscripten's read-only filesystem over Blob/File objects: it
 * serves reads straight from the File with FileReaderSync rather than copying
 * it into memory, so nothing is ever held except the bytes HDF5 actually asks
 * for. That is what lifts the 2 GB ceiling — a 9.3 GB atlas opens in
 * milliseconds because opening it reads only the superblock and the metadata.
 *
 * FileReaderSync exists only inside a worker, so this must run in one.
 */
/** Bytes fetched per block. Big enough that the per-read cost disappears. */
const BLOCK = 4 << 20
/** Blocks kept resident: 48 x 4 MB, so the cache is bounded at ~192 MB. */
const MAX_BLOCKS = 48

/**
 * Present a File to HDF5 as a file it can read, without copying it anywhere.
 *
 * Emscripten ships WORKERFS for exactly this, and it is unusable here: it
 * answers every read with a fresh Blob slice through FileReaderSync, and that
 * call costs about 4.4 ms however few bytes it returns. HDF5 walks its metadata
 * in thousands of small reads, so a scan that takes 2 s against a local file
 * had still not finished after ten minutes through WORKERFS. Measured on the
 * 9.3 GB atlas: 2000 reads of 4 KB took 8.8 s - 0.9 MB/s - while one large read
 * ran at 20 MB/s. The cost is per call, not per byte.
 *
 * So reads are served out of 4 MB blocks with an LRU behind them. HDF5's small
 * reads then land in memory, and the 4.4 ms is paid once per block instead of
 * once per read.
 *
 * FileReaderSync exists only inside a worker, so this must run in one.
 */
function mountBlockCached(FS: H5, file: File, dir: string, name: string): string {
  const fr = new FileReaderSync()
  const cache = new Map<number, Uint8Array>()

  const blockAt = (i: number): Uint8Array => {
    const hit = cache.get(i)
    if (hit) {
      cache.delete(i)        // reinsert, so iteration order is least-recent first
      cache.set(i, hit)
      return hit
    }
    const start = i * BLOCK
    const slice = file.slice(start, Math.min(file.size, start + BLOCK))
    const block = new Uint8Array(fr.readAsArrayBuffer(slice) as ArrayBuffer)
    cache.set(i, block)
    if (cache.size > MAX_BLOCKS) cache.delete(cache.keys().next().value as number)
    return block
  }

  try { FS.mkdir(dir) } catch { /* left from a previous file */ }
  try { FS.unlink(`${dir}/${name}`) } catch { /* first time */ }

  const node = FS.createFile(dir, name, {}, true, false)
  // stat() reports `usedBytes`, not `size`. Setting only the latter leaves HDF5
  // looking at what it believes is an empty file — which it does not report as
  // an error, it just never finds a superblock and reads nothing, forever.
  Object.defineProperty(node, 'usedBytes', { get: () => file.size })
  node.size = file.size
  node.stream_ops = {
    llseek(stream: H5, offset: number, whence: number): number {
      const base = whence === 1 ? stream.position : whence === 2 ? node.size : 0
      const pos = base + offset
      if (pos < 0) throw new FS.ErrnoError(28)
      return pos
    },
    read(_s: H5, buffer: Uint8Array, offset: number, length: number, position: number): number {
      if (position >= node.size) return 0
      const end = Math.min(node.size, position + length)
      let pos = position, done = 0
      while (pos < end) {
        const i = Math.floor(pos / BLOCK)
        const block = blockAt(i)
        const from = pos - i * BLOCK
        const n = Math.min(end - pos, block.length - from)
        if (n <= 0) break
        buffer.set(block.subarray(from, from + n), offset + done)
        done += n
        pos += n
      }
      return done
    },
  }
  return `${dir}/${name}`
}

/**
 * Open a File lazily and scan it, at any size.
 *
 * Nothing is held but the blocks HDF5 actually asks for, so the 2 GB ceiling of
 * the copy-it-all path does not apply here.
 */
export async function scanH5adFile(h5wasm: H5Module, file: File): Promise<{
  scan: Scan; handle: H5
}> {
  const { FS } = await h5wasm.ready
  if (typeof FileReaderSync === 'undefined') {
    throw new H5adError('a file this large can only be read inside a worker, and this is not one.')
  }
  const safe = file.name.replace(/[^\w.-]/g, '_')
  const handle = openOrFail(h5wasm, mountBlockCached(FS, file, '/lazy', safe))
  return { scan: scanOpenH5ad(handle, file.name), handle }
}
