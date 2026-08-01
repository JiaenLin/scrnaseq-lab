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
  categorical, fromFactor, internGeneSet, maybeDiscrete, numeric,
  type Column, type EmbeddingRef, type MatrixRef, type Scan,
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

const isGroup = (o: H5): boolean => o?.constructor?.name === 'Group'
const isDataset = (o: H5): boolean => o?.constructor?.name === 'Dataset'

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
export async function scanH5ad(
  h5wasm: H5Module, bytes: Uint8Array, filename: string,
): Promise<Scan> {
  const { FS } = await h5wasm.ready
  // Written relative to the emulated working directory: the root of h5wasm's
  // in-memory filesystem is not writable in its Node build.
  const path = `scan_${filename.replace(/[^\w.-]/g, '_')}`
  try { FS.unlink(path) } catch { /* first run */ }
  FS.writeFile(path, bytes)

  let f: H5
  try {
    f = new h5wasm.File(path, 'r')
  } catch {
    throw new H5adError('this file is not readable as HDF5. An .h5ad is an HDF5 file — if this came from Seurat, save it as .rds instead.')
  }

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
  const varSet = internGeneSet(geneSets, readVar(f.get('var'), 'var'), 'var')
  const rawVarDs = f.get('raw/var') ?? f.get('raw.var')
  let rawSet: string | null = null
  if (rawVarDs) {
    try {
      rawSet = internGeneSet(geneSets, readVar(rawVarDs, 'raw.var'), 'raw')
    } catch { /* no raw genes */ }
  }

  // ---- candidate matrices ----------------------------------------------
  const matrices: MatrixRef[] = []
  const add = (key: string, m: H5, geneSet: string, cells: number, genes: number) => {
    if (!m) return
    try {
      matrices.push({
        key, geneSet, nCells: cells, nGenes: genes,
        kind: classify(sampleMatrix(m)),
        load: () => readMatrix(m, cells, genes, key),
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
      embeddings.push({
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
      })
    })
  } else if (isGroup(obsm)) {
    for (const k of obsm.keys()) {
      const ds = obsm.get(k)
      const s = ds?.shape as number[] | undefined
      if (!s || s.length < 2 || s[1] < 2 || s[0] !== nCells) continue
      embeddings.push({
        key: k, nDims: s[1],
        load: () => {
          const flat = toF32(ds.slice([[0, nCells], [0, 2]]))
          return flat.length === nCells * 2 ? flat : toF32(ds.value).filter((_v, i) => i % s[1] < 2)
        },
      })
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
    geneSets, provenance, notes,
  }
}
