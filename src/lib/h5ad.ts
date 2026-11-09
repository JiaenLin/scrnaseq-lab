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
  type Column, type EmbeddingRef, type GeneAlias, type MatrixRef, type RunRecord, type RunSection, type Scan,
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
 * A string column in anndata's nullable encoding: a group of `values` + `mask`.
 *
 * This is what anndata (0.11 on) writes for pandas' own string dtype, and on
 * pandas 3 that dtype is the default — so an object written with default
 * settings stores its obs and var names this way, and every string column it
 * did not turn into a categorical. The family's writers opt out of it; anything
 * written with the defaults does not, and the reader used to fail on the first
 * one with a property read on undefined. A masked entry is missing: its slot in
 * `values` is a placeholder, never the value.
 *
 * Null for anything else, so callers fall through to the encodings they knew.
 */
function nullableStrings(o: H5): { values: string[]; mask: Uint8Array } | null {
  if (!isGroup(o) || text(attr(o, 'encoding-type') ?? '') !== 'nullable-string-array') return null
  const v = o.get('values')
  if (!isDataset(v)) return null
  const values = toStrings(v.value)
  const mask = new Uint8Array(values.length)
  const m = o.get('mask')
  if (isDataset(m)) {
    const raw = m.value as ArrayLike<number | bigint | boolean>
    for (let i = 0; i < mask.length; i++) mask[i] = raw[i] ? 1 : 0
  }
  return { values, mask }
}

/** Masked entries as empty strings — where every caller treats an empty name as no name. */
const unmasked = (s: { values: string[]; mask: Uint8Array }): string[] =>
  s.values.map((v, i) => (s.mask[i] ? '' : v))

/** A nullable string column as a grouping: masked cells are NA, never a level of their own. */
function maskedColumn(name: string, s: { values: string[]; mask: Uint8Array }): Column {
  const index = new Map<string, number>()
  const levels: string[] = []
  const codes = new Int32Array(s.values.length)
  for (let i = 0; i < s.values.length; i++) {
    if (s.mask[i]) { codes[i] = -1; continue }
    let at = index.get(s.values[i])
    if (at === undefined) {
      at = levels.length
      levels.push(s.values[i])
      index.set(s.values[i], at)
    }
    codes[i] = at
  }
  return fromFactor(name, codes, levels)
}

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
        const strings = nullableStrings(o)
        if (strings) {
          if (strings.values.length !== nCells) {
            note(`obs/${key} has ${strings.values.length} values for ${nCells} cells; skipped`)
          } else {
            out.push(maskedColumn(key, strings))
          }
          continue
        }
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
    const strings = nullableStrings(ds)
    if (strings) return unmasked(strings)
    if (isDataset(ds)) return toStrings(ds.value)
    if (ds) {
      throw new H5adError(`${label}/${name} is stored as "${text(attr(ds, 'encoding-type') ?? 'an unnamed group')}", which this reader cannot read as gene names`)
    }
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
      const strings = nullableStrings(o)
      if (strings) {
        if (strings.values.length === nGenes) out.push({ name: key, values: unmasked(strings) })
        continue
      }
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


/**
 * How much of the gene list a CSR matrix ever touches.
 *
 * Read the same way values are sampled — a few contiguous blocks spread across
 * the index vector — because a single block covers only a handful of cells and
 * would understate a matrix whose detected genes vary from cell to cell. That
 * is the whole distinction being measured: a real expression matrix touches
 * different genes in different cells and adds up to most of the list, while a
 * model fitted on selected features touches the same few in every cell however
 * far you read.
 */
function geneCoverage(m: H5, nGenes: number): number {
  const idx = m.get('indices')
  const n = (idx.shape as number[])[0]
  if (!n || !nGenes) return 1
  const blocks = 24
  const per = Math.max(1, Math.min(200_000, Math.floor(n / blocks)))
  const seen = new Set<number>()
  for (let b = 0; b < blocks; b++) {
    const at = Math.max(0, Math.min(n - per, Math.floor((n / blocks) * b)))
    const s = idx.slice([[at, at + per]]) as ArrayLike<number | bigint>
    for (let i = 0; i < s.length; i++) seen.add(Number(s[i]))
    if (seen.size >= nGenes) break
  }
  return seen.size / nGenes
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

// ---------------------------------------------------------------------------
// The run's own record.
//
// An object says in uns what made it, and every tool says it in its own place:
// the integration stage under uns/<tool>/provenance, the annotator as
// uns/<column>_provenance, the QC stage as a group whose schema names
// provenance, and the harness's run record (when there is one) as a JSON
// string in uns/harness_run. The reader cannot know every tool, so it names
// none: a record is what calls itself provenance, by its key or by its own
// schema, at the top of uns or one level down. Deeper than that is inside some
// other record, not a record of the run. What it finds is copied, not
// summarised — a summary would decide which of a tool's fields matter.

/** Where the reader looks for the run's record. Written into the bundle beside what it found. */
export const RUN_SEARCHED = 'uns/harness_run; every entry of uns, and every child of a dict directly under uns, named provenance or *_provenance; every dict directly under uns whose own schema names provenance'

/** Arrays longer than this are named in a record rather than read into it. */
const RUN_ARRAY_CAP = 10_000

/**
 * The most one record may take in meta.json, in bytes as written; past it the
 * record is named in place with its size, not copied.
 *
 * The run section is written into every part of a split object, so a record
 * is paid for once per part. A tool's provenance is a few hundred bytes to a
 * few tens of thousands (an annotator's per-label counts on an atlas); a
 * quarter of a megabyte holds those with room to spare, and stops a record
 * that is really a result (a log, a table) from being copied forty times.
 *
 * Bytes as written, not compact JSON: meta.json is indented, one space a
 * level, and an array of numbers written one per line takes several times
 * its compact size.
 */
export const RUN_RECORD_BUDGET = 262_144

/**
 * The most all the written records may take together, in bytes as written.
 * Records that do not fit are left out, and one entry, `leftOut`, counts them.
 */
export const RUN_SECTION_BUDGET = 1_048_576

/**
 * Groups deeper than this below their record's own entry are named in place.
 * Reading a chain of groups costs more at every level (800 levels took a
 * minute and a half, 2,000 did not finish), and a run record is not that deep.
 */
export const RUN_MAX_DEPTH = 16

/**
 * The most groups and datasets read for one record; past it the rest is named
 * in place. Depth alone does not bound a read: hard links make a few groups
 * into exponentially many paths (17 groups of two links each are 2^17 paths in
 * a 265 KB file, which read for minutes and never finished).
 */
export const RUN_MAX_NODES = 10_000

/** Records read in full before the section is filled; beyond this, a record is read again if it is kept. */
const RUN_HOLD = 8 * RUN_SECTION_BUDGET

/**
 * One record's walk: the HDF5 objects reached so far, by identity, and how
 * many groups and datasets have been read. An object reached again — a hard
 * link back to an ancestor, or a second link to one already read — is named
 * in place, never read again.
 */
interface Walk {
  /** The record's own path, for the notes. */
  record: string
  seen: Map<string, string>
  visited: number
  relinked: Set<string>
  crowded: Set<string>
}

/**
 * An HDF5 object's identity: the object reference h5wasm makes for it, which
 * is its address in the file, whatever path reached it. Null where the
 * reader cannot make one, and then that object is simply not tracked.
 */
function identity(o: H5): string | null {
  try {
    const bytes = o.create_reference?.()?.ref_data as ArrayLike<number> | undefined
    return bytes ? Array.from(bytes).join('.') : null
  } catch {
    return null
  }
}

/**
 * HDF5 type classes JSON has no form for, named in place rather than copied:
 * a reference's bytes mean nothing outside its file, and opaque bytes, a
 * bitfield or a time have no agreed reading.
 */
const UNHELD: Record<number, string> = {
  2: 'HDF5 time values', 4: 'an HDF5 bitfield', 5: 'opaque HDF5 bytes', 7: 'HDF5 references',
}

/**
 * A dict in anndata's layout, or a plain group from a writer that names no encoding.
 *
 * (Tailwind reads this directory for class names, so a variable here is never
 * called `t` and negated: the app has a `t` class, and the scan took the
 * negation for its important variant and shipped it as CSS.)
 */
const isDict = (g: H5): boolean => {
  const encoding = text(attr(g, 'encoding-type') ?? '')
  return isGroup(g) && (encoding === '' || encoding === 'dict')
}

/**
 * One stored value as JSON can hold it: int64 as a number (its exact digits as
 * a string past 2**53), NaN as "NaN", and a list — a compound row, an array
 * member, a ragged row — value by value. A row of a structured array arrives
 * as a list holding BigInts, and copying it as it came crashed the whole
 * conversion at JSON.stringify.
 */
function jsonValue(v: unknown): unknown {
  if (typeof v === 'bigint') {
    return v >= BigInt(Number.MIN_SAFE_INTEGER) && v <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(v) : v.toString()
  }
  if (typeof v === 'number') return Number.isFinite(v) ? v : String(v)
  if (Array.isArray(v) || ArrayBuffer.isView(v)) return Array.from(v as ArrayLike<unknown>, jsonValue)
  if (v !== null && typeof v === 'object') return '(not copied: a value JSON has no form for)'
  return v ?? null
}

/** A flat list back into the dataset's shape, row-major. */
function nest(flat: unknown[], shape: number[]): unknown[] {
  if (shape.length <= 1) return flat
  const step = shape.slice(1).reduce((a, b) => a * b, 1)
  const out: unknown[] = []
  for (let i = 0; i < flat.length; i += step) out.push(nest(flat.slice(i, i + step), shape.slice(1)))
  return out
}

/**
 * A dataset as JSON: a scalar as itself, an array as a list (nested by shape),
 * a structured array as one object per row under its field names.
 *
 * What cannot be copied is named in place, so one odd leaf never costs the
 * rest of its record.
 */
function plainDataset(d: H5, path: string): unknown {
  const info = d.metadata ?? {}
  const unheld = UNHELD[info.type as number]
  if (unheld) return `(not copied: ${path} holds ${unheld})`
  if (info.type === 1 && info.size === 2) return `(not copied: ${path} is float16, which this reader cannot read)`
  const shape = (d.shape as number[] | undefined) ?? []
  const size = shape.reduce((a, b) => a * b, 1)
  if (size > RUN_ARRAY_CAP) return `(not copied: ${path} holds ${size.toLocaleString('en-US')} values)`
  let v: unknown
  try {
    v = d.value
  } catch (e) {
    return `(not copied: ${path} could not be read: ${(e as Error).message})`
  }
  if (v == null) return `(not copied: ${path} is of a type this reader cannot read)`
  // h5py stores a numpy bool as an HDF5 enum of FALSE and TRUE, and h5wasm
  // hands back the integer. It was a boolean when it was written; any other
  // enum value is its member's name.
  const members = info.enum_type?.members as Record<string, number> | undefined
  const bool = !!members && Object.keys(members).length === 2 && 'TRUE' in members && 'FALSE' in members
  const fields = info.compound_type?.members as { name: string }[] | undefined
  const one = (x: unknown): unknown => {
    if (bool) return Number(x) === members!.TRUE
    if (members) return Object.keys(members).find(k => members[k] === Number(x)) ?? jsonValue(x)
    if (fields && Array.isArray(x)) return Object.fromEntries(fields.map((f, i) => [f.name, jsonValue(x[i])]))
    return jsonValue(x)
  }
  if (!shape.length) {
    // A scalar comes back as itself, or wrapped in a one-element list; a
    // compound scalar is its row of fields, possibly wrapped the same way.
    if (fields) return one(Array.isArray(v) && Array.isArray(v[0]) ? v[0] : v)
    return one(typeof v === 'object' && 'length' in v ? (v as ArrayLike<unknown>)[0] : v)
  }
  return nest(Array.from(v as ArrayLike<unknown>, one), shape)
}

/**
 * An uns entry as JSON: a dict as an object, in the order the file lists it.
 *
 * Three things bound the read, each naming what it stops in place: `depth`
 * counts groups below the record's own entry (past RUN_MAX_DEPTH a group is
 * not read); the walk counts the groups and datasets read (past
 * RUN_MAX_NODES the rest are not); and an object already reached by another
 * link is not read again.
 */
function plain(o: H5, path: string, walk: Walk, depth = 0): unknown {
  if (!isDataset(o) && !isGroup(o)) return null
  if (walk.visited >= RUN_MAX_NODES) {
    walk.crowded.add(walk.record)
    return `(not copied: ${path} lies past the ${RUN_MAX_NODES.toLocaleString('en-US')} groups and datasets read for one record)`
  }
  walk.visited++
  const id = identity(o)
  if (id !== null) {
    const first = walk.seen.get(id)
    if (first !== undefined) {
      walk.relinked.add(walk.record)
      return `(not copied: ${path} is ${first} again, reached by another link)`
    }
    walk.seen.set(id, path)
  }
  if (isDataset(o)) return plainDataset(o, path)
  if (depth > RUN_MAX_DEPTH) return `(not copied: ${path} is nested more than ${RUN_MAX_DEPTH} groups below its record)`
  const strings = nullableStrings(o)
  if (strings) {
    return strings.values.length > RUN_ARRAY_CAP
      ? `(not copied: ${path} holds ${strings.values.length.toLocaleString('en-US')} values)`
      : strings.values.map((s, i) => (strings.mask[i] ? null : s))
  }
  if (!isDict(o)) return `(not copied: ${path} is a ${text(attr(o, 'encoding-type'))})`
  const out: Record<string, unknown> = {}
  const keys = o.keys() as string[]
  for (let i = 0; i < keys.length; i++) {
    if (walk.visited >= RUN_MAX_NODES) {
      // One sentence for the rest of this group, not one per entry.
      walk.crowded.add(walk.record)
      out['(not copied)'] = `${(keys.length - i).toLocaleString('en-US')} more entries of ${path}, past the ${RUN_MAX_NODES.toLocaleString('en-US')} groups and datasets read for one record`
      break
    }
    out[keys[i]] = plain(o.get(keys[i]), `${path}/${keys[i]}`, walk, depth + 1)
  }
  return out
}

/**
 * The bytes a string takes in meta.json: UTF-8, quoted, escaped as
 * JSON.stringify escapes it (a lone surrogate as \uXXXX).
 */
function stringBytes(s: string): number {
  let b = 2
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i)
    if (c === 0x22 || c === 0x5c) b += 2
    else if (c < 0x20) b += c === 8 || c === 9 || c === 10 || c === 12 || c === 13 ? 2 : 6
    else if (c < 0x80) b += 1
    else if (c < 0x800) b += 2
    else if (c >= 0xd800 && c <= 0xdbff) {
      const d = s.charCodeAt(i + 1)
      if (d >= 0xdc00 && d <= 0xdfff) { b += 4; i++ } else b += 6
    } else if (c >= 0xdc00 && c <= 0xdfff) b += 6
    else b += 3
  }
  return b
}

/**
 * The bytes a value takes in meta.json, which the builder writes as
 * JSON.stringify(meta, null, 1): at nesting `depth`, each member on its own
 * line indented depth + 1 spaces, the closing bracket on a line indented
 * `depth`. Counted, not written out, so a 50 MB record is measured without a
 * second 50 MB string. A record sits at depth 4 (meta, run, records, entry).
 */
function writtenSize(v: unknown, depth: number): number {
  if (Array.isArray(v)) {
    if (!v.length) return 2
    let s = v.length * (depth + 3) + depth + 2
    for (const x of v) s += writtenSize(x, depth + 1)
    return s
  }
  if (v !== null && typeof v === 'object') {
    const keys = Object.keys(v).filter(k => (v as Record<string, unknown>)[k] !== undefined)
    if (!keys.length) return 2
    let s = keys.length * (depth + 5) + depth + 2
    for (const k of keys) s += stringBytes(k) + writtenSize((v as Record<string, unknown>)[k], depth + 1)
    return s
  }
  if (typeof v === 'string') return stringBytes(v)
  // Numbers, booleans and null are ASCII; anything else JSON cannot hold
  // throws here, which fails its record and never the conversion.
  const json = JSON.stringify(v ?? null)
  if (json === undefined) throw new Error(`a ${typeof v} has no JSON form`)
  return json.length
}

/**
 * What one entry of run.records adds to meta.json, given what its record
 * takes at depth 4: the path, and 40 bytes of indentation and punctuation.
 */
const entryBytes = (path: string, recordBytes: number) => 40 + stringBytes(path) + recordBytes

/** Paths for a sentence: the first ten, then a count. */
function some(paths: string[]): string {
  const shown = paths.slice(0, 10).join(', ')
  return paths.length > 10 ? `${shown} and ${(paths.length - 10).toLocaleString('en-US')} more` : shown
}

/** Does this group name provenance in its own `schema` string? */
function schemaSaysProvenance(o: H5): boolean {
  if (!isDict(o)) return false
  try {
    const s = o.get('schema')
    return isDataset(s) && /provenance/i.test(text(s.value))
  } catch {
    return false
  }
}

/**
 * Every record of the run the object carries.
 *
 * Found in uns order, then read, measured as meta.json will hold it, and
 * written in uns order. What is written is decided kind by kind: the harness's
 * run record first, then each `provenance`, then each `*_provenance`, then the
 * records known by their schema. Within a kind the records are taken in uns
 * order, first fit: one that does not fit in what is left of
 * RUN_SECTION_BUDGET is left out, and the next is still tried. So a record of
 * an earlier kind is never pushed out by a later kind; within one kind,
 * earlier records in uns can still fill the section before a tool's own
 * record, and then that record is left out and named first, because the
 * records left out are named largest first. A record over RUN_RECORD_BUDGET is
 * written as a sentence naming its size. One entry, `leftOut`, and one note
 * cover every record left out.
 *
 * Nothing found here can stop the conversion. An uns entry that cannot even be
 * inspected, a record that cannot be read, and a value JSON cannot hold are
 * each said in one note, never thrown.
 */
function readRun(uns: H5, note: (s: string) => void): RunSection {
  const records: RunRecord[] = []
  if (!isGroup(uns)) return { searched: RUN_SEARCHED, records }
  const named = (k: string) => k === 'provenance' || k.endsWith('_provenance')
  const kindOf = (k: string) => (k === 'provenance' ? 1 : 2)
  const n = (x: number) => x.toLocaleString('en-US')
  const relinked = new Set<string>()
  const crowded = new Set<string>()
  const walkOf = (record: string): Walk => ({ record, seen: new Map(), visited: 0, relinked, crowded })

  // ---- found, in uns order; each entry inspected on its own -------------
  type Found = { path: string; kind: number; read: (say: (s: string) => void) => unknown }
  const found: Found[] = []
  const uninspected: string[] = []
  let why = ''
  for (const k of uns.keys() as string[]) {
    const path = `uns/${k}`
    try {
      const o = uns.get(k)
      if (!o) continue
      if (k === 'harness_run') {
        // The harness's record is one JSON string. Parsed, so that it reads
        // as the structure it is; kept as its text when it does not parse.
        found.push({ path, kind: 0, read: say => {
          const v = plain(o, path, walkOf(path))
          if (typeof v !== 'string') return v
          try {
            return JSON.parse(v)
          } catch {
            say(`${path} is not JSON, so it travels in run as its text`)
            return v
          }
        } })
      } else if (named(k)) {
        found.push({ path, kind: kindOf(k), read: () => plain(o, path, walkOf(path)) })
      } else if (schemaSaysProvenance(o)) {
        found.push({ path, kind: 3, read: () => plain(o, path, walkOf(path)) })
      } else if (isDict(o)) {
        for (const c of o.keys() as string[]) {
          const at = `${path}/${c}`
          if (named(c)) found.push({ path: at, kind: kindOf(c), read: () => plain(o.get(c), at, walkOf(at)) })
        }
      }
    } catch (e) {
      // An attribute h5wasm cannot read (float16), a broken link: this entry
      // cannot be classified, and the conversion goes on without it.
      uninspected.push(path)
      why ||= (e as Error).message
    }
  }

  // ---- read and measured ---------------------------------------------------
  type Read = { at: number; size: number; record?: unknown; failed?: string }
  const read: Read[] = []
  let held = 0
  found.forEach((f, at) => {
    try {
      const record = f.read(note) ?? null
      const size = writtenSize(record, 4)
      // Held for the filling below while it stays small; a record not held is
      // read again if it is written.
      const hold = size <= RUN_RECORD_BUDGET && held + size <= RUN_HOLD
      if (hold) held += size
      read.push({ at, size, ...(hold ? { record } : {}) })
    } catch (e) {
      read.push({ at, size: 0, failed: (e as Error).message })
    }
  })

  // ---- filled: kind by kind, first fit in uns order --------------------------
  const order = [...read].sort((a, b) => found[a.at].kind - found[b.at].kind || a.at - b.at)
  const written = new Map<number, unknown>()
  const failed: string[] = []
  const tooBig: string[] = []
  const leftOut: Read[] = []
  let spent = 0
  for (const r of order) {
    const path = found[r.at].path
    const over = r.failed === undefined && r.size > RUN_RECORD_BUDGET
    const sentence = `(not copied: ${n(r.size)} bytes as written, over the ${n(RUN_RECORD_BUDGET)}-byte budget for one record)`
    const cost = entryBytes(path, r.failed !== undefined ? 4 : over ? stringBytes(sentence) : r.size)
    if (spent + cost > RUN_SECTION_BUDGET) { leftOut.push(r); continue }
    let record: unknown = null
    if (over) record = sentence
    else if (r.failed === undefined) {
      try {
        record = 'record' in r ? r.record : (found[r.at].read(() => {}) ?? null)
      } catch (e) {
        r.failed = (e as Error).message
      }
    }
    spent += cost
    written.set(r.at, record)
    if (r.failed !== undefined) failed.push(path)
    else if (over) tooBig.push(path)
  }
  for (const r of read) if (written.has(r.at)) records.push({ path: found[r.at].path, record: written.get(r.at) })
  // Largest first, so that a tool's own record, which is rarely the smallest,
  // is named among the first few whenever it is left out.
  const outPaths = [...leftOut].sort((a, b) => b.size - a.size || a.at - b.at).map(r => found[r.at].path)

  // ---- said, once each -------------------------------------------------------
  if (uninspected.length) {
    note(`${some(uninspected)} could not be inspected (${why}), so the run section cannot tell whether ${uninspected.length > 1 ? 'they hold records' : 'it holds a record'}`)
  }
  if (relinked.size) {
    note(`${some([...relinked])} ${relinked.size > 1 ? 'reach' : 'reaches'} an HDF5 object by more than one link; the object is copied once, and each other link to it is named in place`)
  }
  if (crowded.size) {
    note(`${some([...crowded])} ${crowded.size > 1 ? 'hold' : 'holds'} more than ${n(RUN_MAX_NODES)} groups and datasets; past that, the rest is named in place, unread`)
  }
  if (failed.length) {
    note(`${some(failed)} could not be read; run lists ${failed.length > 1 ? 'them' : 'it'} with no content`)
  }
  if (tooBig.length) {
    note(`${some(tooBig)} ${tooBig.length > 1 ? 'are' : 'is'} over the ${n(RUN_RECORD_BUDGET)}-byte budget for one record as meta.json writes it, so run names ${tooBig.length > 1 ? 'each' : 'it'} with its size instead of copying it into every part`)
  }
  if (leftOut.length) {
    note(`${n(leftOut.length)} record${leftOut.length > 1 ? 's were' : ' was'} left out of run, largest first: ${some(outPaths)}. The section had reached its ${n(RUN_SECTION_BUDGET)}-byte budget, filled kind by kind (the harness's run record, then each provenance, then *_provenance, then records known by their schema), each kind in uns order`)
  }
  return {
    searched: RUN_SEARCHED, records,
    ...(leftOut.length ? { leftOut: { count: leftOut.length, paths: outPaths.slice(0, 10) } } : {}),
  }
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
          // Only CSR: in a CSC matrix `indices` are cell ids, and counting
          // distinct cells would answer a different question entirely.
          coverage: () => geneCoverage(m, genes),
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
  const run = readRun(uns, note)

  return {
    format: 'h5ad', source: filename, nCells, columns, matrices, embeddings,
    geneSets, geneAliases, provenance, run, notes,
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
