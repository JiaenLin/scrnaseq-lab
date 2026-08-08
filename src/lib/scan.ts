// What a scan of an object looks like, before any choice has been made.
//
// Both readers produce this and nothing else, so everything downstream — the
// questions, the warnings, the builder — is written once against one shape.

import type { CellMajor, MatrixKind } from './matrix.ts'

/** One column of obs / meta.data, already reduced to codes and levels. */
export interface Column {
  name: string
  kind: 'categorical' | 'numeric'
  /** Level names, in first-appearance or factor order — never sorted. */
  levels: string[]
  /** Level index per cell. Empty for numeric columns. */
  codes: Int32Array
  /** Values per cell. Empty for categorical columns. */
  values: Float32Array
  /** How many cells hold NA. */
  missing: number
}

export interface MatrixRef {
  /** How it is addressed in the file, e.g. `raw.X` or `assays$RNA@counts`. */
  key: string
  nGenes: number
  nCells: number
  /** Which gene list this matrix's columns index. */
  geneSet: string
  kind: MatrixKind
  /** Read it in full. Deferred: most objects hold more matrices than we need. */
  load: () => CellMajor
  /**
   * Nonzeros per cell, without reading the values.
   *
   * For a CSR matrix this is a difference of the row offsets — a few hundred KB
   * for an atlas — so a split can be costed exactly before any of it is read.
   */
  nnzPerCell?: () => Int32Array
  /**
   * Read several disjoint cell subsets in one forward pass.
   *
   * Random access is not an option on a real atlas: the payload is chunked and
   * compressed, so a scattered row read decompresses a whole chunk for ~2500
   * values. Measured on a 736 M-nonzero file, per-row reads run at 20 ms each —
   * 100 minutes for one pass — while a forward scan sustains 11 M nonzeros/s
   * and finishes in about a minute. So everything here is sequential.
   */
  loadGroups?: (groups: Int32Array[], onProgress?: (frac: number) => void) => CellMajor[]
}

export interface EmbeddingRef {
  key: string
  nDims: number
  /** Interleaved x,y — the first two dimensions. */
  load: () => Float32Array
}

/**
 * Everything about a scan that survives postMessage.
 *
 * The reader lives in a worker — it holds an open HDF5 handle and closures over
 * it — but the UI only ever needs the description. So the worker keeps `Scan`
 * and sends this, and every question the interface asks is answered from it.
 */
export interface ScanInfo {
  format: 'h5ad' | 'seurat'
  source: string
  nCells: number
  columns: Column[]
  matrices: { key: string; geneSet: string; nGenes: number; nCells: number; kind: MatrixKind }[]
  embeddings: { key: string; nDims: number }[]
  geneSets: Record<string, string[]>
  provenance: Record<string, string | null>
  notes: string[]
  /** Nonzeros per cell, when the source could supply them cheaply. */
  nnzPerCell?: Int32Array | null
}

export interface Scan {
  format: 'h5ad' | 'seurat'
  /** Filename, used as the default label. */
  source: string
  nCells: number
  columns: Column[]
  matrices: MatrixRef[]
  embeddings: EmbeddingRef[]
  /** Gene names per gene set named by `MatrixRef.geneSet`. */
  geneSets: Record<string, string[]>
  provenance: Record<string, string | null>
  /** Everything the reader had to decide, in the order it decided it. */
  notes: string[]
}

/** The three questions the lab asks. */
export interface Choices {
  cluster: string
  sample: string | null
  condition: string | null
  embedding: string
  label: string
}

// ---------------------------------------------------------------------------
// Turning raw column values into a Column.

/** Distinct values in first-appearance order. Sorting would reorder a time course. */
export function categorical(name: string, values: string[]): Column {
  const index = new Map<string, number>()
  const levels: string[] = []
  const codes = new Int32Array(values.length)
  let missing = 0
  for (let i = 0; i < values.length; i++) {
    const v = values[i]
    if (v === 'NA' || v === '' || v == null) missing++
    let at = index.get(v)
    if (at === undefined) {
      at = levels.length
      levels.push(v)
      index.set(v, at)
    }
    codes[i] = at
  }
  return { name, kind: 'categorical', levels, codes, values: new Float32Array(0), missing }
}

export function fromFactor(name: string, codes: Int32Array, levels: string[]): Column {
  let missing = 0
  for (const c of codes) if (c < 0) missing++
  return { name, kind: 'categorical', levels, codes, values: new Float32Array(0), missing }
}

export function numeric(name: string, values: Float32Array): Column {
  let missing = 0
  for (const v of values) if (!Number.isFinite(v)) missing++
  return { name, kind: 'numeric', levels: [], codes: new Int32Array(0), values, missing }
}

/**
 * A numeric column with few distinct values is a grouping, whatever its dtype.
 *
 * Cluster ids arrive as 0..8 integers more often than as strings, and a column
 * the user cannot pick as a cluster because it happens to be stored as int is a
 * column they will go back to R to convert.
 */
export function maybeDiscrete(col: Column, nCells: number): Column {
  if (col.kind !== 'numeric') return col
  const seen = new Set<number>()
  for (const v of col.values) {
    if (!Number.isFinite(v)) continue
    if (!Number.isInteger(v)) return col
    seen.add(v)
    if (seen.size > 100) return col
  }
  if (seen.size < 2 || seen.size > Math.max(20, nCells / 10)) return col
  const order = [...seen].sort((a, b) => a - b)
  const at = new Map(order.map((v, i) => [v, i]))
  const codes = new Int32Array(col.values.length)
  for (let i = 0; i < col.values.length; i++) codes[i] = at.get(col.values[i]) ?? -1
  return {
    name: col.name, kind: 'categorical', levels: order.map(String), codes,
    values: new Float32Array(0), missing: col.missing,
  }
}

// ---------------------------------------------------------------------------
// Ranking the candidates. These names are conventions, not a standard, so a
// match only ever moves a column up the list — it never picks one silently.

const CLUSTER_NAMES = [
  'cell_type', 'celltype', 'cell.type', 'cell_types', 'seurat_annotations',
  'annotation', 'annotations', 'cell_ontology_class', 'louvain', 'leiden',
  'seurat_clusters', 'clusters', 'cluster', 'ident',
]
const SAMPLE_NAMES = [
  'sample', 'sample_id', 'sampleid', 'orig.ident', 'donor', 'donor_id',
  'patient', 'subject', 'animal', 'mouse', 'library', 'batch', 'replicate',
]
const CONDITION_NAMES = [
  'condition', 'treatment', 'group', 'genotype', 'stim', 'status', 'disease',
  'timepoint', 'time_point', 'time', 'day', 'age', 'tissue', 'region', 'sex',
]

const rank = (name: string, list: string[]): number => {
  const n = name.toLowerCase()
  const exact = list.indexOf(n)
  if (exact >= 0) return exact
  const partial = list.findIndex(k => n.includes(k))
  return partial >= 0 ? 100 + partial : 1000
}

export type Role = 'cluster' | 'sample' | 'condition'

/** Is this column usable in this role at all, and if not, why not. */
export function unusable(col: Column, role: Role, nCells: number): string | null {
  if (col.kind !== 'categorical') return 'not a grouping — it holds a measurement per cell'
  const k = col.levels.length
  if (k < 2) return `only one value (${col.levels[0] ?? 'none'})`
  // A column with roughly one value per cell is a barcode or a UMI, and it is
  // not a grouping in any role. The threshold matches what the message claims:
  // an earlier version said "nearly one per cell" while rejecting anything past
  // a third, which threw out real annotations on small objects.
  if (nCells >= 4 && k >= nCells * 0.9) {
    return `${k} values across ${nCells} cells — nearly one per cell`
  }
  if (role === 'cluster' && k > 300) {
    return `${k} distinct values — that is an identifier, not an annotation`
  }
  if (role === 'sample' && k > 200) return `${k} distinct values — too many to be samples`
  if (role === 'condition' && k > 50) return `${k} distinct values — too many to be groups`
  return null
}

/**
 * Do two columns cut the cells the same way, whatever they call the pieces?
 *
 * `seurat_clusters` holding 0..8, `RNA_snn_res.0.5` holding 0..8 and
 * `seurat_annotations` holding "B", "NK" are routinely the same partition
 * three times over. Comparing codes directly would miss it — the level order
 * differs — so this checks for a bijection between the two labellings.
 */
export function samePartition(a: Column, b: Column): boolean {
  if (a.codes.length !== b.codes.length || a.levels.length !== b.levels.length) return false
  const fwd = new Map<number, number>()
  const back = new Map<number, number>()
  for (let i = 0; i < a.codes.length; i++) {
    const x = a.codes[i], y = b.codes[i]
    const seen = fwd.get(x)
    if (seen === undefined) fwd.set(x, y)
    else if (seen !== y) return false
    const rev = back.get(y)
    if (rev === undefined) back.set(y, x)
    else if (rev !== x) return false
  }
  return true
}

/** For each column, the earlier column it duplicates — so the UI can say so. */
export function duplicateOf(columns: Column[]): Map<string, string> {
  const out = new Map<string, string>()
  const cats = columns.filter(c => c.kind === 'categorical')
  for (let i = 0; i < cats.length; i++) {
    for (let j = 0; j < i; j++) {
      if (out.has(cats[j].name)) continue
      if (samePartition(cats[i], cats[j])) { out.set(cats[i].name, cats[j].name); break }
    }
  }
  return out
}

/**
 * Candidates for a role, best first.
 *
 * `exclude` drops columns already spoken for by another role. A column cannot
 * be both the thing cells are split by and the thing they are compared across,
 * and offering it in both places invites a bundle whose condition and cluster
 * are the same partition — which produces comparisons with nothing to compare.
 */
export function candidates(scan: ScanInfo, role: Role, exclude: (string | null)[] = []): Column[] {
  const names = role === 'cluster' ? CLUSTER_NAMES
    : role === 'sample' ? SAMPLE_NAMES : CONDITION_NAMES
  const taken = scan.columns.filter(c => exclude.includes(c.name))
  return scan.columns
    .filter(c => !unusable(c, role, scan.nCells))
    .filter(c => !taken.some(t => t.name === c.name || samePartition(t, c)))
    .map(c => ({ c, r: rank(c.name, names), k: c.levels.length, op: opaque(c) }))
    // Opaque codings sink in every role. After that the roles want opposite
    // things: a cell annotation is the finest *named* grouping in an object,
    // while a sample or a condition is a coarse one — so cluster prefers more
    // levels and the others fewer. Breaking cluster ties by fewest levels made
    // a 5-value tissue-ontology column outrank 133 named cell types on a real
    // atlas, and the bundle came out labelled UBERON:0000033.
    .sort((a, b) => a.r - b.r || Number(a.op) - Number(b.op) ||
      (role === 'cluster' ? b.k - a.k : a.k - b.k))
    .map(x => x.c)
}

/**
 * Are these levels identifiers rather than names?
 *
 * `UBERON:0000955`, `EFO:0009899` and bare integers all group cells correctly
 * and label a figure uselessly. They stay selectable — some objects really are
 * annotated this way — but they should never win by default over something a
 * reader can understand.
 */
function opaque(c: Column): boolean {
  const named = c.levels.filter(l => /[A-Za-z]/.test(l) && !/^[A-Za-z]+[:_]\d+$/.test(l))
  return named.length * 2 < c.levels.length
}

/**
 * The cluster column we would pick if asked to guess.
 *
 * Prefers a named annotation over a numbered one: `seurat_clusters` holding
 * 0..8 and `seurat_annotations` holding "B", "NK" describe the same grouping,
 * and only one of them makes a readable figure.
 */
export function guessCluster(scan: ScanInfo): string | null {
  const cands = candidates(scan, 'cluster')
  if (!cands.length) return null
  const named = cands.find(c =>
    rank(c.name, CLUSTER_NAMES) < 1000 && c.levels.some(l => !/^-?\d+$/.test(l)))
  return (named ?? cands[0]).name
}

export function guessRole(scan: ScanInfo, role: 'sample' | 'condition'): string | null {
  const names = role === 'sample' ? SAMPLE_NAMES : CONDITION_NAMES
  const hit = candidates(scan, role).find(c => rank(c.name, names) < 1000)
  return hit?.name ?? null
}

const EMBED_ORDER = ['umap', 'tsne', 'draw_graph', 'pca']

/** UMAP, then t-SNE, then a force layout, then PCA. */
export function guessEmbedding(scan: ScanInfo): string | null {
  const usable = scan.embeddings.filter(e => e.nDims >= 2)
  if (!usable.length) return null
  for (const want of EMBED_ORDER) {
    const hit = usable.find(e => e.key.toLowerCase().includes(want))
    if (hit) return hit.key
  }
  return usable[0].key
}

export const column = (scan: ScanInfo, name: string | null): Column | undefined =>
  name ? scan.columns.find(c => c.name === name) : undefined

/**
 * Register a gene list, reusing an identical one already registered.
 *
 * Two matrices share a gene set when their gene names match, not when they sit
 * in the same slot. Keying by slot name made `@counts` and `@data` — the same
 * 13 714 features either way — look like different gene sets, and pseudobulk
 * was then disabled to avoid a mismatch that did not exist.
 */
export function internGeneSet(
  sets: Record<string, string[]>, genes: string[], preferred: string,
): string {
  for (const [key, have] of Object.entries(sets)) {
    if (have.length !== genes.length) continue
    let same = true
    for (let i = 0; i < have.length; i++) {
      if (have[i] !== genes[i]) { same = false; break }
    }
    if (same) return key
  }
  let key = preferred
  for (let i = 2; key in sets; i++) key = `${preferred}#${i}`
  sets[key] = genes
  return key
}

/**
 * Does each sample sit in exactly one condition?
 *
 * If not, the pairing is not a design — it is two independent labels, and the
 * bundle's one-condition-per-sample model would keep whichever came first.
 */
export function crossCheck(
  sample: Column | undefined, condition: Column | undefined,
): string | null {
  if (!sample || !condition) return null
  const seen = new Map<number, number>()
  for (let i = 0; i < sample.codes.length; i++) {
    const s = sample.codes[i], c = condition.codes[i]
    const prev = seen.get(s)
    if (prev === undefined) seen.set(s, c)
    else if (prev !== c) {
      return `sample "${sample.levels[s]}" appears in both "${condition.levels[prev]}" and "${condition.levels[c]}". Each sample must belong to one group, so pick a different pair — or the bundle will record only the first.`
    }
  }
  return null
}

/** Replicates per condition — what decides whether pseudobulk is offered. */
export function replicatesPerCondition(
  sample: Column | undefined, condition: Column | undefined,
): Record<string, number> {
  if (!condition) return {}
  const out: Record<string, Set<number>> = {}
  for (const l of condition.levels) out[l] = new Set()
  for (let i = 0; i < condition.codes.length; i++) {
    const c = condition.levels[condition.codes[i]]
    if (c !== undefined) out[c]?.add(sample ? sample.codes[i] : 0)
  }
  return Object.fromEntries(Object.entries(out).map(([k, v]) => [k, v.size]))
}
