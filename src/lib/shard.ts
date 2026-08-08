// Splitting an atlas into bundles that a browser can actually open.
//
// A 292 k-cell, 736 M-nonzero object cannot be one bundle: the studio holds the
// matrix resident (every marker and DE view scans all genes), at 8 bytes per
// stored value. So the lab splits it along a column the scan already found —
// one bundle per tissue, per class, per age — and each piece opens with the
// full feature set on its own cells.
//
// The split is costed exactly before anything is read, from the row offsets
// alone, so you are told the shape of the result rather than discovering it.

import type { CellMajor } from './matrix.ts'
import type { Column, EmbeddingRef, MatrixRef, Scan } from './scan.ts'

/**
 * What the studio can comfortably hold.
 *
 * The bundle stores 4 bytes of cell index and 4 bytes of value per nonzero, and
 * the studio keeps both resident. 50 M nonzeros is ~400 MB of typed arrays,
 * which a desktop tab manages; past that it becomes unreliable rather than
 * merely slow, and unreliable is worse.
 */
export const COMFORTABLE_NNZ = 50_000_000

/** Nonzeros held in memory at once while reading. Bounds the pass count. */
export const PASS_BUDGET_NNZ = 120_000_000

export interface Shard {
  /** Level name — becomes part of the file name. */
  key: string
  cells: Int32Array
  /** Exact, computed from the row offsets before any values are read. */
  nnz: number
}

export interface SplitPlan {
  column: string
  shards: Shard[]
  /** Shards grouped into one forward pass each, respecting the memory budget. */
  passes: Shard[][]
  totalNnz: number
  /** Shards still too large for the studio even after splitting. */
  oversized: Shard[]
}

/**
 * Cost a split without reading a single value.
 *
 * `nnzPerCell` comes from the CSR row offsets, so the sizes here are exact, not
 * estimates — which matters, because the whole point is to know before spending
 * a minute per pass whether the result will open.
 */
export function planSplit(
  col: Column, nnzPerCell: Int32Array | null, nCells: number,
  budget = PASS_BUDGET_NNZ,
): SplitPlan {
  const buckets: number[][] = col.levels.map(() => [])
  for (let i = 0; i < col.codes.length; i++) {
    const c = col.codes[i]
    if (c >= 0 && c < buckets.length) buckets[c].push(i)
  }

  // Without row offsets (a dense or CSC source) fall back to the mean density,
  // which is honest enough to order the shards but is labelled as an estimate
  // wherever it is shown.
  const mean = nnzPerCell ? 0 : 1
  const shards: Shard[] = buckets
    .map((cells, i) => {
      let nnz = 0
      if (nnzPerCell) for (const c of cells) nnz += nnzPerCell[c]
      else nnz = cells.length * mean
      return { key: col.levels[i], cells: Int32Array.from(cells), nnz }
    })
    .filter(s => s.cells.length > 0)
    .sort((a, b) => b.nnz - a.nnz)

  // Greedy first-fit into passes: each pass holds only what fits in memory,
  // and every pass costs one traversal of the file.
  const passes: Shard[][] = []
  for (const s of shards) {
    const fits = passes.find(p => p.reduce((n, x) => n + x.nnz, 0) + s.nnz <= budget)
    if (fits) fits.push(s)
    else passes.push([s])
  }

  void nCells
  return {
    column: col.name,
    shards,
    passes,
    totalNnz: shards.reduce((n, s) => n + s.nnz, 0),
    oversized: shards.filter(s => s.nnz > COMFORTABLE_NNZ),
  }
}

/** Every categorical column, ranked by how well it splits this object. */
export function splitCandidates(
  columns: Column[], nnzPerCell: Int32Array | null, nCells: number,
): { column: Column; plan: SplitPlan }[] {
  return columns
    .filter(c => c.kind === 'categorical' && c.levels.length > 1 && c.levels.length <= 64)
    .map(column => ({ column, plan: planSplit(column, nnzPerCell, nCells) }))
    // Best split first: fewest shards the studio would choke on, then the
    // smallest worst-case shard, then the fewest files to keep track of.
    .sort((a, b) =>
      a.plan.oversized.length - b.plan.oversized.length ||
      (a.plan.shards[0]?.nnz ?? 0) - (b.plan.shards[0]?.nnz ?? 0) ||
      a.plan.shards.length - b.plan.shards.length)
}

/** Take one column down to a subset of its cells. */
export function subsetColumn(col: Column, cells: Int32Array): Column {
  if (col.kind === 'numeric') {
    const values = new Float32Array(cells.length)
    for (let i = 0; i < cells.length; i++) values[i] = col.values[cells[i]]
    return { ...col, values, codes: new Int32Array(0) }
  }
  // Levels are kept whole rather than compacted to those present: dropping the
  // absent ones would renumber the codes, and a shard whose "cluster 3" means
  // something different from every other shard's is a trap.
  const codes = new Int32Array(cells.length)
  let missing = 0
  for (let i = 0; i < cells.length; i++) {
    codes[i] = col.codes[cells[i]]
    if (codes[i] < 0) missing++
  }
  return { ...col, codes, values: new Float32Array(0), missing }
}

/**
 * A Scan describing one shard, so the existing builder can write it unchanged.
 *
 * The matrix is already in hand — read by the streaming pass — so `load` just
 * hands it over. Everything else is the parent scan with its cells subset.
 */
export function shardScan(
  parent: Scan, shard: Shard, matrices: Map<string, CellMajor>, sourceLabel: string,
): Scan {
  const cells = shard.cells
  return {
    ...parent,
    source: sourceLabel,
    nCells: cells.length,
    columns: parent.columns.map(c => subsetColumn(c, cells)),
    matrices: parent.matrices
      .filter(m => matrices.has(m.key))
      .map((m): MatrixRef => ({
        key: m.key, geneSet: m.geneSet, kind: m.kind,
        nGenes: m.nGenes, nCells: cells.length,
        load: () => matrices.get(m.key)!,
      })),
    embeddings: parent.embeddings.map((e): EmbeddingRef => ({
      key: e.key, nDims: e.nDims,
      load: () => {
        const all = e.load()
        const xy = new Float32Array(cells.length * 2)
        for (let i = 0; i < cells.length; i++) {
          xy[2 * i] = all[2 * cells[i]]
          xy[2 * i + 1] = all[2 * cells[i] + 1]
        }
        return xy
      },
    })),
  }
}

/** A file name that survives a round trip through a file system. */
export const shardFileName = (label: string, key: string): string =>
  `${label}__${key}`.replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 120) + '.zip'
