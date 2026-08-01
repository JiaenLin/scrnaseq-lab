// Sparse matrices, and deciding what a matrix actually contains.

/**
 * Cell-major sparse: one run of (gene, value) pairs per cell.
 *
 * Both input formats already look like this — an .h5ad CSR is cells × genes,
 * and a Seurat dgCMatrix is genes × cells stored by column, which is the same
 * runs in the same order. So this is the shape data arrives in, and the single
 * transpose to gene-major happens once, at the end.
 */
export interface CellMajor {
  nCells: number
  nGenes: number
  /** nCells + 1 */
  indptr: Int32Array
  /** gene index of each stored value */
  indices: Int32Array
  data: Float32Array
}

/** Gene-major CSC over cells — the orientation the bundle stores. */
export interface GeneMajor {
  indptr: Int32Array
  indices: Int32Array
  data: Float32Array
}

export type MatrixKind = 'scaled' | 'counts' | 'log-counts' | 'lognorm' | 'linear'

/** numpy's allclose, so the ported thresholds mean what they meant there. */
const allClose = (v: ArrayLike<number>, round: (x: number) => number,
                  atol = 1e-8, rtol = 1e-5): boolean => {
  for (let i = 0; i < v.length; i++) {
    const b = round(v[i])
    if (!(Math.abs(v[i] - b) <= atol + rtol * Math.abs(b))) return false
  }
  return true
}

/**
 * What kind of matrix is this?
 *
 * Slot names lie. An .h5ad routinely keeps scaled values in X, log1p(counts)
 * in .raw and nothing in layers; a Seurat object keeps scaled values in
 * `scale.data` but not always. So classify by the numbers, never the name.
 *
 *   scaled      has negatives; unusable for any expression view
 *   counts      integer
 *   log-counts  expm1 is integer — logged, but never depth-normalized
 *   lognorm     small positive non-integers; already log1p(CP10K) or similar
 *   linear      large positive non-integers; normalized but not logged
 *
 * The sample is strided across the whole vector rather than taken from the
 * front. A matrix that is scaled only in its later genes is not a thing that
 * happens by accident, but reading the first 20 000 values of a 200-million
 * value array and calling it representative is a guess, and striding costs
 * nothing.
 */
export function classify(sample: ArrayLike<number>): MatrixKind {
  let min = Infinity, max = -Infinity
  const nz: number[] = []
  for (let i = 0; i < sample.length; i++) {
    const v = sample[i]
    if (v === 0 || !Number.isFinite(v)) continue
    nz.push(v)
    if (v < min) min = v
    if (v > max) max = v
  }
  if (nz.length === 0) return 'lognorm'
  if (min < 0) return 'scaled'
  if (allClose(nz, Math.round)) return 'counts'
  if (allClose(nz.map(Math.expm1), Math.round, 1e-3)) return 'log-counts'
  return max < 50 ? 'lognorm' : 'linear'
}

/** Take at most `want` values spread across a vector of length `n`. */
export function stridedSample(
  n: number, want: number, read: (i: number) => number,
): Float64Array {
  const take = Math.min(n, want)
  const out = new Float64Array(take)
  const step = Math.max(1, Math.floor(n / take))
  for (let i = 0; i < take; i++) out[i] = read(Math.min(n - 1, i * step))
  return out
}

/** log1p(CP10K), in place — the transform every view in the studio assumes. */
export function lognormalize(m: CellMajor): CellMajor {
  const data = new Float32Array(m.data.length)
  for (let c = 0; c < m.nCells; c++) {
    const a = m.indptr[c], b = m.indptr[c + 1]
    let total = 0
    for (let k = a; k < b; k++) total += m.data[k]
    const f = total > 0 ? 1e4 / total : 0
    for (let k = a; k < b; k++) data[k] = Math.log1p(m.data[k] * f)
  }
  return { ...m, data }
}

/** expm1 and round — recovering counts that were only hiding behind a log. */
export function expm1Round(m: CellMajor): CellMajor {
  const data = new Float32Array(m.data.length)
  for (let i = 0; i < m.data.length; i++) data[i] = Math.round(Math.expm1(m.data[i]))
  return { ...m, data }
}

/** log1p, for a matrix that is normalized but not logged. */
export function log1p(m: CellMajor): CellMajor {
  const data = new Float32Array(m.data.length)
  for (let i = 0; i < m.data.length; i++) data[i] = Math.log1p(m.data[i])
  return { ...m, data }
}

/**
 * Cell-major to gene-major, by counting sort.
 *
 * Two passes and no comparisons: count each gene's nonzeros, prefix-sum into
 * offsets, then scatter. Because cells are visited in order, each gene's cell
 * indices come out ascending — which the studio relies on and would otherwise
 * have to sort per gene on load.
 */
export function toGeneMajor(m: CellMajor): GeneMajor {
  const nnz = m.indptr[m.nCells]
  const indptr = new Int32Array(m.nGenes + 1)
  for (let k = 0; k < nnz; k++) indptr[m.indices[k] + 1]++
  for (let g = 0; g < m.nGenes; g++) indptr[g + 1] += indptr[g]

  const cursor = indptr.slice(0, m.nGenes)
  const indices = new Int32Array(nnz)
  const data = new Float32Array(nnz)
  for (let c = 0; c < m.nCells; c++) {
    for (let k = m.indptr[c]; k < m.indptr[c + 1]; k++) {
      const at = cursor[m.indices[k]]++
      indices[at] = c
      data[at] = m.data[k]
    }
  }
  return { indptr, indices, data }
}

/** Drop stored zeros — an explicit zero would be counted as a detected cell. */
export function dropZeros(m: CellMajor): CellMajor {
  let nnz = 0
  for (let i = 0; i < m.data.length; i++) if (m.data[i] !== 0) nnz++
  if (nnz === m.data.length) return m
  const indptr = new Int32Array(m.nCells + 1)
  const indices = new Int32Array(nnz)
  const data = new Float32Array(nnz)
  let at = 0
  for (let c = 0; c < m.nCells; c++) {
    indptr[c] = at
    for (let k = m.indptr[c]; k < m.indptr[c + 1]; k++) {
      if (m.data[k] === 0) continue
      indices[at] = m.indices[k]
      data[at] = m.data[k]
      at++
    }
  }
  indptr[m.nCells] = at
  return { ...m, indptr, indices, data }
}

/** Per-cell total and number of detected genes, for the QC block. */
export function cellQC(m: CellMajor): { total: Float32Array; genes: Float32Array } {
  const total = new Float32Array(m.nCells)
  const genes = new Float32Array(m.nCells)
  for (let c = 0; c < m.nCells; c++) {
    let s = 0, g = 0
    for (let k = m.indptr[c]; k < m.indptr[c + 1]; k++) {
      if (m.data[k] === 0) continue
      s += m.data[k]
      g++
    }
    total[c] = s
    genes[c] = g
  }
  return { total, genes }
}

/** Summed counts for a set of cells — one pseudobulk column. */
export function sumCells(m: CellMajor, cells: Int32Array): Int32Array {
  const out = new Int32Array(m.nGenes)
  for (const c of cells) {
    for (let k = m.indptr[c]; k < m.indptr[c + 1]; k++) out[m.indices[k]] += m.data[k]
  }
  return out
}
