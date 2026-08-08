// Building the bundle.
//
// The format is scrnaseq-studio/bundle@1, documented in the studio's
// tools/BUNDLE.md and written identically by tools/export_h5ad.py and
// tools/export_seurat.R. This is the third writer of that format, so it is
// pinned by tests that read a bundle back and compare against the matrix it
// was built from.

import { zipSync, strToU8, type ZippableFile } from 'fflate'
import { CHUNK_GENES, writeChunked } from './chunked.ts'
import {
  cellQC, expm1Round, log1p, lognormalize, sumCells, toGeneMajor,
  type CellMajor, type MatrixKind,
} from './matrix.ts'
import { column, type Choices, type Column, type MatrixRef, type Scan } from './scan.ts'

export const SCHEMA = 'scrnaseq-studio/bundle@1'

export class BuildError extends Error {}

export interface BundleMeta {
  schema: string
  label: string
  source: string
  nCells: number
  nGenes: number
  nnz: number
  clusters: string[]
  samples: { id: string; condition: string }[]
  conditions: string[]
  embedding: string
  expression: string
  hasRawCounts: boolean
  /**
   * Genes per block in `expr.chunk.bin`.
   *
   * Readers must use this number rather than the CHUNK_GENES constant: a bundle
   * written by an older or a differently-tuned lab still has to open.
   */
  chunkGenes: number
  provenance: Record<string, string | null>
  notes: string[]
}

export interface BuildResult {
  zip: Uint8Array
  meta: BundleMeta
  /** Sample × cluster columns written, or 0 when there are no counts. */
  pseudobulkColumns: number
}

type Progress = (stage: string) => void

/**
 * Which matrix to show, and which to sum for pseudobulk.
 *
 * Ported from the Python exporter so all three writers agree. The order of
 * preference matters more than it looks: an object with scaled values in X and
 * log1p(counts) in .raw must come out using .raw, and picking by slot name
 * would pick X every time.
 */
export function chooseMatrices(scan: Scan, note: (s: string) => void): {
  display: CellMajor
  counts: CellMajor | null
  geneSet: string
  expression: string
} {
  const kinds = scan.matrices
  for (const m of kinds) note(`${m.key} looks like ${m.kind}`)

  const first = (k: MatrixKind): MatrixRef | undefined => kinds.find(m => m.kind === k)

  // Raw counts, recovered if they are only hiding behind a log.
  let counts: CellMajor | null = null
  let countsSet: string | null = null
  const asCounts = first('counts')
  const asLogCounts = first('log-counts')
  if (asCounts) {
    counts = asCounts.load()
    countsSet = asCounts.geneSet
    note(`raw counts taken from ${asCounts.key}`)
  } else if (asLogCounts) {
    counts = expm1Round(asLogCounts.load())
    countsSet = asLogCounts.geneSet
    note(`${asLogCounts.key} is log1p(counts); exponentiated back to integer counts, so pseudobulk is available`)
  }

  // Display matrix: prefer one that is genuinely log-normalized already.
  const asLognorm = first('lognorm')
  if (asLognorm) {
    note(`expression taken from ${asLognorm.key}, already log-normalized`)
    const display = asLognorm.load()
    if (counts && countsSet !== asLognorm.geneSet) {
      note('the raw counts index a different gene list from the expression matrix; pseudobulk disabled rather than mismatched')
      counts = null
    }
    return { display, counts, geneSet: asLognorm.geneSet, expression: 'log1p(CP10K)' }
  }

  if (counts && countsSet) {
    note('expression computed here as log1p(CP10K) from the counts')
    return { display: lognormalize(counts), counts, geneSet: countsSet, expression: 'log1p(CP10K)' }
  }

  const asLinear = first('linear')
  if (asLinear) {
    note(`${asLinear.key} is normalized but not logged; log1p applied`)
    return { display: log1p(asLinear.load()), counts: null, geneSet: asLinear.geneSet, expression: 'log1p(CP10K)' }
  }

  throw new BuildError(
    `every matrix in this object holds scaled data (${kinds.map(m => m.key).join(', ')}), and the studio cannot plot z-scores. Re-save with the log-normalized or raw-count matrix included — in Seurat that is the @counts or @data slot, in Scanpy .raw or layers["counts"].`)
}

/** Codes and level names for a role, falling back to one level for the whole object. */
function grouping(col: Column | undefined, fallback: string, nCells: number): {
  codes: Int32Array; levels: string[]
} {
  if (!col || col.kind !== 'categorical') {
    return { codes: new Int32Array(nCells), levels: [fallback] }
  }
  // NA becomes its own level rather than silently joining level 0 — a cell with
  // no annotation is a fact about the object, not a member of the first cluster.
  const NA = col.levels.length
  const raw = new Int32Array(col.codes.length)
  let hasNA = false
  for (let i = 0; i < col.codes.length; i++) {
    if (col.codes[i] < 0) { raw[i] = NA; hasNA = true } else raw[i] = col.codes[i]
  }
  const all = hasNA ? [...col.levels, 'NA'] : col.levels

  // Levels with no cells are dropped. A bundle should describe the cells it
  // actually holds: when an atlas is split, most shards contain only a few of
  // the parent's 93 samples and a handful of its 20 timepoints, and carrying
  // the empty ones through made the studio report "0 samples per group" and a
  // minimum replicate count of zero for groups that were simply not there.
  const used = new Uint8Array(all.length)
  for (const c of raw) used[c] = 1
  let n = 0
  const remap = new Int32Array(all.length).fill(-1)
  const levels: string[] = []
  for (let i = 0; i < all.length; i++) {
    if (!used[i]) continue
    remap[i] = n++
    levels.push(all[i])
  }
  if (levels.length === all.length) return { codes: raw, levels: all }
  const codes = new Int32Array(raw.length)
  for (let i = 0; i < raw.length; i++) codes[i] = remap[raw[i]]
  return { codes, levels }
}

/** A numeric obs column by any of several conventional names. */
function qcColumn(scan: Scan, names: string[]): Float32Array | null {
  for (const n of names) {
    const c = scan.columns.find(x => x.name.toLowerCase() === n.toLowerCase())
    if (c?.kind === 'numeric' && c.values.length === scan.nCells) return c.values
  }
  return null
}

export function buildBundle(
  scan: Scan, choices: Choices, onProgress: Progress = () => {},
): BuildResult {
  const notes = [...scan.notes]
  const note = (s: string) => notes.push(s)
  const n = scan.nCells

  const clusterCol = column(scan, choices.cluster)
  if (!clusterCol || clusterCol.kind !== 'categorical') {
    throw new BuildError(`"${choices.cluster}" is not a usable cell annotation`)
  }

  onProgress('Reading the expression matrix')
  const { display, counts, geneSet, expression } = chooseMatrices(scan, note)

  const genes = scan.geneSets[geneSet]
  if (!genes) throw new BuildError('the chosen matrix has no gene names')
  if (display.nGenes !== genes.length) {
    throw new BuildError(`the matrix has ${display.nGenes} genes but the gene list has ${genes.length}`)
  }
  if (display.nCells !== n) {
    throw new BuildError(`the matrix has ${display.nCells} cells but obs has ${n}`)
  }

  const cluster = grouping(clusterCol, 'all cells', n)
  const sampleCol = column(scan, choices.sample)
  const condCol = column(scan, choices.condition)
  const sample = grouping(sampleCol, 'all cells', n)

  if (!sampleCol) {
    note('no sample column — every cell is treated as one sample. Composition cannot show between-animal spread, and pseudobulk needs several samples per group, so only the per-cell test will be available')
  }

  // One condition per sample: the bundle's model, so it is resolved here rather
  // than left for the studio to discover.
  let conditions: string[]
  let sampleCond: string[]
  if (condCol && condCol.kind === 'categorical') {
    const cond = grouping(condCol, 'all cells', n)
    conditions = cond.levels
    sampleCond = sample.levels.map((_id, si) => {
      for (let i = 0; i < n; i++) if (sample.codes[i] === si) return conditions[cond.codes[i]]
      return conditions[0]
    })
  } else {
    note('no condition column — this object opens as a single group, and every comparison tab stays empty rather than inventing a contrast')
    conditions = ['all cells']
    sampleCond = sample.levels.map(() => 'all cells')
  }

  const emb = scan.embeddings.find(e => e.key === choices.embedding)
  if (!emb) throw new BuildError(`no embedding called "${choices.embedding}" in this object`)
  onProgress('Reading the embedding')
  const xy = emb.load()
  if (/pca/i.test(emb.key)) {
    note('using the first two principal components, which is a much coarser picture than a UMAP')
  }

  // ---- QC ---------------------------------------------------------------
  onProgress('Computing QC')
  let total = qcColumn(scan, ['n_counts', 'total_counts', 'nCount_RNA'])
  let nGene = qcColumn(scan, ['n_genes', 'n_genes_by_counts', 'nFeature_RNA'])
  let mito = qcColumn(scan, ['percent_mito', 'pct_counts_mt', 'percent.mt', 'percent_mt'])

  if (!total || !nGene) {
    const q = cellQC(counts ?? display)
    if (!total) { total = q.total; note('no total-count column; recomputed from the matrix') }
    if (!nGene) { nGene = q.genes; note('no detected-gene column; recomputed from the matrix') }
  }
  if (!mito) {
    mito = new Float32Array(n)
    note('no mitochondrial fraction in the metadata — the QC panel will show a flat zero rather than a made-up number')
  } else {
    let max = 0
    for (const v of mito) if (Number.isFinite(v) && v > max) max = v
    if (max <= 1) {
      const scaled = new Float32Array(n)
      for (let i = 0; i < n; i++) scaled[i] = mito[i] * 100
      mito = scaled
      note('mitochondrial fraction was a proportion; converted to a percentage')
    }
  }

  // ---- pack -------------------------------------------------------------
  onProgress('Transposing to gene-major')
  const csc = toGeneMajor(display)

  // The same matrix a second time, cut into per-gene blocks, so the studio can
  // plot one gene by slicing a few kilobytes out of the file instead of
  // inflating the whole matrix. See chunked.ts: it costs about half again on
  // top of the flat entries, and it is what makes a large part openable at all.
  onProgress('Indexing genes for lazy reading')
  const chunked = writeChunked(csc.indptr, csc.indices, csc.data, CHUNK_GENES)

  const qc = new Float32Array(3 * n)
  for (let i = 0; i < n; i++) {
    qc[3 * i] = total[i]
    qc[3 * i + 1] = nGene[i]
    qc[3 * i + 2] = mito[i]
  }

  const clusterCodes = new Uint16Array(n)
  const sampleCodes = new Uint16Array(n)
  for (let i = 0; i < n; i++) {
    clusterCodes[i] = cluster.codes[i]
    sampleCodes[i] = sample.codes[i]
  }
  if (cluster.levels.length > 65535 || sample.levels.length > 65535) {
    throw new BuildError('more than 65 535 clusters or samples — that is not an annotation')
  }

  const meta: BundleMeta = {
    schema: SCHEMA,
    label: choices.label || scan.source,
    source: `${scan.source} (${scan.format === 'h5ad' ? 'AnnData' : 'Seurat'}, converted in scRNA-seq Lab)`,
    nCells: n,
    nGenes: genes.length,
    nnz: csc.data.length,
    clusters: cluster.levels,
    samples: sample.levels.map((id, i) => ({ id, condition: sampleCond[i] })),
    conditions,
    embedding: emb.key,
    expression,
    hasRawCounts: counts != null,
    chunkGenes: CHUNK_GENES,
    provenance: { ...scan.provenance, clustering: choices.cluster },
    notes,
  }

  const files: Record<string, ZippableFile> = {
    'meta.json': strToU8(JSON.stringify(meta, null, 1)),
    'genes.txt': strToU8(genes.join('\n')),
    'cluster.u16': new Uint8Array(clusterCodes.buffer),
    'sample.u16': new Uint8Array(sampleCodes.buffer),
    'embed.f32': new Uint8Array(xy.buffer, xy.byteOffset, xy.byteLength),
    'qc.f32': new Uint8Array(qc.buffer),
    'expr.indptr.i32': new Uint8Array(csc.indptr.buffer),
    'expr.indices.i32': new Uint8Array(csc.indices.buffer),
    'expr.data.f32': new Uint8Array(csc.data.buffer),
    'expr.chunkptr.i32': new Uint8Array(chunked.ptr.buffer),
    // STORED, never deflated. Each chunk is already its own deflate stream; the
    // whole point is that chunk k occupies a known contiguous byte range of the
    // zip, so a reader can slice it out. Compressing the entry would put every
    // chunk behind one stream that has to be read from its start — which is the
    // problem this format exists to solve.
    'expr.chunk.bin': [chunked.bin, { level: 0 }],
  }

  // ---- pseudobulk -------------------------------------------------------
  let pbColumns = 0
  if (counts && counts.nGenes === genes.length) {
    onProgress('Summing pseudobulk')
    const cols: Int32Array[] = []
    const names: string[] = []
    for (let si = 0; si < sample.levels.length; si++) {
      for (let ci = 0; ci < cluster.levels.length; ci++) {
        const sel: number[] = []
        for (let i = 0; i < n; i++) {
          if (sample.codes[i] === si && cluster.codes[i] === ci) sel.push(i)
        }
        if (!sel.length) continue
        cols.push(sumCells(counts, Int32Array.from(sel)))
        names.push(`${sample.levels[si]}||${cluster.levels[ci]}||${sel.length}`)
      }
    }
    if (cols.length) {
      const lines: string[] = [`gene\t${names.join('\t')}`]
      for (let gi = 0; gi < genes.length; gi++) {
        const row = new Array<string>(cols.length + 1)
        row[0] = genes[gi]
        for (let c = 0; c < cols.length; c++) row[c + 1] = String(cols[c][gi])
        lines.push(row.join('\t'))
      }
      files['pseudobulk.tsv'] = strToU8(lines.join('\n'))
      pbColumns = cols.length
    }
  }

  onProgress('Compressing')
  const zip = zipSync(files, { level: 6 })
  return { zip, meta, pseudobulkColumns: pbColumns }
}
