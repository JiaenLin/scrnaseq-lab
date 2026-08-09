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
import {
  column, geneNameKind,
  type Choices, type Column, type GeneIdKind, type MatrixRef, type Scan,
} from './scan.ts'

export const SCHEMA = 'scrnaseq-studio/bundle@1'

export class BuildError extends Error {}

/** One 2D embedding inside the bundle. */
export interface EmbeddingEntry {
  /** The name the object gave it, e.g. `X_UMAP`. */
  key: string
  /** Entry holding interleaved x,y — float32, 2 × nCells. */
  file: string
}

/**
 * The gene names the matrix rows are *not* indexed by.
 *
 * Written when the object carried both namings. Aligned with `genes.txt`, one
 * line per matrix row, so the studio converts by index and never by lookup.
 */
export interface GeneAliasMeta {
  /** What `gene_alias.txt` holds — the opposite of `geneIdKind`. */
  kind: 'symbol' | 'accession'
  /** The var / feature column it was read from, for the record. */
  column: string
  /** Entry name. */
  file: string
  /**
   * Rows the object had no alias for. Those lines repeat `genes.txt`, so the
   * file is always a usable display name and never a blank.
   */
  missing: number
  /**
   * Rows whose alias is shared with another row.
   *
   * Nothing is merged: several accessions really do map to one symbol, and
   * collapsing them would silently sum two genes' expression. The studio has to
   * decide how to show a repeated name; the bundle keeps the rows apart.
   */
  duplicated: number
}

/**
 * One further categorical column, carried per cell under its own name.
 *
 * The bundle has always had exactly three of these — cluster, sample,
 * condition — because those are the three every view depends on. This is the
 * rest of what the object knows: a dissection beside an age, a Class above a
 * Subclass. Nothing here is a role, so nothing is renamed: `key` is the obs /
 * meta.data column the levels came out of, and that is what the studio's menus
 * say. A reader that does not know about this field opens the bundle unchanged.
 */
export interface ExtraColumn {
  /** The column's own name in the object, e.g. `dissection`. */
  key: string
  /** Entry holding one uint16 level index per cell. */
  file: string
  /** Level names, in the object's own order, dropped to the ones this part uses. */
  levels: string[]
}

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
  /**
   * Further categorical columns, in the order they were chosen.
   *
   * Absent in bundles written before today, which is the same as empty.
   */
  extras?: ExtraColumn[]
  /** The default embedding's key — `embeddings[0].key`. */
  embedding: string
  /**
   * Every 2D embedding the object had, the default first.
   *
   * Absent in bundles written before today; a reader without it has exactly one
   * embedding, `embed.f32`, named by `embedding`.
   */
  embeddings?: EmbeddingEntry[]
  /** What `genes.txt` holds. Absent in older bundles, where it is unknown. */
  geneIdKind?: GeneIdKind
  /** The other naming, when the object had one. */
  geneAlias?: GeneAliasMeta | null
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
/**
 * Below this share of the gene list, a matrix is not the object's expression.
 *
 * Deliberately far below anything a real matrix reaches. The case this exists
 * for is not marginal: on the developing mouse atlas `layers["expected"]` — a
 * velocity model's output — holds values in exactly 2 000 of 31 053 genes, 6%,
 * while `X` reaches most of the list.
 */
const COVERAGE_FLOOR = 0.25

/**
 * Drop matrices that only ever touch a corner of the gene list.
 *
 * Classification by value cannot catch this. A model's expected counts are
 * small positive non-integers, which is precisely the signature of
 * log-normalized expression, so it wins the preference order and the studio
 * then reports 29 053 of 31 053 genes as expressed in no cell anywhere — every
 * marker panel silently empty, every gene set matching nothing that plots.
 *
 * A candidate is dropped only when it is far below the floor *and* another
 * candidate covers at least twice as much, so an object whose only matrix is
 * sparse keeps it: the point is to choose between matrices, never to refuse the
 * object.
 */
export function usableCoverage(
  matrices: MatrixRef[], note: (s: string) => void,
): MatrixRef[] {
  const cov = new Map<MatrixRef, number>()
  for (const m of matrices) {
    if (!m.coverage) continue
    try {
      cov.set(m, m.coverage())
    } catch { /* an unreadable index is not a reason to reject a matrix */ }
  }
  if (cov.size < 2) return matrices
  const best = Math.max(...cov.values())
  const out = matrices.filter(m => {
    const c = cov.get(m)
    if (c === undefined || c >= COVERAGE_FLOOR || c * 2 > best) return true
    note(`${m.key} has values for only ${(c * 100).toFixed(0)}% of the genes, against ${(best * 100).toFixed(0)}% elsewhere — it is a model over selected features, not this object's expression, so it is not used`)
    return false
  })
  return out.length ? out : matrices
}

export function chooseMatrices(scan: Scan, note: (s: string) => void): {
  display: CellMajor
  counts: CellMajor | null
  geneSet: string
  expression: string
} {
  const kinds = usableCoverage(scan.matrices, note)
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


/**
 * Which matrix the split path should stream.
 *
 * `chooseMatrices` decides by loading candidates, which a 700 M-value object
 * cannot afford, so this makes the same decision from the classifications
 * alone. It has to agree: the split path used to take the first matrix that
 * could be streamed, which on the commonest atlas layout — a scaled `X` next to
 * raw counts in `.raw` — is the scaled one. That either ships z-scores as
 * expression or dies at the end of a twenty-minute conversion complaining that
 * every matrix is scaled, on a file that plainly is not.
 *
 * Preference matches `chooseMatrices`: an already log-normalized matrix first,
 * then counts it can build one from. Only one matrix is streamed, so when the
 * object keeps expression and counts apart, the display matrix wins and
 * pseudobulk is left out rather than paying a second pass over the whole file.
 */
export function streamableMatrix(scan: Scan, note: (s: string) => void = () => {}): MatrixRef | null {
  const usable = usableCoverage(scan.matrices, note)
    .filter(m => m.loadGroups && m.kind !== 'scaled')
  for (const kind of ['lognorm', 'counts', 'log-counts', 'linear'] as const) {
    const hit = usable.find(m => m.kind === kind)
    if (hit) return hit
  }
  return null
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

/** An entry name that survives a zip and a file system. */
const safeEntry = (s: string) =>
  (s || 'x').replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '').slice(0, 40) || 'x'

/**
 * The alias list as it goes into the bundle, and what had to be decided.
 *
 * Two things the object does to this list, both of which have to be said out
 * loud rather than papered over:
 *
 *   - an empty symbol. The row keeps its accession, so every line of the file
 *     is a name something can be labelled with.
 *   - a symbol shared by several accessions. Real: read-through transcripts and
 *     paralogue annotations both do it. The rows stay separate — merging them
 *     would add two genes' expression together under one name — and the count
 *     is recorded so the studio can disambiguate rather than pick one.
 */
function alignAlias(genes: string[], alias: { names: string[] }): {
  names: string[]; missing: number; duplicated: number
} {
  const names = new Array<string>(genes.length)
  let missing = 0
  const seen = new Map<string, number>()
  for (let i = 0; i < genes.length; i++) {
    const raw = (alias.names[i] ?? '').trim()
    const ok = raw && raw !== 'NA' && raw !== 'nan'
    if (!ok) missing++
    names[i] = ok ? raw : genes[i]
    seen.set(names[i], (seen.get(names[i]) ?? 0) + 1)
  }
  let duplicated = 0
  for (const nm of names) if ((seen.get(nm) ?? 0) > 1) duplicated++
  return { names, missing, duplicated }
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

  // ---- the columns beyond the three roles --------------------------------
  // Through the same `grouping` as cluster and sample, which is the whole point:
  // a shard drops the levels it holds no cells for and renumbers its codes. A
  // part that inherited the parent's level list would make its "dissection 3"
  // a different region from the next part's, and every figure would still draw.
  const roles = [choices.cluster, choices.sample, choices.condition]
  const extras: ExtraColumn[] = []
  const extraCodes: { file: string; codes: Uint16Array }[] = []
  const extraFiles = new Set<string>()
  for (const name of choices.extras ?? []) {
    const col = column(scan, name)
    if (roles.includes(name)) continue     // already carried, under its role
    if (!col || col.kind !== 'categorical') {
      note(`"${name}" is not a categorical column, so it is not carried as an extra grouping`)
      continue
    }
    const g = grouping(col, 'all cells', n)
    if (g.levels.length > 65535) {
      note(`${col.name} has ${g.levels.length} levels, too many to store as a column`)
      continue
    }
    let file = `extra.${safeEntry(col.name)}.u16`
    for (let i = 2; extraFiles.has(file); i++) file = `extra.${safeEntry(col.name)}-${i}.u16`
    extraFiles.add(file)
    const codes = new Uint16Array(n)
    for (let i = 0; i < n; i++) codes[i] = g.codes[i]
    extras.push({ key: col.name, file, levels: g.levels })
    extraCodes.push({ file, codes })
    // The parent's level count, not this part's: subsetting a column keeps its
    // levels whole, so this sentence reads the same in every part of a split.
    note(`${col.name} travels with the cells as an extra grouping — ${col.levels.length} levels the studio can break a figure down by`)
  }

  const emb = scan.embeddings.find(e => e.key === choices.embedding)
  if (!emb) throw new BuildError(`no embedding called "${choices.embedding}" in this object`)
  onProgress('Reading the embedding')
  const xy = emb.load()
  if (/pca/i.test(emb.key)) {
    note('using the first two principal components, which is a much coarser picture than a UMAP')
  }

  // Every other 2D embedding as well. An object routinely holds a UMAP and a
  // t-SNE of the same cells, and until now the lab kept one and threw the rest
  // away — 8 bytes a cell to carry, against a matrix that costs 8 bytes a
  // *stored value*. The chosen one stays first and stays in `embed.f32`, so a
  // reader that knows nothing about this sees no change at all.
  const others = scan.embeddings.filter(e => e.key !== emb.key && e.nDims >= 2)
  const embedFiles: EmbeddingEntry[] = [{ key: emb.key, file: 'embed.f32' }]
  const extraXY: { file: string; xy: Float32Array }[] = []
  const usedNames = new Set(['embed.f32'])
  for (const e of others) {
    let file = `embed.${safeEntry(e.key)}.f32`
    for (let i = 2; usedNames.has(file); i++) file = `embed.${safeEntry(e.key)}-${i}.f32`
    usedNames.add(file)
    onProgress(`Reading the ${e.key} embedding`)
    let coords: Float32Array
    try {
      coords = e.load()
    } catch (err) {
      // One unreadable extra embedding must not lose the conversion: the
      // default one is already in hand and is what the studio opens with.
      note(`${e.key} could not be read (${(err as Error).message}); it is not in the bundle`)
      continue
    }
    if (coords.length !== 2 * n) {
      note(`${e.key} has ${coords.length / 2} points for ${n} cells; it is not in the bundle`)
      continue
    }
    extraXY.push({ file, xy: coords })
    embedFiles.push({ key: e.key, file })
  }
  if (extraXY.length) {
    note(`${embedFiles.length} embeddings carried (${embedFiles.map(e => e.key).join(', ')}); ${emb.key} opens by default`)
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

  // ---- gene names -------------------------------------------------------
  // Which naming the matrix rows carry, and the other one beside it. The
  // built-in gene sets in the studio are symbols; on this atlas the rows are
  // Ensembl accessions, which is the entire reason those sets matched nothing.
  const geneIdKind = geneNameKind(genes)
  const rawAlias = scan.geneAliases?.[geneSet]
  let geneAlias: GeneAliasMeta | null = null
  let aliasText: string | null = null
  if (rawAlias && rawAlias.names.length === genes.length) {
    const a = alignAlias(genes, rawAlias)
    geneAlias = {
      kind: rawAlias.kind, column: rawAlias.column, file: 'gene_alias.txt',
      missing: a.missing, duplicated: a.duplicated,
    }
    aliasText = a.names.join('\n')
    note(`gene names are ${geneIdKind === 'accession' ? 'accessions' : geneIdKind === 'symbol' ? 'symbols' : 'a mix of accessions and symbols'}; var/${rawAlias.column} carries the ${rawAlias.kind}s and travels with them`)
    if (a.missing) {
      note(`${a.missing} of ${genes.length} genes have no ${rawAlias.kind} in the object; those keep the name the matrix is indexed by`)
    }
    if (a.duplicated) {
      note(`${a.duplicated} genes share a ${rawAlias.kind} with another gene; the rows are kept separate rather than summed`)
    }
  } else if (geneIdKind === 'accession') {
    note('gene names are accessions and this object holds no symbol column, so gene sets and marker lists written as symbols will not match')
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
    extras,
    embedding: emb.key,
    embeddings: embedFiles,
    geneIdKind,
    geneAlias,
    expression,
    hasRawCounts: counts != null,
    chunkGenes: CHUNK_GENES,
    // Which column each role was read from. `clustering` has always been here;
    // the other two are recorded for the same reason — the studio's menus say
    // "Group" because they have to say something, and an object that calls it
    // Age should have the menu say Age. null means the object had no such
    // column, which is a different fact from not knowing.
    provenance: {
      ...scan.provenance,
      clustering: choices.cluster,
      condition: condCol?.kind === 'categorical' ? condCol.name : null,
      sample: sampleCol?.kind === 'categorical' ? sampleCol.name : null,
    },
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

  for (const e of extraXY) {
    files[e.file] = new Uint8Array(e.xy.buffer, e.xy.byteOffset, e.xy.byteLength)
  }
  for (const e of extraCodes) files[e.file] = new Uint8Array(e.codes.buffer)
  if (aliasText != null) files['gene_alias.txt'] = strToU8(aliasText)

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
