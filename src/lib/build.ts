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
  type Choices, type Column, type GeneIdKind, type MatrixRef, type RunSection, type Scan,
} from './scan.ts'
import {
  FactError, factsOf, measured, notRecorded, placeOf, recorded,
  type Fact, type MetaFacts,
} from './facts.ts'
import {
  COVERAGE_FLOOR, choiceOf, measureExpression, measureObject, pickMatrices, plan, planInfo, rowsOf, sampleRuns,
  type Measurements, type Plan, type PlanInfo,
} from './plan.ts'

export const SCHEMA = 'scrnaseq-studio/bundle@1'

export class BuildError extends Error {}

/**
 * The modification time stamped on every entry of a bundle.
 *
 * Fixed, so that converting one object twice gives the same bytes and a sha256
 * can say that two conversions agree. Left to the default, fflate stamps the
 * current time, and two conversions of the same object differed in nothing
 * but that.
 *
 * Built from LOCAL date parts, not from an epoch, because fflate encodes the
 * DOS date and time from local-time getters: 1980-01-01 12:00 constructed in
 * any time zone encodes to the same two words. Noon, so that no zone's
 * offset change near midnight can move it to another day. 1980 is the first
 * year a zip entry can carry.
 */
const ZIP_MTIME = new Date(1980, 0, 1, 12, 0, 0)

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
  /**
   * Cluster levels that are not cell types: an annotator's labels for cells it
   * excluded or could not resolve. They stay in `clusters` — the cells are
   * real and are drawn — and are named here so a reader can leave them out of
   * composition, markers and tests.
   *
   * Written by the command-line converter, which always says (an empty list is
   * "none"); absent in bundles from the app and from before it existed.
   */
  sentinels?: string[]
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
   * The embeddings a reader can switch between, the default first.
   *
   * Besides the default, only two-column embeddings, and each set of
   * coordinates once: a latent space (a PCA, a corrected PCA) cut to its first
   * two columns is not a picture of the cells, and a byte copy of another
   * embedding is the same picture under a second name. The notes name what
   * was left out and why.
   *
   * Absent in bundles written before more than one was carried; a reader
   * without it has exactly one embedding, `embed.f32`, named by `embedding`.
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
  /**
   * What the bundle states about its object, each with its source: recorded
   * (with its place), measured (of what, and how) or chosen (by whom). Written
   * only by walking the table in facts.ts; the ids and shapes are the Studio's
   * (its tools/BUNDLE.md, "facts"). bundle@1 readers ignore it.
   */
  facts?: MetaFacts
  /**
   * The run's identity, not its record: from each record the object's uns holds
   * (a tool's provenance, the harness's run record), the fields tool, version,
   * commit, stage, run_key and job it stores as strings of 256 characters or
   * fewer, and nothing else. `records: []` means the converter looked and the
   * object carries none.
   *
   * Absent in bundles from a Seurat object, which has no uns, and in bundles
   * written before it existed. bundle@1 readers ignore it.
   */
  run?: RunSection
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
  /** The matrices by their keys in the object, and what the expression was made from. */
  displayKey: string
  countsKey: string | null
  from: 'lognorm' | 'counts' | 'log-counts' | 'linear'
  /** What the counts were read as: stored counts, or log1p(counts) exponentiated back. */
  countsKind: 'counts' | 'log-counts' | null
  /** The linear matrix as stored, before the log1p the expression takes; only when `from` is linear. */
  stored: CellMajor | null
} {
  const kinds = usableCoverage(scan.matrices, note)
  for (const m of kinds) note(`${m.key} looks like ${m.kind}`)

  const first = (k: MatrixKind): MatrixRef | undefined => kinds.find(m => m.kind === k)

  // The plan's rule (plan.ts, pickMatrices) decides; this loads what it decided.
  const decided = pickMatrices(kinds)
  const agree = (key: string, cKey: string | null) => {
    if (decided.display?.key !== key || (decided.counts?.key ?? null) !== cKey) {
      throw new BuildError(`the plan and the build chose different matrices (${decided.display?.key} and ${key}): this is a bug in the converter`)
    }
  }

  // Raw counts, recovered if they are only hiding behind a log.
  let counts: CellMajor | null = null
  let countsSet: string | null = null
  const asCounts = first('counts')
  const asLogCounts = first('log-counts')
  const countsRef = asCounts ?? asLogCounts
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
    agree(asLognorm.key, counts ? countsRef!.key : null)
    return { display, counts, geneSet: asLognorm.geneSet, expression: 'log1p(CP10K)',
      displayKey: asLognorm.key, countsKey: counts ? countsRef!.key : null, from: 'lognorm',
      countsKind: counts ? countsRef!.kind as 'counts' | 'log-counts' : null, stored: null }
  }

  if (counts && countsSet) {
    note('expression computed here as log1p(CP10K) from the counts')
    agree(countsRef!.key, countsRef!.key)
    return { display: lognormalize(counts), counts, geneSet: countsSet, expression: 'log1p(CP10K)',
      displayKey: countsRef!.key, countsKey: countsRef!.key, from: countsRef!.kind as 'counts' | 'log-counts',
      countsKind: countsRef!.kind as 'counts' | 'log-counts', stored: null }
  }

  const asLinear = first('linear')
  if (asLinear) {
    note(`${asLinear.key} is normalized but not logged; log1p applied`)
    agree(asLinear.key, null)
    const stored = asLinear.load()
    return { display: log1p(stored), counts: null, geneSet: asLinear.geneSet, expression: 'log1p(normalized)',
      displayKey: asLinear.key, countsKey: null, from: 'linear', countsKind: null, stored }
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

/** Codes and the level names they index, every code a level a cell holds. */
type Grouping = { codes: Int32Array; levels: string[] }

/** "a", "a and b", "a, b and c" — at most ten named, then a count. */
function list(names: string[], max = 10): string {
  const shown = names.length > max ? [...names.slice(0, max), `${names.length - max} more`] : names
  return shown.length < 2 ? shown.join('') : `${shown.slice(0, -1).join(', ')} and ${shown[shown.length - 1]}`
}

/**
 * Split each sample that holds more than one condition into one sample per
 * condition, named `sample|condition`.
 *
 * Only a split sample is renamed, so one animal in two groups does not rename
 * every other. New names come in the sample's level order, then the
 * condition's. A name some other sample already has gets `#2`, `#3`, … rather
 * than merging two samples. With `everyOne`, every sample is split by
 * condition whether it spans or not: that is the no-sample-column case, where
 * the one sample is the whole object and a part of a split object may hold
 * only one condition, yet must name its cells as every other part does.
 *
 * Null when nothing is split.
 */
function splitByCondition(sample: Grouping, cond: Grouping, everyOne = false): (Grouping & {
  /** The samples that were split, by their old names. */
  split: string[]
  /** The names they became. */
  made: string[]
}) | null {
  const n = sample.codes.length
  const nS = sample.levels.length
  const first = new Int32Array(nS).fill(-1)
  const spans = new Uint8Array(nS)
  for (let i = 0; i < n; i++) {
    const s = sample.codes[i], c = cond.codes[i]
    if (first[s] < 0) first[s] = c
    else if (first[s] !== c) spans[s] = 1
  }
  if (everyOne) for (let s = 0; s < nS; s++) if (first[s] >= 0) spans[s] = 1
  if (!spans.some(x => x)) return null

  // The conditions each split sample holds, in the condition's level order.
  const held = new Map<number, Uint8Array>()
  for (let i = 0; i < n; i++) {
    const s = sample.codes[i]
    if (!spans[s]) continue
    let h = held.get(s)
    if (!h) held.set(s, h = new Uint8Array(cond.levels.length))
    h[cond.codes[i]] = 1
  }
  const taken = new Set(sample.levels.filter((_l, s) => !spans[s]))
  const levels: string[] = []
  const keptAt = new Int32Array(nS).fill(-1)
  const splitAt = new Map<number, Int32Array>()
  const split: string[] = []
  const made: string[] = []
  for (let s = 0; s < nS; s++) {
    if (!spans[s]) {
      keptAt[s] = levels.length
      levels.push(sample.levels[s])
      continue
    }
    split.push(sample.levels[s])
    const at = new Int32Array(cond.levels.length).fill(-1)
    const h = held.get(s)!
    for (let c = 0; c < cond.levels.length; c++) {
      if (!h[c]) continue
      const base = `${sample.levels[s]}|${cond.levels[c]}`
      let id = base
      for (let k = 2; taken.has(id); k++) id = `${base}#${k}`
      taken.add(id)
      at[c] = levels.length
      levels.push(id)
      made.push(id)
    }
    splitAt.set(s, at)
  }
  const codes = new Int32Array(n)
  for (let i = 0; i < n; i++) {
    const s = sample.codes[i]
    codes[i] = spans[s] ? splitAt.get(s)![cond.codes[i]] : keptAt[s]
  }
  return { codes, levels, split, made }
}

/**
 * One condition per sample, which is what the bundle records.
 *
 * The bundle stores a condition per SAMPLE, not per cell. A sample whose cells
 * sit in two conditions — a hashed lane carrying two treatments, a donor
 * sampled before and after — used to have every cell given its first cell's
 * condition: per-group totals that still looked plausible, and nothing to say
 * so. Such a sample is two samples. Splitting it keeps every cell's condition,
 * needs no change to the format, and is the rule the studio's own exporter
 * applies (tools/export_h5ad.py), so the writers of bundle@1 agree.
 *
 * Decided on the WHOLE object. The studio merges the parts of a split object
 * by sample id, and a part that held only one condition of a sample would, on
 * its own, leave it unsplit: "s1" in one part, "s1|ctrl" in the next, one
 * animal counted as two replicates. So the split paths call this on the parent
 * before cutting it into parts, and the parts inherit the result. buildBundle
 * calls it too, which costs nothing on a scan already resolved and covers every
 * caller that converts a whole object at once.
 *
 * The sample column is replaced under its own name, so `provenance.sample`
 * still names the column the samples came from. Returns the same Scan when
 * every sample already sits in one condition.
 */
export function oneConditionPerSample(scan: Scan, choices: Choices): Scan {
  const sampleCol = column(scan, choices.sample)
  const condCol = column(scan, choices.condition)
  if (sampleCol?.kind !== 'categorical' || condCol?.kind !== 'categorical') return scan
  const sample = grouping(sampleCol, 'all cells', scan.nCells)
  const cond = grouping(condCol, 'all cells', scan.nCells)
  const s = splitByCondition(sample, cond)
  if (!s) return scan
  const replaced: Column = {
    name: sampleCol.name, kind: 'categorical', levels: s.levels, codes: s.codes,
    values: new Float32Array(0), missing: 0,
  }
  return {
    ...scan,
    columns: scan.columns.map(c => (c === sampleCol ? replaced : c)),
    notes: [...scan.notes,
      `${list(s.split)} ${s.split.length > 1 ? 'hold' : 'holds'} cells in more than one level of ${condCol.name}, so each was split into one sample per level (${list(s.made)}): ${sample.levels.length} samples became ${s.levels.length}. The bundle records one condition per sample, and giving those cells their sample's first condition would have put them in the wrong group with no visible sign`],
  }
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

/** A numeric obs column at one of several conventional keys: read there, and named. */
function qcColumn(scan: Scan, names: string[]): { name: string; values: Float32Array } | null {
  for (const n of names) {
    const c = scan.columns.find(x => x.name.toLowerCase() === n.toLowerCase())
    if (c?.kind === 'numeric' && c.values.length === scan.nCells) return { name: c.name, values: c.values }
  }
  return null
}

/** A per-cell number of the QC block, and the fact that says what it is: made together, in one place. */
export interface PerCell { values: Float32Array; fact: Fact | undefined }

/** The matrices chooseMatrices loaded, as perCellQC needs them. */
export type Summable = Pick<ReturnType<typeof chooseMatrices>, 'display' | 'counts' | 'displayKey' | 'countsKey' | 'countsKind' | 'from' | 'stored'>

/**
 * The per-cell total and genes detected: THE one place that computes them, and
 * the only place their facts' words are made, from the branch that made the
 * numbers. Each is an obs column read at a standard key (recorded, at its
 * place), or computed here from ONE matrix, and then the fact names that matrix
 * and what was done to its stored values before they were summed or counted:
 * the counts as stored; the counts a log1p(counts) matrix holds, exponentiated
 * back and rounded; or, with no counts, the expression matrix as STORED (a
 * linear one before the log1p the bundle's expression takes; a logged one as
 * the logs it holds).
 */
export function perCellQC(scan: Scan, m: Summable, note: (s: string) => void): { total: PerCell; genes: PerCell } {
  const totalCol = qcColumn(scan, ['n_counts', 'total_counts', 'nCount_RNA'])
  const nGeneCol = qcColumn(scan, ['n_genes', 'n_genes_by_counts', 'nFeature_RNA'])
  const src = m.counts
    ? { matrix: m.counts, of: placeOf(m.countsKey!),
      done: m.countsKind === 'log-counts' ? ', after expm1 of each stored value (its log1p(counts)), rounded to an integer' : '' }
    : { matrix: m.from === 'linear' && m.stored ? m.stored : m.display, of: placeOf(m.displayKey), done: '' }
  let q: { total: Float32Array; genes: Float32Array } | null = null
  const computed = () => (q ??= cellQC(src.matrix))
  const total: PerCell = totalCol
    ? { values: totalCol.values, fact: factOr(() => recorded(totalCol.name, `obs/${totalCol.name}`), note, `the total column's name`) }
    : { values: computed().total, fact: measured(null, src.of, `summed per cell over every gene${src.done}`) }
  const genes: PerCell = nGeneCol
    ? { values: nGeneCol.values, fact: factOr(() => recorded(nGeneCol.name, `obs/${nGeneCol.name}`), note, `the detected-gene column's name`) }
    : { values: computed().genes, fact: measured(null, src.of, `values above zero counted per cell${src.done}`) }
  if (!totalCol) note('no total-count column; recomputed from the matrix')
  if (!nGeneCol) note('no detected-gene column; recomputed from the matrix')
  return { total, genes }
}

/** A fact, or none and a note when the object's own text cannot be carried as one line. */
function factOr(make: () => Fact, note: (s: string) => void, what: string): Fact | undefined {
  try {
    return make()
  } catch (e) {
    if (e instanceof FactError) { note(`${what} is not one line of text a reader can show, so the bundle states no fact for it`); return undefined }
    throw e
  }
}

/** The plan and the whole object's measurements, made here when the caller did not (a whole object). */
export function planOf(scan: Scan, choices: Choices): { plan: Plan; measured: Measurements } {
  const p = plan({ nCells: scan.nCells, columns: scan.columns, ...planInfo(scan, null) } as PlanInfo, choices)
  return { plan: p, measured: measureObject(scan, p, choices) }
}

export function buildBundle(
  given: Scan, choices: Choices, onProgress: Progress = () => {},
  /** The parent's plan and measurements: given for every part of a split, made here otherwise. */
  ctx?: { plan: Plan; measured: Measurements },
): BuildResult {
  // A no-op on a part whose parent was already resolved; see the function.
  const scan = oneConditionPerSample(given, choices)
  const notes = [...scan.notes]
  const note = (s: string) => notes.push(s)
  const n = scan.nCells
  const { plan: p, measured: ms } = ctx ?? planOf(scan, choices)

  const clusterCol = column(scan, choices.cluster)
  if (!clusterCol || clusterCol.kind !== 'categorical') {
    throw new BuildError(`"${choices.cluster}" is not a usable cell annotation`)
  }

  onProgress('Reading the expression matrix')
  const chosenMatrices = chooseMatrices(scan, note)
  const { display, counts, geneSet, displayKey, from } = chosenMatrices
  if ((counts != null) !== p.carried) {
    throw new BuildError('the plan and the build disagree on whether the raw counts travel: this is a bug in the converter')
  }
  const expr = ms.expression ?? measureExpression(from, displayKey, from === 'lognorm' ? rowsOf(display, sampleRuns(display.nCells)) : null)
  const expression = expr.label

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
  const cond = condCol?.kind === 'categorical' ? grouping(condCol, 'all cells', n) : null
  let sample: Grouping = grouping(sampleCol, 'all cells', n)

  if (sampleCol?.kind !== 'categorical' && cond) {
    // With no sample column the one sample is the whole object, and it holds
    // every condition: as one sample, every cell used to be given the first
    // cell's condition and the other groups came out empty. So the cells of
    // each condition are one sample. Every one is split, not only when it
    // spans, because a part of a split object can hold a single condition and
    // must still name its cells as the other parts do.
    sample = splitByCondition(sample, cond, true) ?? sample
    note(`no sample column — the cells of each condition are treated as one sample (${list(sample.levels)}). Composition cannot show between-animal spread, and pseudobulk needs several samples per group, so only the per-cell test will be available`)
  } else if (!sampleCol) {
    note('no sample column — every cell is treated as one sample. Composition cannot show between-animal spread, and pseudobulk needs several samples per group, so only the per-cell test will be available')
  }

  // One condition per sample: the bundle's model, so it is resolved here rather
  // than left for the studio to discover. Every sample now sits in one
  // condition, so its first cell speaks for all of them.
  let conditions: string[]
  let sampleCond: string[]
  if (cond) {
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
  // Said from its width, never from its name: a key holding "pca" is not a record of a PCA.
  if (emb.nDims > 2) {
    note(`${emb.key} has ${emb.nDims} columns and only its first two are drawn: the first two dimensions of a latent space, a much coarser picture than a UMAP`)
  }

  // The other embeddings a reader can switch to. An object routinely holds a
  // UMAP and a t-SNE of the same cells, and the lab once kept one and threw
  // the rest away — 8 bytes a cell to carry, against a matrix that costs 8
  // bytes a *stored value*. Two kinds are left out, and named. A latent space
  // (a PCA, a corrected PCA, anything wider than two columns) cut to its first
  // two columns is not a picture of the cells, whatever it is called. And a
  // byte copy of an embedding already offered is the same picture under a
  // second name: an integration stage writes its chosen method's UMAP twice,
  // as X_umap and X_umap_<method>. Copies are found on the whole object's
  // coordinates, so every part of a split decides the same way. The chosen one
  // stays first and stays in `embed.f32`, so a reader that knows nothing about
  // this sees no change at all.
  // Which travel is the plan's (plan.ts): copies were found there, on the whole
  // object's coordinates, before any part was cut.
  const embedFiles: EmbeddingEntry[] = [{ key: emb.key, file: 'embed.f32' }]
  const extraXY: { file: string; xy: Float32Array }[] = []
  const usedNames = new Set(['embed.f32'])
  const latent: string[] = []
  const copies: string[] = []
  for (const e of scan.embeddings) {
    if (e.key === emb.key || e.nDims < 2) continue
    if (e.nDims > 2) {
      latent.push(`${e.key} (${e.nDims} columns)`)
      continue
    }
    const twin = p.embeddings.copies.find(c => c.key === e.key)
    if (twin) {
      copies.push(`${e.key} (the same coordinates as ${twin.of})`)
      continue
    }
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
    let file = `embed.${safeEntry(e.key)}.f32`
    for (let i = 2; usedNames.has(file); i++) file = `embed.${safeEntry(e.key)}-${i}.f32`
    usedNames.add(file)
    extraXY.push({ file, xy: coords })
    embedFiles.push({ key: e.key, file })
  }
  if (extraXY.length) {
    note(`${embedFiles.length} embeddings carried (${embedFiles.map(e => e.key).join(', ')}); ${emb.key} opens by default`)
  }
  if (latent.length || copies.length) {
    const why = [
      latent.length ? `${list(latent)}, ${latent.length > 1 ? 'latent spaces' : 'a latent space'}: the first two columns of one are not a picture of the cells` : '',
      copies.length ? `${list(copies)}, ${copies.length > 1 ? 'copies' : 'a copy'} offered once` : '',
    ].filter(Boolean).join('; ')
    note(`not offered as embeddings: ${why}.${latent.length ? ' A latent space can still be chosen as the one that opens first, and is then drawn by its first two columns' : ''}`)
  }

  // What was measured on the whole object (plan.ts): the counts left out, the
  // cells not called and their coordinates. The same in every part.
  for (const s of ms.notes) note(s)

  // The run's identity, read from uns by the reader. Said in the notes too,
  // because the notes are what the studio shows.
  if (scan.run) {
    const paths = scan.run.records.map(r => r.path)
    note(paths.length
      ? `run carries the run's identity, not its record: from ${list(paths, 20)}, each record's tool, version, commit, stage, run_key and job where it holds them as strings of 256 characters or fewer, and nothing else`
      : "the object's uns holds no run record the lab can carry (none with a tool, version, commit, stage, run_key or job stored as a string of 256 characters or fewer), so this bundle cannot say which run made it")
  }

  // ---- QC ---------------------------------------------------------------
  onProgress('Computing QC')
  const MITO = ['percent_mito', 'pct_counts_mt', 'percent.mt', 'percent_mt']
  const mitoCol = qcColumn(scan, MITO)
  // The total and the genes detected, each with the words of what was summed,
  // from the one function that sums them.
  const perCell = perCellQC(scan, chosenMatrices, note)
  const total = perCell.total.values
  const nGene = perCell.genes.values
  let mito = mitoCol?.values ?? null
  const qcFacts = {
    total: perCell.total.fact,
    genes_detected: perCell.genes.fact,
    mito: mitoCol
      ? factOr(() => recorded(mitoCol.name, `obs/${mitoCol.name}`), note, `the mitochondrial column's name`)
      : notRecorded(MITO.map(k => `obs/${k}`).join(', ')),
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
    // Only the ones this part has cells for, like `clusters` itself.
    ...(choices.sentinels ? { sentinels: choices.sentinels.filter(s => cluster.levels.includes(s)) } : {}),
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
    // Walked from the table in facts.ts, never written by hand.
    facts: factsOf({
      recorded: scan.recorded,
      measured: {
        expression: expr.fact, counts: ms.counts, ...qcFacts,
        ...(ms.gene_ids ? { gene_ids: ms.gene_ids } : {}),
      },
      chosen: {
        cell_type_column: choiceOf(choices.cluster, choices.by?.cluster, note, 'the cell type column'),
        embedding: choiceOf(choices.embedding, choices.by?.embedding, note, 'the embedding'),
        // The columns provenance.sample and .condition name, and who chose them.
        sample_column: choiceOf(sampleCol?.kind === 'categorical' ? sampleCol.name : null, choices.by?.sample, note, 'the sample column'),
        condition_column: choiceOf(condCol?.kind === 'categorical' ? condCol.name : null, choices.by?.condition, note, 'the condition column'),
        ...(choices.sentinels
          ? { not_called: choiceOf(choices.sentinels.length ? choices.sentinels.join(', ') : null, choices.by?.sentinels, note, 'the labels of the cells not called') }
          : {}),
      },
    }),
    ...(scan.run ? { run: scan.run } : {}),
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
  const zip = zipSync(files, { level: 6, mtime: ZIP_MTIME })
  return { zip, meta, pseudobulkColumns: pbColumns }
}
