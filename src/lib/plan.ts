// One plan for the page and the build.
//
// Every sentence the lab's page shows before converting (whether the object is
// held in parts, what travels, whether raw counts will be in the bundle,
// whether pseudobulk will be possible, how many embeddings) is computed here,
// by `plan`, and the command and the worker run the same plan: which matrix is
// read, which embeddings are carried and which are not, whether the counts
// travel. The page used to predict from the scan apart from the build, and
// promised what the build then did not do ("Pseudobulk DESeq2: yes" on an
// object whose counts a split leaves out; "All 11 travel" where 5 did).
//
// `plan` reads only the scan's description (what the page has, ScanInfo);
// `measureObject` measures the facts that need the values (the per-cell totals
// behind the expression, the counts, the gene names, the cells with no
// coordinates), on the WHOLE object before any part is cut from it, so that
// every part of a split states them as the whole conversion does.

import {
  FactError, chosen, measured, placeOf,
  type By, type Chosen, type Fact, type Measured,
} from './facts.ts'
import type { CellMajor, MatrixKind } from './matrix.ts'
import {
  geneAccessions, replicatesPerCondition,
  type Choices, type Column, type EmbeddingRef, type Scan, type ScanInfo,
} from './scan.ts'
import { COMFORTABLE_NNZ, chooseSplit, type SplitDecision } from './shard.ts'

/** Samples a group needs before pseudobulk DESeq2 is offered: the Studio's rule. */
export const MIN_REPS_PB = 4

/**
 * Below this share of the gene list, a matrix is not the object's expression.
 * See usableCoverage in build.ts, which says so in the notes.
 */
export const COVERAGE_FLOOR = 0.25

/** What the plan reads: the scan's description, as the page has it. */
export interface PlanInfo {
  nCells: number
  columns: Column[]
  matrices: { key: string; geneSet: string; kind: MatrixKind; covered?: number; streamable?: boolean }[]
  /** The matrices the scan could not read: what they hold is not known, so nothing is concluded of them. */
  unreadable?: { key: string; why: string }[]
  /** `same`: the first embedding, in the object's order, whose two columns are these bit for bit. */
  embeddings: { key: string; nDims: number; same?: string }[]
  nnzPerCell?: Int32Array | null
}

/** Which matrix the expression comes from, and which is summed for pseudobulk. */
export interface MatrixPick {
  display: { key: string; geneSet: string; from: 'lognorm' | 'counts' | 'log-counts' | 'linear' } | null
  counts: { key: string; kind: 'counts' | 'log-counts' } | null
  /** Why the object's counts are not carried; null when they are, or when it holds none. */
  dropped: string | null
}

export interface Plan {
  split: SplitDecision | null
  /** The matrix read in parts when the object is split. */
  stream: string | null
  pick: MatrixPick
  /** The object's raw counts, measured on its matrices' sampled values: the first matrix that holds them. */
  countsIn: { key: string; kind: 'counts' | 'log-counts' } | null
  /**
   * With no counts among the matrices read: the places of those that could not
   * be read. Then whether the object holds raw counts is not known, and nothing
   * is said of it but that (a failed read states no fact).
   */
  countsUnread: string[] | null
  carried: boolean
  embeddings: {
    /** In the bundle, the one that opens first. */
    carried: string[]
    latent: { key: string; nDims: number }[]
    copies: { key: string; of: string }[]
  }
  /** The chosen embedding is wider than two columns: only its first two are drawn. */
  wide: { key: string; nDims: number } | null
  pseudobulk: { possible: boolean; because: string }
}

/** Matrices that reach most of the gene list, or every one when none can be costed. The rule of usableCoverage. */
export function usable<T extends { covered?: number }>(ms: T[]): T[] {
  const cov = ms.filter(m => m.covered !== undefined).map(m => m.covered as number)
  if (cov.length < 2) return ms
  const best = Math.max(...cov)
  const out = ms.filter(m => m.covered === undefined || m.covered >= COVERAGE_FLOOR || m.covered * 2 > best)
  return out.length ? out : ms
}

/** The preference of chooseMatrices, from the classifications alone. */
export function pickMatrices(kinds: { key: string; kind: MatrixKind; geneSet: string }[]): MatrixPick {
  const first = (k: MatrixKind) => kinds.find(m => m.kind === k)
  const c = first('counts') ?? first('log-counts')
  const ln = first('lognorm')
  const counts = c ? { key: c.key, kind: c.kind as 'counts' | 'log-counts' } : null
  if (ln) {
    if (c && c.geneSet !== ln.geneSet) {
      return { display: { key: ln.key, geneSet: ln.geneSet, from: 'lognorm' }, counts: null,
        dropped: 'they index a different gene list from the expression matrix' }
    }
    return { display: { key: ln.key, geneSet: ln.geneSet, from: 'lognorm' }, counts, dropped: null }
  }
  if (c) return { display: { key: c.key, geneSet: c.geneSet, from: c.kind as 'counts' | 'log-counts' }, counts, dropped: null }
  const lin = first('linear')
  if (lin) return { display: { key: lin.key, geneSet: lin.geneSet, from: 'linear' }, counts: null, dropped: null }
  return { display: null, counts: null, dropped: null }
}

/** The matrix a split reads in parts: the preference of streamableMatrix. */
export function streamPick<T extends { kind: MatrixKind; streamable?: boolean; covered?: number }>(ms: T[]): T | null {
  const ok = usable(ms).filter(m => m.streamable && m.kind !== 'scaled')
  for (const kind of ['lognorm', 'counts', 'log-counts', 'linear'] as const) {
    const hit = ok.find(m => m.kind === kind)
    if (hit) return hit
  }
  return null
}

/** A column by its name, or none. */
const column = (o: { columns: Column[] }, name: string | null | undefined): Column | undefined =>
  (name ? o.columns.find(c => c.name === name) : undefined)

/** "a", "a and b", "a, b and c". */
const and = (xs: string[]) => (xs.length < 2 ? xs.join('') : `${xs.slice(0, -1).join(', ')} and ${xs[xs.length - 1]}`)

export function plan(info: PlanInfo, choices: Choices, comfortable = COMFORTABLE_NNZ): Plan {
  const split = chooseSplit(info.columns, info.nnzPerCell ?? null, info.nCells, comfortable)
  const kinds = usable(info.matrices)
  const whole = pickMatrices(kinds)
  const c = kinds.find(m => m.kind === 'counts') ?? kinds.find(m => m.kind === 'log-counts')
  const countsIn = c ? { key: c.key, kind: c.kind as 'counts' | 'log-counts' } : null
  const lost = info.unreadable ?? []
  const countsUnread = !countsIn && lost.length ? lost.map(x => placeOf(x.key)) : null

  let pick = whole
  let stream: string | null = null
  if (split) {
    const s = streamPick(info.matrices)
    stream = s ? s.key : null
    pick = s ? pickMatrices([s]) : { display: null, counts: null, dropped: null }
    if (countsIn && !pick.counts) {
      pick = { ...pick, dropped: `the object is converted in ${split.plan.shards.length} parts, and one matrix is read for them, ${s ? placeOf(s.key) : 'none'}, the expression` }
    }
  }
  const carried = pick.counts !== null

  // ---- embeddings: the default first, then every two-column one once ------
  const def = info.embeddings.find(e => e.key === choices.embedding)
  const classOf = (key: string) => info.embeddings.find(e => e.key === key)?.same ?? key
  const carriedKeys = [choices.embedding]
  const latent: { key: string; nDims: number }[] = []
  const copies: { key: string; of: string }[] = []
  for (const e of info.embeddings) {
    if (e.key === choices.embedding || e.nDims < 2) continue
    if (e.nDims > 2) { latent.push({ key: e.key, nDims: e.nDims }); continue }
    const twin = carriedKeys.find(k => classOf(k) === classOf(e.key))
    if (twin) copies.push({ key: e.key, of: twin })
    else carriedKeys.push(e.key)
  }

  // ---- pseudobulk: the samples a group, and the counts carried -------------
  const sampleCol = column(info, choices.sample)
  const condCol = column(info, choices.condition)
  const reps = Object.values(replicatesPerCondition(sampleCol, condCol))
  const enough = condCol != null && reps.length > 1 && reps.every(v => v >= MIN_REPS_PB)
  const because = countsUnread ? `${and(countsUnread)} could not be read, so whether this object holds raw counts is not known, and none are summed`
    : !countsIn ? 'this object holds no raw counts to sum'
    : !carried ? `the raw counts in ${placeOf(countsIn.key)} will not be in the bundle: ${pick.dropped}`
    : !condCol ? 'needs a condition to compare'
    : reps.length < 2 ? 'needs at least two groups'
    : enough ? `${reps.join(' v ')} samples per group`
    : `${reps.join(' v ')} samples per group — fewer than ${MIN_REPS_PB}, so only the per-cell Wilcoxon is offered`

  return {
    split, stream, pick, countsIn, countsUnread, carried,
    embeddings: { carried: carriedKeys, latent, copies },
    wide: def && def.nDims > 2 ? { key: def.key, nDims: def.nDims } : null,
    pseudobulk: { possible: enough && carried, because },
  }
}

/** What the page says before converting: each sentence from the plan the build runs. */
export interface PageSays {
  parts: string | null
  travel: string
  notOffered: string | null
  wide: string | null
  counts: string
  pseudobulk: { yes: boolean; because: string }
}

export function pageSays(p: Plan): PageSays {
  const e = p.embeddings
  const n = e.carried.length
  const not = [
    ...e.latent.map(x => `${x.key} (${x.nDims} columns, a latent space)`),
    ...e.copies.map(x => `${x.key} (the same coordinates as ${x.of})`),
  ]
  return {
    parts: p.split ? `${p.split.reason}.` : null,
    travel: `${n} embedding${n === 1 ? '' : 's'} travel${n === 1 ? 's' : ''} in the file: ${e.carried.join(', ')}; ${e.carried[0]} opens first.`,
    notOffered: not.length ? `Not offered beside ${n === 1 ? 'it' : 'them'}: ${and(not)}.` : null,
    wide: p.wide ? `${p.wide.key} has ${p.wide.nDims} columns, and only its first two are drawn: a much coarser picture than a UMAP.` : null,
    counts: p.carried && p.pick.counts ? `Raw counts travel: ${placeOf(p.pick.counts.key)}, summed per sample and cell type.`
      : p.countsIn ? `Your object holds raw counts at ${placeOf(p.countsIn.key)}; they will not be in the bundle: ${p.pick.dropped}.`
      : p.countsUnread ? `${and(p.countsUnread)} could not be read, so whether this object holds raw counts is not known; the bundle carries none.`
      : 'This object holds no raw counts: no matrix’s sampled values are integers.',
    pseudobulk: { yes: p.pseudobulk.possible, because: p.pseudobulk.because },
  }
}

// ---------------------------------------------------------------------------
// What the plan reads, from a scan.

/** The same float32 values, bit for bit — what the bundle would store twice. */
export function sameCoordinates(a: Float32Array, b: Float32Array): boolean {
  if (a.length !== b.length) return false
  const x = new Uint32Array(a.buffer, a.byteOffset, a.length)
  const y = new Uint32Array(b.buffer, b.byteOffset, b.length)
  for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) return false
  return true
}

/** Each embedding's class of identical coordinates, on the whole object. */
function sameAs(es: EmbeddingRef[]): Map<string, string> {
  const reps: { key: string; at: Float32Array }[] = []
  const out = new Map<string, string>()
  for (const e of es) {
    if (e.nDims < 2) continue
    let at: Float32Array
    try {
      at = e.whole ? e.whole() : e.load()
    } catch {
      continue
    }
    const r = reps.find(x => sameCoordinates(x.at, at))
    out.set(e.key, r ? r.key : e.key)
    if (!r) reps.push({ key: e.key, at })
  }
  return out
}

/** The scan described for the plan; the worker sends this to the page. */
export function planInfo(scan: Scan, nnzPerCell: Int32Array | null): Pick<ScanInfo, 'matrices' | 'embeddings'> & { nnzPerCell: Int32Array | null; unreadable: { key: string; why: string }[] } {
  const same = sameAs(scan.embeddings)
  return {
    nnzPerCell,
    unreadable: scan.unreadableMatrices ?? [],
    matrices: scan.matrices.map(m => {
      let covered: number | undefined
      try { covered = m.coverage ? m.coverage() : undefined } catch { covered = undefined }
      return {
        key: m.key, geneSet: m.geneSet, nGenes: m.nGenes, nCells: m.nCells, kind: m.kind,
        streamable: !!m.loadGroups, ...(covered === undefined ? {} : { covered }),
      }
    }),
    embeddings: scan.embeddings.map(e => ({ key: e.key, nDims: e.nDims, ...(same.has(e.key) ? { same: same.get(e.key) } : {}) })),
  }
}

// ---------------------------------------------------------------------------
// The facts measured on the whole object.

/** Cells read for the per-cell totals: up to 2,000, in 8 runs spread over the object. */
export function sampleRuns(nCells: number, want = 2000, runs = 8): [number, number][] {
  if (nCells <= want) return [[0, nCells]]
  const per = Math.floor(want / runs)
  const out: [number, number][] = []
  for (let r = 0; r < runs; r++) {
    const at = Math.min(nCells - per, Math.floor((nCells / runs) * r))
    out.push([at, at + per])
  }
  return out
}

/** The rows of the runs, from a matrix already in hand. */
export function rowsOf(m: CellMajor, runs: [number, number][]): CellMajor {
  let nnz = 0, n = 0
  for (const [a, b] of runs) { nnz += m.indptr[b] - m.indptr[a]; n += b - a }
  const indptr = new Int32Array(n + 1), indices = new Int32Array(nnz), data = new Float32Array(nnz)
  let c = 0, at = 0
  for (const [a, b] of runs) {
    for (let r = a; r < b; r++) {
      indptr[c++] = at
      for (let k = m.indptr[r]; k < m.indptr[r + 1]; k++) { indices[at] = m.indices[k]; data[at++] = m.data[k] }
    }
  }
  indptr[n] = at
  return { nCells: n, nGenes: m.nGenes, indptr, indices, data }
}

/** The expression's label for readers that know no facts (`meta.expression`), and its fact. */
export interface Expression { label: string; fact: Measured }

/**
 * What the expression values are, measured. On a log-normalized matrix, expm1
 * summed per cell over the sampled rows: "scaled to N per cell" only when every
 * sum is within 0.1% of N. Built here from counts or from linear values, it is
 * what this converter did to them.
 */
export function measureExpression(from: NonNullable<MatrixPick['display']>['from'], key: string, rows: CellMajor | null): Expression {
  const of = placeOf(key)
  if (from === 'counts' || from === 'log-counts') {
    return { label: 'log1p(CP10K)', fact: measured('log1p of counts scaled to 10,000 per cell', of,
      `computed at conversion from the counts in ${of}${from === 'log-counts' ? ', recovered by expm1 of its stored log1p(counts)' : ''}`) }
  }
  if (from === 'linear') {
    return { label: 'log1p(normalized)', fact: measured('log1p of the stored normalized values', of,
      'its sampled stored values are non-integers above the log scale; log1p applied at conversion') }
  }
  if (!rows || rows.nCells === 0) {
    return { label: 'log-normalized', fact: measured('log-normalized; the per-cell total was not measured', of,
      'the matrix is not stored by cell, so its per-cell sums were not read') }
  }
  const sums: number[] = []
  for (let c = 0; c < rows.nCells; c++) {
    let s = 0
    for (let k = rows.indptr[c]; k < rows.indptr[c + 1]; k++) s += Math.expm1(rows.data[k])
    sums.push(s)
  }
  const sorted = [...sums].sort((a, b) => a - b)
  const median = sorted[Math.floor(sorted.length / 2)]
  const within = (n: number) => n > 0 && sums.every(s => Math.abs(s - n) <= 1e-3 * n)
  const nice = Number(median.toPrecision(2))
  const N = within(nice) ? nice : within(Math.round(median)) ? Math.round(median) : null
  const span = `${sums.length.toLocaleString('en-US')} cells in ${rows.nCells < 2000 ? 'all' : '8 runs spread over the object'}`
  const lo = sorted[0], hi = sorted[sorted.length - 1]
  const r = (x: number) => x.toLocaleString('en-US', { maximumFractionDigits: 3 })
  if (N === null) {
    return { label: 'log-normalized', fact: measured('log-normalized, no constant per-cell total', of,
      `expm1 summed per cell over ${span}: the sums run from ${r(lo)} to ${r(hi)}`) }
  }
  const label = N === 10_000 ? 'log1p(CP10K)' : N === 1_000_000 ? 'log1p(CPM)' : `log1p(counts scaled to ${N})`
  return { label, fact: measured(`log1p of counts scaled to ${N.toLocaleString('en-US')} per cell`, of,
    `expm1 summed per cell over ${span}: every sum is within 0.1% of ${N.toLocaleString('en-US')}`) }
}

/** The facts the plan's object holds, measured on the whole of it, and the notes that say what was measured. */
export interface Measurements {
  expression?: Expression
  counts: Measured
  gene_ids?: Measured
  notes: string[]
}

export function measureObject(scan: Scan, p: Plan, choices: Choices): Measurements {
  const notes: string[] = []
  const d = p.pick.display
  const m = d ? scan.matrices.find(x => x.key === d.key) : undefined

  // The expression, from whole rows read out of the file where they can be.
  let expression: Expression | undefined
  if (d && (d.from !== 'lognorm' || m?.readRuns)) {
    expression = measureExpression(d.from, d.key, d.from === 'lognorm' ? m!.readRuns!(sampleRuns(scan.nCells)) : null)
  }

  // The counts, said of the object, and whether the bundle carries them.
  const all = scan.matrices.map(x => placeOf(x.key))
  const ci = p.countsIn
  const counts = ci
    ? measured(true, placeOf(ci.key), ci.kind === 'counts'
      ? 'its sampled stored values (up to 20,000, in 8 blocks spread over it) are integers'
      : 'its sampled stored values (up to 20,000, in 8 blocks spread over it) are log1p of integers',
      p.carried ? { carried: true } : { carried: false, why: p.pick.dropped ?? 'they are not read' })
    : p.countsUnread
      // A failed read states no fact: the value is not known, and the words say why.
      ? measured(null, p.countsUnread.join(', ').slice(0, 256),
        `could not be read, so whether the object holds raw counts is not known; no matrix that could be read has sampled stored values that are integers or log1p of integers`,
        { carried: false, why: 'what the object holds there could not be read' })
      : measured(false, all.length === 1 ? all[0] : `every matrix of the object (${all.join(', ')})`.slice(0, 256),
        `no matrix's sampled stored values (up to 20,000 a matrix, in 8 blocks) are integers, nor log1p of integers`)
  if (ci && !p.carried) {
    notes.push(`the object holds raw counts at ${placeOf(ci.key)}, and this bundle does not carry them: ${p.pick.dropped}. Pseudobulk sums cannot be saved from it`)
  }

  // The gene names, every one of them.
  let gene_ids: Measured | undefined
  const genes = d ? scan.geneSets[d.geneSet] : undefined
  if (genes) {
    const acc = geneAccessions(genes)
    const n = (x: number) => x.toLocaleString('en-US')
    gene_ids = measured(`${n(acc)} of ${n(genes.length)} names are accessions`,
      d!.geneSet === 'raw' ? 'raw/var' : 'var',
      'every gene name tested against the accession forms (Ensembl, FlyBase, WormBase, TAIR, LOC and numeric ids); the rest are taken as symbols')
  }

  // The cells not called, and which have no coordinates: counted on every embedding carried.
  const sentinels = choices.sentinels ?? []
  const cl = column(scan, choices.cluster)
  if (sentinels.length && cl?.kind === 'categorical') {
    const codes = new Set(sentinels.map(s => cl.levels.indexOf(s)).filter(i => i >= 0))
    let total = 0
    for (const c of cl.codes) if (codes.has(c)) total++
    const missing: string[] = []
    const fine: string[] = []
    for (const key of p.embeddings.carried) {
      const e = scan.embeddings.find(x => x.key === key)
      if (!e) continue
      let xy: Float32Array
      try { xy = e.whole ? e.whole() : e.load() } catch { continue }
      let k = 0
      for (let i = 0; i < cl.codes.length; i++) {
        if (codes.has(cl.codes[i]) && !(Number.isFinite(xy[2 * i]) && Number.isFinite(xy[2 * i + 1]))) k++
      }
      if (k) missing.push(`${k.toLocaleString('en-US')} have no finite coordinates on ${key}`)
      else fine.push(key)
    }
    const named = sentinels.length < 2 ? sentinels.join('') : `${sentinels.slice(0, -1).join(', ')} or ${sentinels[sentinels.length - 1]}`
    notes.push(missing.length
      ? `of the ${total.toLocaleString('en-US')} cells labelled ${named}, measured on each embedding carried: ${missing.join('; ')}${fine.length ? `; every one has them on ${and(fine)}` : ''}`
      : `the ${total.toLocaleString('en-US')} cells labelled ${named} have finite coordinates on every embedding carried (${and(fine)})`)
  }
  return { ...(expression ? { expression } : {}), counts, ...(gene_ids ? { gene_ids } : {}), notes }
}

/** A chosen fact, or none when the caller did not say who chose, or the value is not one line. */
export function choiceOf(value: string | null, by: By | undefined, note: (s: string) => void, what: string): Chosen | undefined {
  if (!by) return undefined
  try {
    return chosen(value, by)
  } catch (e) {
    if (e instanceof FactError) { note(`${what} is not one line of text a reader can show, so the bundle states no fact for it`); return undefined }
    throw e
  }
}

export type { Fact }
