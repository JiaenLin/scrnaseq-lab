// What a bundle states about its object, each statement with its source.
//
// One closed table, FACTS, and one type, Fact. `meta.facts` is built only by
// walking the table, and a Fact cannot be made without its source: a recorded
// one without the place it was read at (or, when the place holds nothing, where
// the converter looked), a measured one without the matrix it was measured on
// and how, a chosen one without who chose. The constructors below refuse the
// rest at run time, and the type refuses it at compile time
// (src/lib/facts.typecheck.ts).
//
// The kinds and their words are the Studio's (its src/lib/facts.ts and
// tools/BUNDLE.md, section "facts"): only a RECORDED fact with its place is
// worded "as recorded in the object". A fact is never decided from the look of
// a name: an embedding called X_harmony is not a record that harmony was run, a
// column called predicted_doublet is not a record of a doublet tool, an uns key
// called log1p_cp10k is not a record of a normalisation. Where the object
// records nothing at the place below, the fact is "not recorded", with that
// place, and nothing fills it.

/** Who chose: a flag of the command, the command's guess from a name, the lab's page. */
export type By = 'flag' | 'guess' | 'app'

export type Recorded = { kind: 'recorded'; value: string; at: string }
export type NotRecorded = { kind: 'recorded'; value: null; searched: string }
export type Measured = {
  kind: 'measured'; value: string | number | boolean | null; of: string; how: string
  /** Only on `counts`: whether the bundle carries them, and why not. */
  carried?: boolean; why?: string
}
export type Chosen = { kind: 'chosen'; value: string | null; by: By }
export type Fact = Recorded | NotRecorded | Measured | Chosen

export class FactError extends Error {}

/** The reader's limits (the Studio's src/lib/facts.ts): a value, a place, a sentence. */
export const MAX_VALUE = 64
export const MAX_PLACE = 256
export const MAX_HOW = 300

/** One short line: no line break or other control character, no markup. */
export function isLine(v: unknown, max: number): v is string {
  if (typeof v !== 'string' || v.length === 0 || v.length > max) return false
  for (let i = 0; i < v.length; i++) {
    const c = v.charCodeAt(i)
    if (c < 32 || c === 127 || c === 60 || c === 62) return false
  }
  return true
}

const need = (v: unknown, max: number, what: string): string => {
  if (!isLine(v, max)) throw new FactError(`${what} is not one line of ${max} characters or fewer`)
  return v
}

/** Read at a named place in the object, verbatim. */
export function recorded(value: string, at: string): Recorded {
  return { kind: 'recorded', value: need(value, MAX_VALUE, 'a recorded value'), at: need(at, MAX_PLACE, 'its place') }
}

/** Nothing at the place: said with where the converter looked. */
export function notRecorded(searched: string): NotRecorded {
  return { kind: 'recorded', value: null, searched: need(searched, MAX_PLACE, 'where the converter looked') }
}

/** Measured from the values: of what, and how. */
export function measured(
  value: string | number | boolean | null, of: string, how: string,
  counts?: { carried: boolean; why?: string },
): Measured {
  const ok = value === null || typeof value === 'boolean' || (typeof value === 'number' && Number.isFinite(value)) ||
    isLine(value, MAX_VALUE)
  if (!ok) throw new FactError('a measured value is not one short line, a number, or yes or no')
  const f: Measured = { kind: 'measured', value, of: need(of, MAX_PLACE, 'what was measured'), how: need(how, MAX_HOW, 'how') }
  if (counts) {
    if (typeof counts.carried !== 'boolean') throw new FactError('carried is not yes or no')
    f.carried = counts.carried
    if (!counts.carried) f.why = need(counts.why, MAX_HOW, 'why they are not carried')
  }
  return f
}

/** A decision of the conversion, and who made it. */
export function chosen(value: string | null, by: By): Chosen {
  if (by !== 'flag' && by !== 'guess' && by !== 'app') throw new FactError('a choice does not say who chose')
  return { kind: 'chosen', value: value === null ? null : need(value, MAX_PLACE, 'a chosen value'), by }
}

// ---------------------------------------------------------------------------
// The closed table.

/** The ids the Studio's reader takes (its src/lib/facts.ts, FactId). */
export type FactId =
  | 'integration' | 'integration_embedding' | 'normalization' | 'doublets' | 'ambient'
  | 'expression' | 'counts' | 'total' | 'genes_detected' | 'mito'
  | 'cell_type_column' | 'embedding' | 'not_called' | 'gene_ids'
  // The lab's own two, which the Studio's reader (at db78d6e) has no id for:
  // it ignores an id it does not know (its tools/BUNDLE.md, "facts"), so they
  // ride in meta.facts unread there until a reader takes them.
  | 'sample_column' | 'condition_column'

/** The facts read from the object's own records, by the .h5ad reader (h5ad.ts, readRecorded). */
export type RecordedId = 'integration' | 'integration_embedding' | 'normalization' | 'doublets' | 'ambient'

/** Everything a conversion knows about its object, each already a Fact or absent. */
export interface Sources {
  /** Read from the object's records. Absent for a Seurat object: no place is read, so none is stated. */
  recorded?: Partial<Record<RecordedId, Recorded | NotRecorded>>
  /** Measured or read by the build and by the plan (build.ts, plan.ts). */
  measured: Partial<Record<'expression' | 'counts' | 'total' | 'genes_detected' | 'mito' | 'gene_ids', Fact>>
  /** The conversion's choices, with who made them; absent when the caller did not say. */
  chosen: Partial<Record<'cell_type_column' | 'embedding' | 'not_called' | 'sample_column' | 'condition_column', Chosen>>
}

interface Row {
  /** The kinds this fact may have. Any other is refused: a method is never "measured". */
  kinds: readonly Fact['kind'][]
  /** Where it comes from. */
  read: (s: Sources) => Fact | undefined
}

const rec = (id: RecordedId): Row => ({ kinds: ['recorded'], read: s => s.recorded?.[id] })

export const FACTS: Record<FactId, Row> = {
  // A dict directly under uns holding default_method (and default_embedding beside it).
  integration: rec('integration'),
  integration_embedding: rec('integration_embedding'),
  // scanpy's record of log1p, uns/log1p.
  normalization: rec('normalization'),
  // A dict directly under uns holding doublet_method; one holding ambient_method.
  doublets: rec('doublets'),
  ambient: rec('ambient'),
  expression: { kinds: ['measured'], read: s => s.measured.expression },
  counts: { kinds: ['measured'], read: s => s.measured.counts },
  total: { kinds: ['recorded', 'measured'], read: s => s.measured.total },
  genes_detected: { kinds: ['recorded', 'measured'], read: s => s.measured.genes_detected },
  mito: { kinds: ['recorded'], read: s => s.measured.mito },
  cell_type_column: { kinds: ['chosen'], read: s => s.chosen.cell_type_column },
  embedding: { kinds: ['chosen'], read: s => s.chosen.embedding },
  not_called: { kinds: ['chosen'], read: s => s.chosen.not_called },
  gene_ids: { kinds: ['measured'], read: s => s.measured.gene_ids },
  // The columns the samples and the conditions were read from, and who chose
  // them: written where the Studio's reader ignores what it does not know.
  sample_column: { kinds: ['chosen'], read: s => s.chosen.sample_column },
  condition_column: { kinds: ['chosen'], read: s => s.chosen.condition_column },
}
export const FACT_IDS = Object.keys(FACTS) as FactId[]

/** Where the converter looks for each recorded fact: what a "not recorded" names. */
export const SEARCHED: Record<RecordedId, string> = {
  integration: 'a dict directly under uns holding default_method',
  integration_embedding: 'a dict directly under uns holding default_embedding',
  normalization: 'uns/log1p',
  doublets: 'a dict directly under uns holding doublet_method',
  ambient: 'a dict directly under uns holding ambient_method',
}

/** The key each recorded fact is read by, in a dict directly under uns. */
export const RECORD_KEYS = {
  integration: 'default_method',
  integration_embedding: 'default_embedding',
  doublets: 'doublet_method',
  ambient: 'ambient_method',
} as const

export type MetaFacts = { v: 1 } & Partial<Record<FactId, Fact>>

/** `meta.facts`: the table walked, in its order; a fact of a kind its row never has is a bug, and stops. */
export function factsOf(s: Sources): MetaFacts {
  const out: MetaFacts = { v: 1 }
  for (const id of FACT_IDS) {
    const f = FACTS[id].read(s)
    if (!f) continue
    if (!FACTS[id].kinds.includes(f.kind)) throw new FactError(`${id} is never ${f.kind}`)
    out[id] = f
  }
  return out
}

/**
 * The four older provenance fields, which the deployed Studio of before the
 * facts reads alone: each the RECORDED value, or null. Never a guess. A reader
 * that knows `facts` ignores them where a fact is given.
 */
export function legacyProvenance(r: Sources['recorded']): Record<'normalization' | 'integration' | 'doublets' | 'ambient', string | null> {
  const v = (id: RecordedId): string | null => {
    const f = r?.[id]
    return f && f.value !== null ? f.value : null
  }
  return { normalization: v('normalization'), integration: v('integration'), doublets: v('doublets'), ambient: v('ambient') }
}

/** A matrix's key as a place in the object: `layers["counts"]` is `layers/counts`, `raw.X` is `raw/X`. */
export function placeOf(key: string): string {
  const layer = /^layers\["(.*)"\]$/.exec(key)
  if (layer) return `layers/${layer[1]}`
  if (key === 'raw.X') return 'raw/X'
  return key
}
