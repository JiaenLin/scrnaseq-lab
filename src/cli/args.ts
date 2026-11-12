// `scrnaseq-lab convert` arguments: parsed, never interpreted.
//
// An absent flag stays `undefined` here, which is the only thing that lets the
// converter tell "guess this" from "I said so". What a value means — whether
// `none` is a column or the absence of one — is decided later, against the
// object, because only the object knows what columns it has.

import { CliError, EXIT } from './errors.ts'

export interface ConvertArgs {
  input: string
  output: string
  cluster?: string
  sample?: string
  condition?: string
  /** Repeatable. `['none']` asks for no extra column. */
  extras?: string[]
  /** Repeatable `a,b[,c]`. An empty list inside means `--interaction none`. */
  interactions?: string[][]
  embedding?: string
  label?: string
  /** Repeatable. `['none']` marks no level as a sentinel. */
  sentinels?: string[]
  help?: boolean
  /**
   * The most stored values one part holds. Never parsed from the command line:
   * set only by a test or a proof script, to split a small object as a large
   * one is split (lib/plan.ts, lib/shard.ts).
   */
  splitAt?: number
}

const SINGLE = ['cluster', 'sample', 'condition', 'embedding', 'label'] as const
const REPEATED = ['extra', 'interaction', 'sentinel'] as const
type Single = typeof SINGLE[number]
type Repeated = typeof REPEATED[number]

const usage = (msg: string) => new CliError(msg, EXIT.usage)

export function parseConvertArgs(argv: string[]): ConvertArgs {
  const paths: string[] = []
  const single: Partial<Record<Single, string>> = {}
  const repeated: Partial<Record<Repeated, string[]>> = {}
  let help = false

  for (let i = 0; i < argv.length; i++) {
    const tok = argv[i]
    if (tok === '--help' || tok === '-h') { help = true; continue }
    if (tok === '--') { paths.push(...argv.slice(i + 1)); break }
    if (!tok.startsWith('--')) { paths.push(tok); continue }

    const eq = tok.indexOf('=')
    const name = tok.slice(2, eq < 0 ? undefined : eq)
    const isSingle = (SINGLE as readonly string[]).includes(name)
    const isRepeated = (REPEATED as readonly string[]).includes(name)
    if (!isSingle && !isRepeated) {
      throw usage(`unknown option --${name}. The options are ${[...SINGLE, ...REPEATED].map(f => `--${f}`).join(' ')} and --help`)
    }
    let value: string | undefined
    if (eq >= 0) value = tok.slice(eq + 1)
    else {
      const next = argv[i + 1]
      if (next !== undefined && !next.startsWith('--')) { value = next; i++ }
    }
    if (!value) throw usage(`--${name} needs a value (--${name}=VALUE when the value starts with --)`)

    if (isSingle) {
      const key = name as Single
      if (single[key] !== undefined) throw usage(`--${name} was given twice ("${single[key]}" and "${value}")`)
      single[key] = value
    } else {
      (repeated[name as Repeated] ??= []).push(value)
    }
  }

  if (help) return { input: '', output: '', help: true }
  if (paths.length < 2) throw usage('convert needs an input .h5ad and an output .zip: scrnaseq-lab convert IN.h5ad OUT.zip [options]')
  if (paths.length > 2) throw usage(`convert takes one input and one output, but was given ${paths.length} paths: ${paths.join(' ')}`)

  let interactions: string[][] | undefined
  if (repeated.interaction) {
    interactions = repeated.interaction.map(v => {
      if (v.trim() === 'none') return []
      const names = v.split(',').map(s => s.trim()).filter(Boolean)
      if (names.length < 2) throw usage(`--interaction ${v}: an interaction needs two or more columns, separated by commas`)
      if (new Set(names).size !== names.length) throw usage(`--interaction ${v} names a column twice`)
      return names
    })
    if (interactions.some(x => !x.length) && interactions.length > 1) {
      throw usage('--interaction none cannot be given beside an interaction')
    }
  }

  return {
    input: paths[0],
    output: paths[1],
    cluster: single.cluster,
    sample: single.sample,
    condition: single.condition,
    extras: repeated.extra,
    interactions,
    embedding: single.embedding,
    label: single.label,
    sentinels: repeated.sentinel,
  }
}
