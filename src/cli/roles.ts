// Who plays which role, from the flags given and guesses for the rest.
//
// The rule is the one a pipeline needs: a flag that was given is used as given
// — checked against the object, never second-guessed — and a flag that was not
// given is guessed and the guess is said out loud, as the exact flag that
// would have produced it. Run with every guess pasted back in, the command
// writes the same bytes.
//
// What is guessed:
//   --cluster      the app's own pick: a named annotation, else the finest grouping
//   --sample       the column named like a sample
//   --condition    a design factor (see design.ts), the one named like a
//                  condition first; without a sample column, the name alone
//   --extra        every other design factor
//   --interaction  the two carried factors' interaction, when they cross
//   --embedding    a UMAP, then t-SNE, a force layout, PCA
//   --sentinel     the annotator's sentinel labels, when the annotation has them
//   --label        the file name
//
// `none` is a value for --sample, --condition, --extra, --interaction and
// --sentinel: "there is none, do not guess". A column that really is called
// "none" wins over the keyword.

import {
  candidates, crossCheck, guessCluster, guessEmbedding, guessRole, namedLike, unusable,
  type Choices, type Column, type Role, type Scan,
} from '../lib/scan.ts'
import { crossed, designFactors, interactionColumn } from '../lib/design.ts'
import type { ConvertArgs } from './args.ts'
import { CliError, EXIT } from './errors.ts'

/**
 * The single-cell profile's sentinel labels: what the annotator writes for a
 * cell it excluded, and for one it could not resolve. They label cells, not
 * cell types.
 */
export const PROFILE_SENTINELS = ['EXCLUDED', 'UNRESOLVED']

export interface Resolved {
  choices: Choices
  /** Columns the conversion adds to the object: one per interaction. */
  added: Column[]
  /** The interactions, as the column names they combine. */
  interactions: string[][]
  sentinels: string[]
  /** The flags that were guessed, in the order they were decided. */
  guessed: string[]
  /** One line per guessed flag, written as the flag itself. */
  said: string[]
  /** Sentences for the bundle's notes. */
  notes: string[]
  /** Said on the terminal only. */
  info: string[]
}

/** A value as it would have to be typed into a shell. */
const shell = (v: string) => (/^[\w.,:@%+=/-]+$/.test(v) ? v : `'${v.replace(/'/g, `'\\''`)}'`)

const ROLE_WORD: Record<Role, string> = {
  cluster: 'the cell annotation', sample: 'the sample column', condition: 'the condition', extra: 'an extra column',
}

export function resolveRoles(scan: Scan, a: ConvertArgs): Resolved {
  const n = scan.nCells
  const guessed: string[] = []
  const said: string[] = []
  const notes: string[] = []
  const info: string[] = []
  const has = (name: string) => scan.columns.some(c => c.name === name)
  const isNone = (v: string | undefined) => v === 'none' && !has('none')
  const col = (name: string | null): Column | undefined =>
    name == null ? undefined : scan.columns.find(c => c.name === name)

  const say = (flag: string, values: string[], why: string) => {
    guessed.push(flag)
    const shown = values.length ? values.map(v => `--${flag} ${shell(v)}`).join(' ') : `--${flag} none`
    said.push(`guessed ${shown}  (${why})`)
  }
  const named = (flag: string, name: string, role: Role, taken: (string | null)[] = []): Column => {
    const c = col(name)
    if (!c) {
      throw new CliError(`--${flag} ${name}: this object has no obs column called "${name}". Its columns are: ${scan.columns.map(x => x.name).join(', ') || '(none)'}`, EXIT.roles)
    }
    const why = unusable(c, role, n)
    if (why) throw new CliError(`--${flag} ${name}: "${name}" cannot be ${ROLE_WORD[role]} — ${why}`, EXIT.roles)
    if (taken.includes(c.name)) throw new CliError(`--${flag} ${name}: "${name}" already plays another role`, EXIT.roles)
    return c
  }

  // ---- the cell annotation ----------------------------------------------
  let cluster: string
  if (a.cluster !== undefined) {
    cluster = named('cluster', a.cluster, 'cluster').name
  } else {
    const g = guessCluster(scan)
    if (!g) throw new CliError('no column can be the cell annotation: none is a grouping of 2 to 300 levels. Name one with --cluster', EXIT.build)
    cluster = g
    say('cluster', [g], namedLike(g, 'cluster') ? 'named like a cell annotation' : 'the finest grouping; none is named like an annotation')
  }

  // ---- samples ------------------------------------------------------------
  let sample: string | null
  if (a.sample !== undefined) {
    sample = isNone(a.sample) ? null : named('sample', a.sample, 'sample', [cluster]).name
  } else {
    sample = guessRole(scan, 'sample', [cluster])
    say('sample', sample ? [sample] : [], sample ? 'named like a sample' : 'no column is named like a sample')
  }
  const sampleCol = col(sample)

  // ---- the design ---------------------------------------------------------
  const factors = designFactors(scan.columns, sampleCol, n, [cluster])
  const isFactor = new Set(factors.map(f => f.name))

  let condition: string | null
  if (a.condition !== undefined) {
    condition = isNone(a.condition) ? null : named('condition', a.condition, 'condition', [cluster, sample]).name
  } else if (sampleCol) {
    // `candidates` ranks by name first, then the coarsest — and keeps obs order
    // among equals — so among the factors this is the conventional choice.
    const pick = candidates(scan, 'condition', [cluster, sample]).find(c => isFactor.has(c.name)) ?? null
    condition = pick?.name ?? null
    say('condition', pick ? [pick.name] : [], !pick ? 'no column takes one value in every sample'
      : namedLike(pick.name, 'condition') ? 'named like a condition, and one value in every sample'
      : 'one value in every sample; no such column is named like a condition')
  } else {
    condition = guessRole(scan, 'condition', [cluster])
    say('condition', condition ? [condition] : [], condition
      ? 'named like a condition; with no sample column, nothing tells a factor from a cell label'
      : 'no sample column, and no column named like a condition')
  }

  let extras: string[]
  if (a.extras !== undefined) {
    extras = a.extras.length === 1 && isNone(a.extras[0]) ? []
      : [...new Set(a.extras.map(x => named('extra', x, 'extra').name))]
        .filter(x => x !== cluster && x !== sample && x !== condition)
  } else {
    extras = factors.map(f => f.name).filter(x => x !== condition)
    say('extra', extras, extras.length ? 'every other column that takes one value in every sample' : 'no other column takes one value in every sample')
  }

  // ---- interactions ---------------------------------------------------------
  let interactions: string[][]
  if (a.interactions !== undefined) {
    interactions = a.interactions.filter(x => x.length)
      .map(list => list.map(x => named('interaction', x, 'extra').name))
  } else {
    // In the object's own column order, so the name does not depend on which
    // factor happened to be chosen as the condition.
    const carried = scan.columns.filter(c => isFactor.has(c.name) && (c.name === condition || extras.includes(c.name)))
    if (carried.length === 2 && crossed(carried)) {
      interactions = [carried.map(c => c.name)]
      say('interaction', [carried.map(c => c.name).join(',')], 'the two carried factors cross')
    } else {
      interactions = []
      say('interaction', [], carried.length < 2 ? 'fewer than two design factors are carried'
        : carried.length > 2 ? `${carried.length} design factors are carried; name the interaction wanted`
        : 'the two carried factors are nested, not crossed')
    }
  }
  const added: Column[] = []
  for (const list of interactions) {
    const inter = interactionColumn(list.map(x => col(x)!))
    if (has(inter.name) || added.some(c => c.name === inter.name)) {
      throw new CliError(`--interaction ${list.join(',')}: the object already has a column called "${inter.name}"`, EXIT.roles)
    }
    if (inter.levels.length > 65535) {
      throw new CliError(`--interaction ${list.join(',')}: ${inter.levels.length} combinations, more than a column can hold`, EXIT.roles)
    }
    added.push(inter)
  }

  // ---- what is not carried is said ----------------------------------------
  const carried = new Set([cluster, sample, condition, ...extras].filter((x): x is string => !!x))
  for (const f of factors) {
    if (carried.has(f.name)) continue
    const inside = interactions.filter(list => list.includes(f.name)).map(list => list.join(':'))
    const used = new Set(Array.from(f.codes).filter(c => c >= 0)).size
    notes.push(inside.length
      ? `${f.name} is a design factor (one value in every sample, ${used} levels) carried only inside ${inside.join(', ')}, not as a column of its own; convert with --extra ${f.name} to group by it alone`
      : `${f.name} is a design factor (one value in every sample, ${used} levels) and is not carried; convert with --extra ${f.name} to carry it`)
  }
  const others = scan.columns.filter(c => c.kind === 'categorical' && !carried.has(c.name) && !isFactor.has(c.name) && !unusable(c, 'extra', n))
  if (others.length) info.push(`also in the object, not carried: ${others.map(c => c.name).join(', ')} (each can be added with --extra)`)

  const crossing = crossCheck(sampleCol, col(condition))
  if (crossing) notes.push(crossing)

  // ---- the embedding --------------------------------------------------------
  let embedding: string
  if (a.embedding !== undefined) {
    const e = scan.embeddings.find(x => x.key === a.embedding && x.nDims >= 2)
    if (!e) {
      throw new CliError(`--embedding ${a.embedding}: this object has no embedding called "${a.embedding}" with two or more dimensions. It has: ${scan.embeddings.map(x => `${x.key} (${x.nDims}-D)`).join(', ') || '(none)'}`, EXIT.roles)
    }
    embedding = e.key
  } else {
    const g = guessEmbedding(scan)
    if (!g) throw new CliError('this object has no embedding of two or more dimensions, and the studio draws every cell on one', EXIT.build)
    embedding = g
    say('embedding', [g], 'a UMAP first, then t-SNE, a force layout, PCA')
  }

  // ---- sentinels ------------------------------------------------------------
  const clusterCol = col(cluster)!
  const present = new Set(clusterCol.levels)
  let sentinels: string[]
  if (a.sentinels !== undefined) {
    const wanted = a.sentinels.length === 1 && a.sentinels[0] === 'none' && !present.has('none') ? [] : [...new Set(a.sentinels)]
    const absent = wanted.filter(s => !present.has(s))
    if (absent.length) info.push(`--sentinel ${absent.join(', ')}: no cell of ${cluster} is labelled so; nothing to mark`)
    sentinels = wanted.filter(s => present.has(s))
  } else {
    sentinels = PROFILE_SENTINELS.filter(s => present.has(s))
    say('sentinel', sentinels, sentinels.length ? "the annotator's labels for cells it excluded or could not resolve" : `no level of ${cluster} is one of ${PROFILE_SENTINELS.join(', ')}`)
  }
  if (sentinels.length) {
    const counts = sentinels.map(s => {
      const code = clusterCol.levels.indexOf(s)
      let k = 0
      for (const c of clusterCol.codes) if (c === code) k++
      return [s, k] as const
    })
    const total = counts.reduce((t, [, k]) => t + k, 0)
    notes.push(`${sentinels.join(' and ')} ${sentinels.length > 1 ? 'are' : 'is'} the annotator's ${sentinels.length > 1 ? 'sentinels' : 'sentinel'}, not cell types: ${total} cells (${counts.map(([s, k]) => `${s} ${k}`).join(', ')}). They stay drawn, and are listed in sentinels so that a reader can leave them out of composition, markers and tests`)
  }

  // ---- the label --------------------------------------------------------------
  let label: string
  if (a.label !== undefined) label = a.label
  else {
    label = (a.input.split(/[\\/]/).pop() ?? a.input).replace(/\.(h5ad|h5)$/i, '')
    say('label', [label], 'the file name')
  }

  return {
    choices: {
      cluster, sample, condition, embedding, label,
      extras: [...extras, ...added.map(c => c.name)],
      sentinels,
    },
    added, interactions, sentinels, guessed, said, notes, info,
  }
}
