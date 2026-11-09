// The design: which of an object's columns describe its samples.
//
// A bundle has one condition, but an experiment can have several factors —
// two crossed treatments, a genotype beside a diet — and an object carries each
// as its own obs column. Choosing the condition by name kept one of them and
// dropped the rest without a word, and an object whose factors had neutral
// names (factor_a, factor_b) came out as a single group.
//
// So a factor is recognised by what it is rather than what it is called: a
// grouping that never varies inside a sample. Cell types vary inside a sample;
// a treatment, a genotype, an age or a batch does not. That test needs only the
// sample column and the codes already in hand, and it holds whatever language
// the columns are named in.
//
// Nothing here reads a file or decides for the user. It answers questions; the
// caller decides what to carry and says what it decided.

import { samePartition, unusable, type Column } from './scan.ts'

/** Does `col` take one value in each sample? NA counts as a value. */
export function sampleLevel(col: Column, sample: Column): boolean {
  if (col.kind !== 'categorical' || sample.kind !== 'categorical') return false
  if (col.codes.length !== sample.codes.length) return false
  const seen = new Map<number, number>()
  for (let i = 0; i < sample.codes.length; i++) {
    const s = sample.codes[i]
    const had = seen.get(s)
    if (had === undefined) seen.set(s, col.codes[i])
    else if (had !== col.codes[i]) return false
  }
  return true
}

/**
 * The design factors, in the object's own column order.
 *
 * A column qualifies when it could be a condition at all (a grouping of 2 to
 * 50 levels, not one value per cell) and takes one value in every sample. The
 * sample column under another name — a donor id that maps one to one onto the
 * samples — is an identifier, not a factor, and is left out. `skip` withdraws
 * columns already given another role.
 *
 * With no sample column there is nothing to test against, and nothing here
 * pretends otherwise: the answer is empty.
 */
export function designFactors(
  columns: Column[], sample: Column | undefined, nCells: number,
  skip: (string | null | undefined)[] = [],
): Column[] {
  if (!sample || sample.kind !== 'categorical') return []
  return columns.filter(c =>
    c.name !== sample.name && !skip.includes(c.name) &&
    !unusable(c, 'condition', nCells) &&
    sampleLevel(c, sample) && !samePartition(c, sample))
}

/** Levels that at least one cell holds. */
function levelsInUse(col: Column): number {
  const used = new Set<number>()
  for (const c of col.codes) if (c >= 0) used.add(c)
  return used.size
}

/**
 * The interaction of two or more groupings, as a grouping of its own.
 *
 * Named in formula notation (`factor_a:factor_b`) with levels `a1:b1`, ordered
 * with the first factor slowest — the order a design table is written in. Only
 * arms that hold cells become levels, so a fractional design gains no empty
 * groups. A cell missing either factor is missing here too.
 */
export function interactionColumn(factors: Column[]): Column {
  const n = factors[0]?.codes.length ?? 0
  const sizes = factors.map(f => f.levels.length)
  // Mixed-radix position of each cell's arm. Float64, not Int32: five factors
  // of fifty levels already pass 2^31.
  const arm = new Float64Array(n)
  let missing = 0
  for (let i = 0; i < n; i++) {
    let at = 0
    for (let k = 0; k < factors.length; k++) {
      const c = factors[k].codes[i]
      if (c < 0) { at = -1; break }
      at = at * sizes[k] + c
    }
    arm[i] = at
    if (at < 0) missing++
  }
  const used = [...new Set(arm)].filter(a => a >= 0).sort((a, b) => a - b)
  const index = new Map(used.map((a, i) => [a, i]))
  const levels = used.map(a => {
    const parts: string[] = new Array(factors.length)
    for (let k = factors.length - 1; k >= 0; k--) {
      parts[k] = factors[k].levels[a % sizes[k]]
      a = Math.floor(a / sizes[k])
    }
    return parts.join(':')
  })
  const codes = new Int32Array(n)
  for (let i = 0; i < n; i++) codes[i] = arm[i] < 0 ? -1 : index.get(arm[i])!
  return {
    name: factors.map(f => f.name).join(':'), kind: 'categorical', levels, codes,
    values: new Float32Array(0), missing,
  }
}

/**
 * Do these factors cross, rather than nest?
 *
 * Crossed factors make more arms together than any of them alone — 2 × 2 is
 * four groups, against two. If one is nested inside another, their interaction
 * is just the finer factor again, and carrying it would add nothing.
 */
export function crossed(factors: Column[]): boolean {
  if (factors.length < 2) return false
  return levelsInUse(interactionColumn(factors)) > Math.max(...factors.map(levelsInUse))
}
