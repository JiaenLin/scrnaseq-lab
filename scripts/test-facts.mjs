// The facts walk: every statement a bundle makes about its object has its source.
//
// It walks FACTS, the closed table of src/lib/facts.ts, so a new row needs a
// RIGHT object and a DECOY here or this test fails. Each row is converted by the
// command from objects written by anndata (scripts/make-facts-fixtures.py):
// the RIGHT one records the fact at the place the lab reads; the DECOY holds
// only a NAME that looks like it (an embedding named like a method, an obs
// column named like a doublet tool, an uns key named like an ambient tool or a
// normalisation, a layer named "counts" holding log values, a two-column key
// holding "pca"). No row may state a value from a name, in the facts or in the
// older provenance fields. It reads only the bundle the command writes, so it
// runs on an older lab and fails there, row by row, for its reason.
//
// Then: the source holds no rule that decides a fact from a name; a fact cannot
// be made without its source; two places recording one fact, and hostile
// values at the place; the page's sentences against the bundle the plan built;
// a split object said, and its facts against the whole conversion's.

import { spawnSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import { mkdtempSync, readFileSync, readdirSync, rmSync, statSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { unzipSync, strFromU8 } from 'fflate'

let failed = 0
const report = (name, ok, detail) => {
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        ${detail}`}`)
}
const check = (name, got, want) => report(name, JSON.stringify(got) === JSON.stringify(want),
  `got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`)
const attempt = (name, fn, want) => {
  try { check(name, fn(), want) } catch (e) { report(name, false, `threw ${e?.stack?.split('\n').slice(0, 2).join(' ') ?? e}`) }
}
const load = async (rel) => {
  try { return await import(rel) } catch (e) { console.log(`        (could not import ${rel}: ${String(e.message).split('\n')[0]})`); return null }
}
const here = (rel) => fileURLToPath(new URL(rel, import.meta.url))
const sha256 = (b) => createHash('sha256').update(b).digest('hex')
const FIX = (name) => here(`./fixtures/facts-${name}.h5ad`)
const tmp = mkdtempSync(join(tmpdir(), 'lab-facts-'))

/** The command, as a person runs it; the collection it wrote, read back. */
function convert(input, name, args = []) {
  const out = join(tmp, `${name}.zip`)
  const r = spawnSync(process.execPath, [here('../bin/scrnaseq-lab.mjs'), 'convert', input, out, '--label', name, ...args],
    { encoding: 'utf8', maxBuffer: 64 << 20 })
  if (r.status !== 0) { console.log(`        (${name}: exit ${r.status})\n${r.stderr.split('\n').slice(-6).join('\n')}`); return null }
  return { out, bytes: readFileSync(out), ...read(readFileSync(out)) }
}
function read(bytes) {
  const c = unzipSync(bytes)
  const index = JSON.parse(strFromU8(c['collection.json']))
  const parts = index.parts.map(p => { const f = unzipSync(c[p.file]); return { files: f, meta: JSON.parse(strFromU8(f['meta.json'])) } })
  return { index, parts, meta: parts[0].meta }
}
const fact = (b, id) => b?.meta.facts?.[id] ?? null
const legacy = (b, k) => (b ? b.meta.provenance?.[k] ?? null : 'no bundle')

const B = {
  a: convert(FIX('a'), 'a'),
  aFlag: convert(FIX('a'), 'a-flag', ['--cluster', 'cell_type', '--embedding', 'X_umap', '--sentinel', 'EXCLUDED', '--sentinel', 'UNRESOLVED']),
  b: convert(FIX('b-decoy'), 'b'),
  e: convert(FIX('e-records'), 'e'),
  f: convert(FIX('f-two'), 'f'),
  g: convert(FIX('g-lognorm'), 'g'),
  h: convert(FIX('h-hostile'), 'h'),
}

const notRec = (searched) => ({ kind: 'recorded', value: null, searched })
const DICT = (k) => `a dict directly under uns holding ${k}`

// One RIGHT and one DECOY for every row of the table: [the bundle, what it states, what it must state].
const ROWS = {
  integration: {
    right: [B.a, b => [fact(b, 'integration'), legacy(b, 'integration')],
      [{ kind: 'recorded', value: 'scanvi', at: 'uns/scintegrate/default_method' }, 'scanvi']],
    decoy: [B.b, b => [fact(b, 'integration'), legacy(b, 'integration')], [notRec(DICT('default_method')), null]],
  },
  integration_embedding: {
    right: [B.a, b => fact(b, 'integration_embedding'), { kind: 'recorded', value: 'X_scanvi', at: 'uns/scintegrate/default_embedding' }],
    decoy: [B.b, b => fact(b, 'integration_embedding'), notRec(DICT('default_embedding'))],
  },
  normalization: {
    right: [B.e, b => [fact(b, 'normalization'), legacy(b, 'normalization')], [{ kind: 'recorded', value: 'log1p, base 2', at: 'uns/log1p' }, 'log1p, base 2']],
    decoy: [B.b, b => [fact(b, 'normalization'), legacy(b, 'normalization')], [notRec('uns/log1p'), null]],
  },
  doublets: {
    right: [B.e, b => fact(b, 'doublets'), { kind: 'recorded', value: 'scrublet', at: 'uns/qc_stage/doublet_method' }],
    decoy: [B.b, b => [fact(b, 'doublets'), legacy(b, 'doublets')], [notRec(DICT('doublet_method')), null]],
  },
  ambient: {
    right: [B.e, b => fact(b, 'ambient'), { kind: 'recorded', value: 'cellbender', at: 'uns/qc_stage/ambient_method' }],
    decoy: [B.b, b => [fact(b, 'ambient'), legacy(b, 'ambient')], [notRec(DICT('ambient_method')), null]],
  },
  expression: {
    right: [B.a, b => [fact(b, 'expression')?.value, fact(b, 'expression')?.of, b?.meta.expression],
      ['log1p of counts scaled to 10,000 per cell', 'X', 'log1p(CP10K)']],
    // An uns string says "CP10K"; the values do not sum to one total.
    decoy: [B.g, b => [fact(b, 'expression')?.value, b?.meta.expression], ['log-normalized, no constant per-cell total', 'log-normalized']],
  },
  counts: {
    right: [B.a, b => { const f = fact(b, 'counts'); return [f?.kind, f?.value, f?.of, f?.carried, b?.meta.hasRawCounts] },
      ['measured', true, 'layers/counts', true, true]],
    // A layer NAMED counts holding log values is not counts.
    decoy: [B.b, b => { const f = fact(b, 'counts'); return [f?.kind, f?.value, b?.meta.hasRawCounts] }, ['measured', false, false]],
  },
  total: {
    right: [B.a, b => [fact(b, 'total')?.kind, fact(b, 'total')?.of], ['measured', 'layers/counts']],
    decoy: [B.b, b => [fact(b, 'total')?.kind, fact(b, 'total')?.of], ['measured', 'X']],
  },
  genes_detected: {
    right: [B.e, b => fact(b, 'genes_detected'), { kind: 'recorded', value: 'n_genes_by_counts', at: 'obs/n_genes_by_counts' }],
    decoy: [B.b, b => [fact(b, 'genes_detected')?.kind, fact(b, 'genes_detected')?.of], ['measured', 'X']],
  },
  mito: {
    right: [B.e, b => fact(b, 'mito'), { kind: 'recorded', value: 'pct_counts_mt', at: 'obs/pct_counts_mt' }],
    // obs/percent_mito_estimate: named like it, not at a place the lab reads.
    decoy: [B.b, b => fact(b, 'mito')?.value, null],
  },
  cell_type_column: {
    right: [B.aFlag, b => fact(b, 'cell_type_column'), { kind: 'chosen', value: 'cell_type', by: 'flag' }],
    decoy: [B.a, b => fact(b, 'cell_type_column'), { kind: 'chosen', value: 'cell_type', by: 'guess' }],
  },
  embedding: {
    right: [B.aFlag, b => fact(b, 'embedding'), { kind: 'chosen', value: 'X_umap', by: 'flag' }],
    // A two-column key holding "pca": no sentence about principal components.
    decoy: [B.b, b => [fact(b, 'embedding'), b?.meta.notes.some(s => /principal component/i.test(s))],
      [{ kind: 'chosen', value: 'X_umap', by: 'guess' }, false]],
  },
  not_called: {
    right: [B.aFlag, b => [fact(b, 'not_called'), b?.meta.notes.some(s => /20 have no finite coordinates on X_umap\b/.test(s))],
      [{ kind: 'chosen', value: 'EXCLUDED, UNRESOLVED', by: 'flag' }, true]],
    // Labels spelled as the profile's words: chosen by a guess, and never said "drawn".
    decoy: [B.a, b => [fact(b, 'not_called')?.by, b?.meta.notes.some(s => /stay drawn/.test(s))], ['guess', false]],
  },
  gene_ids: {
    right: [B.a, b => [fact(b, 'gene_ids')?.value, fact(b, 'gene_ids')?.of], ['4 of 40 names are accessions', 'var']],
    // A var column named gene_ids holding symbols changes nothing.
    decoy: [B.b, b => [fact(b, 'gene_ids')?.value, fact(b, 'gene_ids')?.of], ['4 of 40 names are accessions', 'var']],
  },
}

console.log('\nTHE TABLE WALKED: EVERY ROW, A RIGHT OBJECT AND A DECOY')
const factsMod = await load('../src/lib/facts.ts')
const ids = factsMod?.FACT_IDS ?? Object.keys(ROWS)
attempt('every row of FACTS has a right object and a decoy here', () => ids.filter(id => !ROWS[id]?.right || !ROWS[id]?.decoy), [])
for (const id of ids) {
  const row = ROWS[id]
  if (!row) continue
  for (const which of ['right', 'decoy']) {
    const [b, said, want] = row[which]
    attempt(`${id}, ${which === 'right' ? 'RIGHT: at its place' : 'DECOY: only a name'}`, () => (b ? said(b) : 'no bundle'), want)
  }
}
attempt('the decoy object states no value as recorded, in the facts or the older fields', () => [
  Object.entries(B.b?.meta.facts ?? {}).filter(([k, f]) => k !== 'v' && f.kind === 'recorded' && f.value !== null).map(([k]) => k),
  ['normalization', 'integration', 'doublets', 'ambient'].map(k => legacy(B.b, k)),
], [[], [null, null, null, null]])
attempt('every bundle carries facts, version 1', () => Object.values(B).map(b => b?.meta.facts?.v ?? null), Object.values(B).map(() => 1))

console.log('\nNO RULE IN THE SOURCE DECIDES A FACT FROM A NAME')
{
  // A regex literal or a string test holding a method's, a tool's or a
  // normalisation's word: the shape every guess this wave removed had.
  const WORD = /harmony|scanvi|scvi|mnn|bbknn|combat|scanorama|liger|integrat|doublet|scrublet|soupx|cellbender|decontx|ambient|pca|log1p|cp10k|normali/i
  const files = []
  const walk = (d) => { for (const n of readdirSync(d)) { const p = join(d, n); if (statSync(p).isDirectory()) walk(p); else if (/\.tsx?$/.test(n)) files.push(p) } }
  walk(here('../src'))
  const hits = []
  for (const f of files) {
    readFileSync(f, 'utf8').split('\n').forEach((line, i) => {
      const code = line.replace(/\/\/.*$/, '').trim()
      if (!code || code.startsWith('*') || code.startsWith('/*')) return
      const literals = [...code.matchAll(/(?:^|[=(,:!&|?\s])\/((?:\\.|[^/\n])+)\/[dgimsuy]*\s*\.test\(/g)].map(m => m[1])
      const strings = [...code.matchAll(/\.(?:includes|startsWith|endsWith)\(\s*['"`]([^'"`]+)['"`]\s*\)/g)].map(m => m[1])
      for (const s of [...literals, ...strings]) if (WORD.test(s)) hits.push(`${f.split('/src/')[1]}:${i + 1}: ${s}`)
    })
  }
  attempt('no regex or string test over a name holds a method, a tool or a normalisation', () => hits, [])
}

console.log('\nA FACT CANNOT BE MADE WITHOUT ITS SOURCE')
if (factsMod) {
  const refuses = (fn) => { try { fn(); return false } catch (e) { return e instanceof factsMod.FactError } }
  attempt('a recorded fact without its place, a measured one without of or how, a chosen one without by: refused', () => [
    refuses(() => factsMod.recorded('scanvi')), refuses(() => factsMod.notRecorded()), refuses(() => factsMod.measured(true, 'X')),
    refuses(() => factsMod.measured(true, undefined, 'how')), refuses(() => factsMod.chosen('cell_type')),
    refuses(() => factsMod.chosen('cell_type', 'inferred')), refuses(() => factsMod.measured(true, 'X', 'how', { carried: false })),
  ], [true, true, true, true, true, true, true])
  attempt('a value the reader would call unreadable is refused: markup, a line break, 65 characters', () => [
    refuses(() => factsMod.recorded('<b>x</b>', 'uns/a')), refuses(() => factsMod.recorded('a\nb', 'uns/a')),
    refuses(() => factsMod.recorded('s'.repeat(65), 'uns/a')), refuses(() => factsMod.recorded('s'.repeat(64), 'uns/a')),
  ], [true, true, true, false])
  attempt('a kind its row never has stops the walk', () => refuses(() => factsMod.factsOf({ measured: {}, chosen: {},
    recorded: { integration: { kind: 'measured', value: 'x', of: 'X', how: 'h' } } })), true)
} else report('the facts module loads', false, 'no src/lib/facts.ts')

console.log('\nTWO PLACES, AND HOSTILE VALUES AT THE PLACE')
attempt('two dicts recording a method: neither is stated, both are named with their places', () => [
  fact(B.f, 'integration'), legacy(B.f, 'integration'),
  B.f?.meta.notes.some(s => s.startsWith('2 places record an integration method: ') && s.includes('uns/scintegrate/default_method (scanvi)')
    && s.includes('uns/other_stage/default_method (harmony)') && s.endsWith("the bundle states none of them as the object's")),
], [null, null, true])
attempt('a group, an array, a number, 10,000 characters, a line break, markup: none read, each named', () => [
  fact(B.h, 'integration'),
  ['as_group', 'as_array', 'as_number', 'as_long', 'as_lines', 'as_markup'].map(k => B.h?.meta.notes.some(s => s.startsWith(`uns/${k}/default_method `))),
], [notRec(DICT('default_method')), [true, true, true, true, true, true]])
attempt('a recorded embedding the object does not hold is said so', () => B.a?.meta.notes.some(s => /default_embedding names/.test(s)), false)

console.log('\nTHE PAGE SAYS WHAT THE BUNDLE THEN HOLDS')
const planMod = await load('../src/lib/plan.ts')
const h5 = await load('../src/lib/h5ad.ts')
const roles = await load('../src/cli/roles.ts')
const cliMod = await load('../src/cli/convert.ts')
const h5wasm = (await import('h5wasm/node')).default
const lib = await h5wasm.ready
/** The page's sentences: the plan made from the description the worker sends the page. */
function pageOf(input, args = {}, splitAt = undefined) {
  const f = h5.openInput(h5wasm, lib, input)
  try {
    const scan = h5.scanOpenH5ad(f, input.split('/').pop())
    const r = roles.resolveRoles(scan, { input, output: 'x.zip', ...args })
    const m = scan.matrices.find(x => x.nnzPerCell)
    const info = { nCells: scan.nCells, columns: [...scan.columns], ...planMod.planInfo(scan, m ? m.nnzPerCell() : null) }
    return planMod.pageSays(planMod.plan(info, r.choices, splitAt))
  } finally { f.close() }
}
/** What a collection holds, in the page's terms. */
function holds(b) {
  const keys = b.meta.embeddings.map(e => e.key)
  const tsv = b.parts.every(p => !!p.files['pseudobulk.tsv'])
  const perGroup = {}
  for (const p of b.parts) for (const s of p.meta.samples) (perGroup[s.condition] ??= new Set()).add(s.id)
  const reps = Object.values(perGroup).map(s => s.size)
  return { keys, parts: b.parts.length, counts: tsv && b.parts.every(p => p.meta.hasRawCounts), pseudobulk: tsv && reps.length > 1 && reps.every(n => n >= 4) }
}
const said = (p) => ({
  keys: /: (.*); /.exec(p.travel)[1].split(', '), parts: p.parts ? Number(/in (\d+) parts/.exec(p.parts)[1]) : 1,
  counts: p.counts.startsWith('Raw counts travel'), pseudobulk: p.pseudobulk.yes,
})
let split = null
if (planMod && h5 && roles && cliMod) {
  for (const [name, input] of [['a', FIX('a')], ['b', FIX('b-decoy')], ['e', FIX('e-records')], ['g', FIX('g-lognorm')]]) {
    attempt(`${name}: the page's embeddings, parts, raw counts and pseudobulk are the bundle's`,
      () => { const p = said(pageOf(input)); const h = holds(B[name]); return [p, h] }, (() => { const h = holds(B[name]); return [h, h] })())
  }
  // An object large enough to be split, made so by the plan's threshold: the command's own function.
  const out = join(tmp, 'a-split.zip')
  const lines = []
  try {
    await cliMod.convert({ input: FIX('a'), output: out, label: 'a', splitAt: 2000 }, s => lines.push(s))
    split = { out, bytes: readFileSync(out), ...read(readFileSync(out)) }
  } catch (e) { console.log(`        (split: ${e.message})`) }
  attempt('a-split: the page\'s sentences are the split bundle\'s', () => [said(pageOf(FIX('a'), {}, 2000)), split && holds(split)],
    split ? [holds(split), holds(split)] : ['split', 'split'])
  attempt('a-split: the page says the counts will not be in the bundle, and why', () => pageOf(FIX('a'), {}, 2000).counts,
    `Your object holds raw counts at layers/counts; they will not be in the bundle: the object is converted in ${split?.parts.length} parts, and one matrix is read for them, X, the expression.`)
  attempt('a whole: the page says the counts travel', () => pageOf(FIX('a')).counts, 'Raw counts travel: layers/counts, summed per sample and cell type.')
  attempt('a: two embeddings travel, the latent spaces and the copy are named as not offered', () => [pageOf(FIX('a')).travel, pageOf(FIX('a')).notOffered], [
    '2 embeddings travel in the file: X_umap, X_umap_harmony; X_umap opens first.',
    'Not offered beside them: X_harmony (10 columns, a latent space), X_scanvi (10 columns, a latent space) and X_umap_scanvi (the same coordinates as X_umap).'])
} else report('the plan, the reader and the command load', false, 'a module did not import')

console.log('\nA SPLIT OBJECT SAYS ITS COUNTS ARE NOT CARRIED; ITS OTHER FACTS ARE THE WHOLE\'S')
attempt('it is split, and no part carries counts', () => [!!split && split.parts.length > 1, split?.parts.every(p => !p.files['pseudobulk.tsv'] && !p.meta.hasRawCounts)], [true, true])
attempt('every part and the collection say the counts are in the object and not carried, with why', () => {
  const said = (notes) => notes.some(s => /^the object holds raw counts at layers\/counts, and this bundle does not carry them: the object is converted in \d+ parts/.test(s))
  return [split?.parts.map(p => said(p.meta.notes)).every(Boolean), said(split?.index.notes ?? []),
    split?.parts.map(p => [p.meta.facts.counts.value, p.meta.facts.counts.carried, /converted in \d+ parts/.test(p.meta.facts.counts.why)]).every(x => x.join() === 'true,false,true')]
}, [true, true, true])
attempt('every other fact equals the whole conversion\'s, but the totals\' matrix (the counts are not read)', () => {
  const strip = (f) => { const { counts: _c, total: _t, genes_detected: _g, ...rest } = f; return JSON.stringify(rest) }
  return [split?.parts.every(p => strip(p.meta.facts) === strip(B.a.meta.facts)), split?.parts.map(p => p.meta.facts.total.of).every(o => o === 'X'), B.a?.meta.facts.total.of]
}, [true, true, 'layers/counts'])

console.log('\nTHE SAME BYTES TWICE')
{
  const again = convert(FIX('a'), 'a')
  attempt('two conversions of one object are byte-identical', () => !!again && sha256(again.bytes) === sha256(B.a.bytes), true)
}

rmSync(tmp, { recursive: true, force: true })
console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll facts tests passed\n')
process.exit(failed ? 1 : 0)
