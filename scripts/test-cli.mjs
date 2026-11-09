// The command-line converter, run for real.
//
// Every conversion here goes through `bin/scrnaseq-lab.mjs` in a child
// process, on the committed fixtures in scripts/fixtures/ — synthetic objects
// shaped like the single-cell family's final object (see
// scripts/make-family-fixture.py). Nothing is mocked between the command and
// the bytes it writes, because the thing being promised is the file a person
// or a pipeline gets back.
//
// What is checked, and why each one exists:
//
//   - anndata's default string encoding is read. On pandas >= 3 it stores
//     obs/var names as nullable-string-array groups, and the reader used to
//     crash on the first one it met.
//   - design factors are found by what they are — constant within every
//     sample — not by what they are called. Neutral names used to give a
//     single-group bundle, and a second factor was dropped without a word.
//   - every design factor is carried, the interaction too, and one that is
//     left out is named in the notes.
//   - the same input gives the same bytes: in another process, at another
//     time, in another time zone.
//   - a role the caller did not name is guessed out loud, and one the caller
//     did name is never second-guessed.
//   - each kind of failure has its own exit code, so a pipeline can tell a
//     typo from a broken file.

import { spawnSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import {
  copyFileSync, existsSync, lstatSync, mkdirSync, mkdtempSync, readFileSync, rmSync, symlinkSync,
} from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { unzipSync, strFromU8 } from 'fflate'
import h5wasm from 'h5wasm/node'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { buildBundle } from '../src/lib/build.ts'
import { categorical } from '../src/lib/scan.ts'
import { readCollectionIndex, readEntry } from '../src/lib/collection.ts'

let failed = 0
const report = (name, ok, detail) => {
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        ${detail}`}`)
}
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  report(name, ok, `got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`)
}
const matches = (name, text, re) => report(name, re.test(String(text)), `no match for ${re} in:\n${String(text).slice(0, 2000)}`)
/** A check whose setup may throw: a throw is a failure of this check, not of the file. */
const attempt = (name, fn, want) => {
  try {
    check(name, fn(), want)
  } catch (e) {
    report(name, false, `threw ${e?.constructor?.name ?? 'error'}: ${e?.message ?? e}`)
  }
}
/** Import a module this file tests, or null — so a missing one fails its checks, not the run. */
const tryImport = async (rel) => {
  try { return await import(rel) } catch (e) { console.log(`        (could not import ${rel}: ${e.message.split('\n')[0]})`); return null }
}

const here = (rel) => fileURLToPath(new URL(rel, import.meta.url))
const FIX = {
  classic: here('./fixtures/family-2x2.h5ad'),
  dflt: here('./fixtures/family-2x2-default.h5ad'),
  mask: here('./fixtures/nullable-mask.h5ad'),
  scaled: here('./fixtures/scaled-only.h5ad'),
}
const BIN = here('../bin/scrnaseq-lab.mjs')
const tmp = mkdtempSync(join(tmpdir(), 'lab-cli-test-'))
const sha256 = (b) => createHash('sha256').update(b).digest('hex')

/** Run the real command. */
function cli(args, env = {}) {
  const r = spawnSync(process.execPath, [BIN, ...args], {
    encoding: 'utf8', env: { ...process.env, ...env }, maxBuffer: 64 << 20,
  })
  let summary = null
  try { summary = JSON.parse(r.stdout) } catch { /* not a summary */ }
  return { code: r.status, out: r.stdout ?? '', err: r.stderr ?? '', summary }
}

/** A finished file, read back the way the studio reads it: index first, then parts. */
async function openOutput(path) {
  const bytes = readFileSync(path)
  const blob = new Blob([bytes])
  const idx = await readCollectionIndex(blob)
  if (!idx) throw new Error(`${path} is not a collection`)
  const parts = []
  for (const info of idx.meta.parts) {
    const raw = await readEntry(blob, idx.entries.get(info.file))
    const files = unzipSync(raw)
    parts.push({ info, raw, files, meta: JSON.parse(strFromU8(files['meta.json'])) })
  }
  return { bytes, meta: idx.meta, parts }
}
const tryOpen = async (path) => { try { return await openOutput(path) } catch { return null } }

/**
 * The checks the studio's bundle reader makes before it draws anything
 * (scrnaseq-studio src/lib/bundle.ts parseBundle, readEmbeddings, readExtras),
 * restated so a bundle that would be refused there is refused here.
 */
function studioProblems(files, meta) {
  const p = []
  const need = (n) => files[n] ?? (p.push(`no ${n}`), new Uint8Array(0))
  const view = (C, n, w) => {
    const b = need(n)
    if (b.length % w) p.push(`${n} is not a whole number of ${w}-byte values`)
    return new C(b.slice(0, b.length - (b.length % w)).buffer)
  }
  if (meta.schema !== 'scrnaseq-studio/bundle@1') p.push(`schema ${meta.schema}`)
  const genes = strFromU8(need('genes.txt')).split('\n').filter(g => g.length)
  if (genes.length !== meta.nGenes) p.push(`genes.txt has ${genes.length}, meta says ${meta.nGenes}`)
  const n = meta.nCells
  const indptr = view(Int32Array, 'expr.indptr.i32', 4)
  const indices = view(Int32Array, 'expr.indices.i32', 4)
  const data = view(Float32Array, 'expr.data.f32', 4)
  const cluster = view(Uint16Array, 'cluster.u16', 2)
  const sample = view(Uint16Array, 'sample.u16', 2)
  const embed = view(Float32Array, 'embed.f32', 4)
  const qc = view(Float32Array, 'qc.f32', 4)
  if (indptr.length !== meta.nGenes + 1) p.push('expr.indptr length')
  if (indptr[indptr.length - 1] !== indices.length) p.push('expr.indptr does not end at nnz')
  if (indices.length !== data.length) p.push('indices and data lengths differ')
  if (cluster.length !== n || sample.length !== n) p.push('cluster/sample length')
  if (embed.length !== 2 * n) p.push('embed.f32 length')
  if (qc.length !== 3 * n) p.push('qc.f32 length')
  if (!meta.clusters?.length) p.push('no clusters')
  if (!meta.samples?.length) p.push('no samples')
  for (const i of indices) if (i < 0 || i >= n) { p.push('a cell index outside the matrix'); break }
  for (let i = 0; i < n; i++) {
    if (cluster[i] > meta.clusters.length - 1) { p.push('cluster code past the clusters'); break }
    if (sample[i] >= meta.samples.length) { p.push('sample code past the samples'); break }
  }
  for (const e of (meta.embeddings ?? []).slice(1)) {
    if (!files[e.file]) p.push(`embedding ${e.key} has no ${e.file}`)
    else if (files[e.file].length !== 8 * n) p.push(`${e.file} length`)
  }
  for (const e of Array.isArray(meta.extras) ? meta.extras : []) {
    const codes = view(Uint16Array, e.file, 2)
    if (codes.length !== n) p.push(`${e.file} length`)
    if (!e.levels?.length) p.push(`${e.key} has no levels`)
    for (const c of codes) if (c >= e.levels.length) { p.push(`${e.file} refers past its levels`); break }
  }
  return p
}

// The fixture's design, as the generator wrote it: samples s1..s8, two per arm,
// arms in the order a1b1, a1b2, a2b1, a2b2.
const DESIGN = {
  s1: ['a1', 'b1'], s2: ['a1', 'b1'], s3: ['a1', 'b2'], s4: ['a1', 'b2'],
  s5: ['a2', 'b1'], s6: ['a2', 'b1'], s7: ['a2', 'b2'], s8: ['a2', 'b2'],
}
const EXPLICIT = [
  '--cluster', 'cell_type', '--sample', 'sample', '--condition', 'factor_a',
  '--extra', 'factor_b', '--interaction', 'factor_a,factor_b', '--embedding', 'X_umap',
  '--sentinel', 'EXCLUDED', '--sentinel', 'UNRESOLVED', '--label', 'family',
]

await h5wasm.ready
/** Scan, run `fn`, close. The scan's matrices read lazily, so `fn` runs while the file is open. */
function withScan(path, fn) {
  const f = new h5wasm.File(path, 'r')
  try { return fn(scanOpenH5ad(f, path.split(/[\\/]/).pop())) } finally { f.close() }
}
const scanFile = (path) => withScan(path, s => s)

// ---------------------------------------------------------------------------
console.log("\nTHE READER TAKES ANNDATA'S DEFAULT STRING ENCODING")
{
  const classic = scanFile(FIX.classic)
  let dflt = null
  attempt('an object written with the defaults (nullable strings) scans', () => {
    dflt = scanFile(FIX.dflt)
    return dflt.nCells
  }, 240)
  attempt('its gene names are the classic copy\'s, in order',
    () => dflt.geneSets[Object.keys(dflt.geneSets)[0]], classic.geneSets[Object.keys(classic.geneSets)[0]])
  attempt('its var gene_id column is found as the accessions',
    () => { const a = Object.values(dflt.geneAliases)[0]; return [a.kind, a.column, a.names.slice(0, 2)] },
    ['accession', 'gene_id', ['ENSXXXG00000000001', 'ENSXXXG00000000002']])
  attempt('its obs columns and levels are the classic copy\'s',
    () => dflt.columns.map(c => [c.name, c.kind, c.levels]),
    classic.columns.map(c => [c.name, c.kind, c.levels]))

  let mask = null
  attempt('a masked nullable-string object scans', () => { mask = scanFile(FIX.mask); return mask.nCells }, 8)
  attempt('a masked obs string is a missing value, never a level of its own',
    () => { const c = mask.columns.find(x => x.name === 'note'); return [c.kind, c.levels, c.missing, c.codes[2], c.codes[5]] },
    ['categorical', ['n0', 'n1', 'n3', 'n4', 'n6', 'n7'], 2, -1, -1])
  attempt('a masked var string is empty, so the alias counts it as missing',
    () => Object.values(mask.geneAliases)[0].names,
    ['ENSXXXG00000000001', '', 'ENSXXXG00000000003', 'ENSXXXG00000000004'])
}

// ---------------------------------------------------------------------------
console.log('\nDESIGN FACTORS ARE TOLD BY WHAT THEY ARE, NOT WHAT THEY ARE CALLED')
{
  const design = await tryImport('../src/lib/design.ts')
  // Eight cells, four samples: f1 and f2 cross, `cell` varies inside a sample,
  // `donor` is the sample column under another name, `nested` sits inside f1.
  const sample = categorical('sample', ['p', 'p', 'q', 'q', 'r', 'r', 't', 't'])
  const f1 = categorical('f1', ['x1', 'x1', 'x1', 'x1', 'x2', 'x2', 'x2', 'x2'])
  const f2 = categorical('f2', ['y1', 'y1', 'y2', 'y2', 'y1', 'y1', 'y2', 'y2'])
  const cell = categorical('cell', ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B'])
  const donor = categorical('donor', ['d1', 'd1', 'd2', 'd2', 'd3', 'd3', 'd4', 'd4'])
  const nested = categorical('nested', ['n1', 'n1', 'n2', 'n2', 'n3', 'n3', 'n3', 'n3'])
  attempt('a column constant within every sample is sample-level', () => design.sampleLevel(f2, sample), true)
  attempt('a column that varies inside a sample is not', () => design.sampleLevel(cell, sample), false)
  attempt('the design factors are the sample-level groupings, in obs order — not the sample under another name, not a cell label',
    () => design.designFactors([sample, cell, f1, donor, f2], sample, 8).map(c => c.name), ['f1', 'f2'])
  attempt('a column named like a condition but varying inside samples is not a factor',
    () => design.designFactors([sample, categorical('condition', ['c1', 'c2', 'c1', 'c2', 'c1', 'c2', 'c1', 'c2']), f1], sample, 8)
      .map(c => c.name),
    ['f1'])
  attempt('no sample column, no way to tell a factor from a label', () => design.designFactors([f1, f2], undefined, 8).length, 0)
  attempt('the interaction is named in formula notation and keeps every arm, first factor slowest',
    () => { const i = design.interactionColumn([f1, f2]); return [i.name, i.levels, Array.from(i.codes)] },
    ['f1:f2', ['x1:y1', 'x1:y2', 'x2:y1', 'x2:y2'], [0, 0, 1, 1, 2, 2, 3, 3]])
  attempt('a missing value in either factor is missing in the interaction',
    () => {
      const a = { ...f1, codes: Int32Array.from([0, -1, 0, 0, 1, 1, 1, 1]), missing: 1 }
      const i = design.interactionColumn([a, f2])
      return [i.codes[1], i.missing]
    }, [-1, 1])
  attempt('an arm with no cells is not a level', () => design.interactionColumn([f1, nested]).levels,
    ['x1:n1', 'x1:n2', 'x2:n3'])
  attempt('crossed factors are crossed', () => design.crossed([f1, f2]), true)
  attempt('a nested pair is not — its interaction is just the finer factor', () => design.crossed([f1, nested]), false)
}

// ---------------------------------------------------------------------------
console.log('\nTHE ARGUMENTS MEAN WHAT THEY SAY')
{
  const args = await tryImport('../src/cli/args.ts')
  attempt('positional input and output, and every role flag',
    () => {
      const o = args.parseConvertArgs(['in.h5ad', 'out.zip', '--cluster', 'ct', '--sample=s',
        '--condition', 'g', '--extra', 'x', '--extra=y', '--interaction', 'g, x', '--embedding', 'X_umap',
        '--label', 'L', '--sentinel', 'EXCLUDED'])
      return [o.input, o.output, o.cluster, o.sample, o.condition, o.extras, o.interactions, o.embedding, o.label, o.sentinels]
    },
    ['in.h5ad', 'out.zip', 'ct', 's', 'g', ['x', 'y'], [['g', 'x']], 'X_umap', 'L', ['EXCLUDED']])
  attempt('an absent flag stays absent, so it can be guessed',
    () => { const o = args.parseConvertArgs(['in.h5ad', 'out.zip']); return [o.cluster, o.extras, o.interactions, o.sentinels] },
    [undefined, undefined, undefined, undefined])
  const usage = (argv) => { try { args.parseConvertArgs(argv); return 'accepted' } catch (e) { return e.exit } }
  attempt('a single-valued flag given twice is a usage error', () => usage(['i', 'o', '--cluster', 'a', '--cluster', 'b']), 2)
  attempt('a flag with no value is a usage error', () => usage(['i', 'o', '--cluster']), 2)
  attempt('an unknown flag is a usage error', () => usage(['i', 'o', '--clutser', 'a']), 2)
  attempt('an interaction needs two factors', () => usage(['i', 'o', '--interaction', 'a']), 2)
  attempt('a third positional is a usage error', () => usage(['i', 'o', 'p']), 2)
}

// ---------------------------------------------------------------------------
console.log('\nTHE COMMAND CONVERTS A FAMILY-SHAPED OBJECT, DESIGN AND ALL')
const OUT1 = join(tmp, 'explicit.zip')
const r1 = cli(['convert', FIX.classic, OUT1, ...EXPLICIT])
check('it exits 0', r1.code, 0)
if (r1.code !== 0) console.log(r1.err.split('\n').slice(-12).join('\n'))
const o1 = await tryOpen(OUT1)
report('the output is a collection the studio\'s index reader opens', !!o1, 'no collection')
attempt('the summary on stdout carries the file\'s own sha256', () => r1.summary.sha256, o1 && sha256(o1.bytes))
const p1 = o1?.parts[0]
attempt('one part, a complete bundle of every cell and gene', () => [o1.parts.length, p1.meta.schema, p1.meta.nCells, p1.meta.nGenes],
  [1, 'scrnaseq-studio/bundle@1', 240, 48])
attempt('the studio\'s own bundle checks find nothing wrong', () => studioProblems(p1.files, p1.meta), [])
attempt('the condition is the first factor, and every sample keeps its own level',
  () => [p1.meta.provenance.condition, p1.meta.conditions, p1.meta.samples.map(s => [s.id, s.condition])],
  ['factor_a', ['a1', 'a2'], Object.entries(DESIGN).map(([s, [a]]) => [s, a])])
attempt('the second factor and the interaction travel as columns the studio can group by',
  () => p1.meta.extras.map(e => [e.key, e.levels]),
  [['factor_b', ['b1', 'b2']], ['factor_a:factor_b', ['a1:b1', 'a1:b2', 'a2:b1', 'a2:b2']]])
attempt('every cell carries its own sample\'s factor_b and arm',
  () => {
    const sample = new Uint16Array(p1.files['sample.u16'].slice().buffer)
    const fb = new Uint16Array(p1.files[p1.meta.extras[0].file].slice().buffer)
    const arm = new Uint16Array(p1.files[p1.meta.extras[1].file].slice().buffer)
    let wrong = 0
    for (let i = 0; i < sample.length; i++) {
      const [a, b] = DESIGN[p1.meta.samples[sample[i]].id]
      if (p1.meta.extras[0].levels[fb[i]] !== b || p1.meta.extras[1].levels[arm[i]] !== `${a}:${b}`) wrong++
    }
    return wrong
  }, 0)
attempt('nothing is said to be left out', () => p1.meta.notes.filter(s => /not carried/.test(s)), [])
attempt('the annotator\'s sentinels are marked, in the index and in the part',
  () => [o1.meta.sentinels, p1.meta.sentinels], [['EXCLUDED', 'UNRESOLVED'], ['EXCLUDED', 'UNRESOLVED']])
attempt('and counted in the notes', () => p1.meta.notes.some(s => /EXCLUDED/.test(s) && /UNRESOLVED/.test(s) && /\d+ cells/.test(s)), true)
attempt('the index keeps the whole object\'s level order for the condition and every extra',
  () => [o1.meta.condOrder, o1.meta.extraOrder], [['a1', 'a2'], { factor_b: ['b1', 'b2'], 'factor_a:factor_b': ['a1:b1', 'a1:b2', 'a2:b1', 'a2:b2'] }])
attempt('the run\'s uns record (nested dicts, string lists, an empty list) was read past without harm',
  () => [r1.code, p1.meta.provenance.integration], [0, 'harmony'])

// ---------------------------------------------------------------------------
console.log("\nANNDATA'S DEFAULT ENCODING CONVERTS TO THE SAME BUNDLE")
const OUT2 = join(tmp, 'default.zip')
const r2 = cli(['convert', FIX.dflt, OUT2, ...EXPLICIT])
check('it exits 0', r2.code, 0)
if (r2.code !== 0) console.log(r2.err.split('\n').slice(-12).join('\n'))
const o2 = await tryOpen(OUT2)
attempt('every entry but meta.json is byte-identical to the classic copy\'s',
  () => {
    const a = o1.parts[0].files, b = o2.parts[0].files
    const names = [...new Set([...Object.keys(a), ...Object.keys(b)])].filter(n => n !== 'meta.json').sort()
    return names.filter(n => !a[n] || !b[n] || sha256(a[n]) !== sha256(b[n]))
  }, [])
attempt('meta.json differs only in naming its source file',
  () => {
    const a = { ...o1.parts[0].meta }, b = { ...o2.parts[0].meta }
    delete a.source; delete b.source
    return JSON.stringify(a) === JSON.stringify(b)
  }, true)

// ---------------------------------------------------------------------------
console.log('\nTHE SAME INPUT GIVES THE SAME BYTES')
const OUT3 = join(tmp, 'again-utc.zip')
const OUT4 = join(tmp, 'again-far.zip')
const r3 = cli(['convert', FIX.classic, OUT3, ...EXPLICIT], { TZ: 'UTC' })
const r4 = cli(['convert', FIX.classic, OUT4, ...EXPLICIT], { TZ: 'Pacific/Kiritimati' })
check('both runs exit 0', [r3.code, r4.code], [0, 0])
attempt('a second conversion, in another process, gives the same sha256',
  () => sha256(readFileSync(OUT3)), o1 && sha256(o1.bytes))
attempt('so does a third in a time zone fourteen hours away',
  () => sha256(readFileSync(OUT4)), o1 && sha256(o1.bytes))
attempt('collection.json lists each part\'s sha256, and it is the part\'s',
  () => o1.parts.map(p => p.info.sha256 === sha256(p.raw) && /^[0-9a-f]{64}$/.test(p.info.sha256)), [true])
{
  // The builder alone, with the clock moved a day between two builds of the
  // same scan: what used to differ was the zip entries' modification time.
  attempt('a bundle built a day later is byte-identical',
    () => withScan(FIX.classic, scan => {
      const choices = { cluster: 'cell_type', sample: 'sample', condition: 'factor_a', embedding: 'X_umap', label: 'L', extras: [] }
      const realNow = Date.now
      try {
        const a = buildBundle(scan, choices).zip
        Date.now = () => realNow() + 86_400_000
        const b = buildBundle(scan, choices).zip
        return sha256(a) === sha256(b)
      } finally {
        Date.now = realNow
      }
    }), true)
}

// ---------------------------------------------------------------------------
console.log('\nA ROLE NOT NAMED IS GUESSED OUT LOUD; A ROLE NAMED IS NEVER SECOND-GUESSED')
const OUT5 = join(tmp, 'guessed.zip')
const r5 = cli(['convert', FIX.classic, OUT5, '--label', 'family'])
check('with no role flags it exits 0', r5.code, 0)
matches('it says which cell annotation it guessed', r5.err, /^guessed --cluster cell_type\b/m)
matches('which sample column', r5.err, /^guessed --sample sample\b/m)
matches('which condition — a neutral name, found by content', r5.err, /^guessed --condition factor_a\b/m)
matches('which other factors it carries', r5.err, /^guessed --extra factor_b\b/m)
matches('which interaction', r5.err, /^guessed --interaction factor_a,factor_b\b/m)
matches('which embedding opens first', r5.err, /^guessed --embedding X_umap\b/m)
matches('which labels are sentinels', r5.err, /^guessed --sentinel EXCLUDED --sentinel UNRESOLVED\b/m)
attempt('the summary lists what was guessed', () => r5.summary.guessed,
  ['cluster', 'sample', 'condition', 'extra', 'interaction', 'embedding', 'sentinel'])
attempt('the guesses are exactly the explicit flags: the same bytes', () => r5.summary.sha256, r1.summary?.sha256)

const OUT6 = join(tmp, 'named.zip')
const r6 = cli(['convert', FIX.classic, OUT6, '--condition', 'factor_b', '--label', 'family'])
check('naming the other factor as the condition exits 0', r6.code, 0)
const o6 = await tryOpen(OUT6)
attempt('the named condition is the one used', () => [o6.parts[0].meta.provenance.condition, o6.parts[0].meta.conditions],
  ['factor_b', ['b1', 'b2']])
report('and it is not reported as a guess', r6.code === 0 && !/^guessed --condition/m.test(r6.err), r6.err.slice(0, 600))
attempt('the factor it displaced is carried instead, and the interaction keeps obs order',
  () => o6.parts[0].meta.extras.map(e => e.key), ['factor_a', 'factor_a:factor_b'])

// ---------------------------------------------------------------------------
console.log('\nA DESIGN FACTOR LEFT OUT IS NAMED, NOT DROPPED')
const OUT7 = join(tmp, 'leftout.zip')
const r7 = cli(['convert', FIX.classic, OUT7, '--condition', 'factor_a', '--extra', 'none', '--interaction', 'none'])
check('carrying one factor only exits 0', r7.code, 0)
const o7 = await tryOpen(OUT7)
attempt('no extra columns were written', () => o7.parts[0].meta.extras, [])
attempt('a note names factor_b as a design factor that is not carried',
  () => o7.parts[0].meta.notes.filter(s => /factor_b/.test(s) && /not carried/.test(s)).length, 1)
matches('and the command says so too', r7.err, /factor_b.*not carried/)
const OUT8 = join(tmp, 'inside.zip')
const r8 = cli(['convert', FIX.classic, OUT8, '--condition', 'factor_a', '--extra', 'none', '--interaction', 'factor_a,factor_b'])
const o8 = await tryOpen(OUT8)
attempt('a factor carried only inside an interaction is said to be',
  () => [r8.code, o8.parts[0].meta.notes.filter(s => /factor_b/.test(s) && /only inside factor_a:factor_b/.test(s)).length], [0, 1])

// ---------------------------------------------------------------------------
console.log('\nEACH KIND OF FAILURE HAS ITS OWN EXIT CODE')
const code = (...args) => cli(args).code
check('no arguments: 2 (usage)', code('convert'), 2)
check('an unknown flag: 2 (usage)', code('convert', FIX.classic, join(tmp, 'x.zip'), '--bogus', 'x'), 2)
check('an unknown command: 2 (usage)', code('transmogrify'), 2)
check('an input that does not exist: 3 (input)', code('convert', join(tmp, 'nope.h5ad'), join(tmp, 'x.zip')), 3)
check('a file that is not HDF5: 3 (input)', code('convert', here('./test-cli.mjs'), join(tmp, 'x.zip')), 3)
const rMissing = cli(['convert', FIX.classic, join(tmp, 'x.zip'), '--condition', 'no_such_column'])
check('a role naming a column the object lacks: 4 (roles)', rMissing.code, 4)
matches('and the message names it and lists what there is', rMissing.err, /no_such_column[\s\S]*factor_a/)
check('an interaction of a column the object lacks: 4 (roles)',
  code('convert', FIX.classic, join(tmp, 'x.zip'), '--interaction', 'factor_a,nope'), 4)
check('an embedding the object lacks: 4 (roles)', code('convert', FIX.classic, join(tmp, 'x.zip'), '--embedding', 'X_nope'), 4)
const rScaled = cli(['convert', FIX.scaled, join(tmp, 'x.zip')])
check('an object whose only matrix is scaled: 5 (cannot convert)', rScaled.code, 5)
matches('and the message says why', rScaled.err, /scaled/)
check('an output directory that does not exist: 6 (output)',
  code('convert', FIX.classic, join(tmp, 'no', 'such', 'dir', 'x.zip'), ...EXPLICIT), 6)
report('a failed conversion leaves no file behind',
  r1.code === 0 && !existsSync(join(tmp, 'x.zip')) && !existsSync(join(tmp, 'x.zip.partial')),
  'x.zip or x.zip.partial exists, or the command never ran')
check('--help: 0', code('--help'), 0)

// ---------------------------------------------------------------------------
console.log('\nTHE INPUT IS NEVER OVERWRITTEN, HOWEVER THE OUTPUT REACHES IT')
{
  // The output is renamed into place, so an output that IS the input — by the
  // same path, or through a link — would replace the object with its own zip.
  // Each case converts its own copy, so one case's damage cannot decide another's.
  const before = sha256(readFileSync(FIX.classic))
  const self = join(tmp, 'self.h5ad')
  copyFileSync(FIX.classic, self)
  const rs = cli(['convert', self, self, ...EXPLICIT])
  check('converting a copy onto itself: 2 (usage), and the copy is untouched',
    [rs.code, sha256(readFileSync(self))], [2, before])
  check('the refusal names the file',
    [rs.code, /^scrnaseq-lab: .*self\.h5ad/m.test(rs.err)], [2, true])

  const linked = join(tmp, 'linked.h5ad')
  copyFileSync(FIX.classic, linked)
  const link = join(tmp, 'link.zip')
  symlinkSync(linked, link)
  const rl = cli(['convert', linked, link, ...EXPLICIT])
  check('onto a symlink to it: 2, the copy untouched and the link still a link',
    [rl.code, sha256(readFileSync(linked)), lstatSync(link).isSymbolicLink()], [2, before, true])

  const real = join(tmp, 'real')
  mkdirSync(real)
  const inReal = join(real, 'obj.h5ad')
  copyFileSync(FIX.classic, inReal)
  symlinkSync(real, join(tmp, 'alias'))
  const rd = cli(['convert', inReal, join(tmp, 'alias', 'obj.h5ad'), ...EXPLICIT])
  check('onto itself through a symlinked directory: 2, and untouched',
    [rd.code, sha256(readFileSync(inReal))], [2, before])

  const partial = join(tmp, 'q.zip.partial')
  copyFileSync(FIX.classic, partial)
  const rp = cli(['convert', partial, join(tmp, 'q.zip'), ...EXPLICIT])
  check("an input that is the output's .partial: 2, and untouched",
    [rp.code, existsSync(partial) && sha256(readFileSync(partial))], [2, before])
}

// ---------------------------------------------------------------------------
console.log('\nA --sentinel LEVEL THAT NO CELL CARRIES IS REFUSED, NOT SKIPPED')
{
  // Skipped, a typo left the real sentinel counted as a cell type and exited 0.
  const out = join(tmp, 'typo.zip')
  const rt = cli(['convert', FIX.classic, out, ...EXPLICIT, '--sentinel', 'NOPE', '--sentinel', 'TYPO'])
  check('a --sentinel no cell carries: 4 (roles), and every such level is named',
    [rt.code, /NOPE/.test(rt.err) && /TYPO/.test(rt.err)], [4, true])
  check('the message lists the levels that cells do carry',
    [rt.code, /TypeA[\s\S]*EXCLUDED[\s\S]*UNRESOLVED/.test(rt.err)], [4, true])
  report('and nothing is written', rt.code === 4 && !existsSync(out), `exit ${rt.code}; ${out} ${existsSync(out) ? 'exists' : 'absent'}`)
}

rmSync(tmp, { recursive: true, force: true })
console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll command-line tests passed\n')
process.exit(failed ? 1 : 0)
