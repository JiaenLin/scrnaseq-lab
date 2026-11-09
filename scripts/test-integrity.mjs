// What a conversion must not change about the object: which condition each
// cell is in, which run made it, and which coordinates are a picture of it.
//
// Three defects, each silent, each found by the platform review:
//
//   - A sample whose cells sit in more than one condition. The bundle records
//     one condition per sample, so every cell of such a sample was given its
//     first cell's condition. With no sample column the whole object is one
//     "all cells" sample, and every cell was given the first cell's condition.
//     The per-group counts still looked plausible; nothing said so.
//   - The run's own record. The object says in uns what made it (the tool, its
//     version and commit, the job), and the bundle dropped it.
//   - Latent spaces drawn as embeddings, and one embedding offered twice. A
//     10- or 30-column space was cut to its first two columns and offered as a
//     view of the cells, and X_umap, a byte copy of X_umap_harmony, was offered
//     beside it as a second view.
//
// Everything here is synthetic: the committed fixtures in scripts/fixtures,
// copies of them edited in this file through h5wasm, and Scans built in
// memory. Conversions go through bin/scrnaseq-lab.mjs in a child process.

import { spawnSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import { copyFileSync, mkdtempSync, readFileSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { unzipSync, strFromU8 } from 'fflate'
import h5wasm from 'h5wasm/node'

let failed = 0
const report = (name, ok, detail) => {
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        ${detail}`}`)
}
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  report(name, ok, `got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`)
}
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
/** JSON with object keys sorted, so two records compare by content alone. */
const canon = (v) => JSON.stringify(v, (_k, x) => (x && typeof x === 'object' && !Array.isArray(x)
  ? Object.fromEntries(Object.keys(x).sort().map(k => [k, x[k]])) : x))

const here = (rel) => fileURLToPath(new URL(rel, import.meta.url))
const FIX = {
  classic: here('./fixtures/family-2x2.h5ad'),
  dflt: here('./fixtures/family-2x2-default.h5ad'),
  mask: here('./fixtures/nullable-mask.h5ad'),
  records: here('./fixtures/uns-records.h5ad'),
  unreadable: here('./fixtures/uns-unreadable.h5ad'),
}
const BIN = here('../bin/scrnaseq-lab.mjs')
const tmp = mkdtempSync(join(tmpdir(), 'lab-integrity-test-'))
const sha256 = (b) => createHash('sha256').update(b).digest('hex')

const build = await tryImport('../src/lib/build.ts')
const scanMod = await tryImport('../src/lib/scan.ts')
const shard = await tryImport('../src/lib/shard.ts')
const { categorical, fromFactor } = scanMod

/** Run the real command; a `timeout` (ms) kills it, and then `code` is null. */
function cli(args, env = {}, timeout = undefined) {
  const t0 = Date.now()
  const r = spawnSync(process.execPath, [BIN, ...args], {
    encoding: 'utf8', env: { ...process.env, ...env }, maxBuffer: 64 << 20, timeout,
  })
  let summary = null
  try { summary = JSON.parse(r.stdout) } catch { /* not a summary */ }
  return { code: r.status, out: r.stdout ?? '', err: r.stderr ?? '', summary, seconds: (Date.now() - t0) / 1000 }
}

/** A finished file, unzipped: the index, and each part's entries and meta. */
function openOutput(path) {
  const bytes = readFileSync(path)
  const outer = unzipSync(bytes)
  const index = JSON.parse(strFromU8(outer['collection.json']))
  const parts = index.parts.map(info => {
    const files = unzipSync(outer[info.file])
    return { info, files, meta: JSON.parse(strFromU8(files['meta.json'])) }
  })
  return { bytes, index, parts }
}
const tryOpen = (path) => { try { return openOutput(path) } catch { return null } }

/** Convert, and open what was written. */
function convert(input, name, args, timeout = undefined) {
  const out = join(tmp, name)
  const r = cli(['convert', input, out, ...args], {}, timeout)
  if (r.code !== 0) console.log(`        (${name}: exit ${r.code})\n${r.err.split('\n').slice(-8).join('\n')}`)
  return { r, o: r.code === 0 ? tryOpen(out) : null, out }
}

/** The condition the bundle gives each cell: its sample's. */
function conditionsOf(part) {
  const sample = new Uint16Array(part.files['sample.u16'].slice().buffer)
  return Array.from(sample, s => part.meta.samples[s]?.condition ?? null)
}

// ---- reading and editing the fixtures --------------------------------------
await h5wasm.ready

/** A categorical obs column's level per cell, read straight from the file. */
function obsLevels(path, name) {
  const f = new h5wasm.File(path, 'r')
  try {
    const g = f.get(`obs/${name}`)
    const cats = g.get('categories').value
    return Array.from(g.get('codes').value, c => cats[Number(c)])
  } finally { f.close() }
}

/** A copy of a committed fixture, edited through h5wasm. */
function copyWith(src, name, edit) {
  const dst = join(tmp, name)
  copyFileSync(src, dst)
  const f = new h5wasm.File(dst, 'a')
  try { edit(f) } finally { f.flush(); f.close() }
  return dst
}
const enc = (o, type, version) => {
  o.create_attribute('encoding-type', type)
  o.create_attribute('encoding-version', version)
  return o
}
const dict = (g, name) => enc(g.create_group(name), 'dict', '0.1.0')
const str = (g, name, value) => enc(g.create_dataset({ name, data: value }), 'string', '0.2.0')
const int = (g, name, value) =>
  enc(g.create_dataset({ name, data: BigInt64Array.from([BigInt(value)]), shape: [] }), 'numeric-scalar', '0.2.0')

/** Add a categorical obs column, in anndata's layout, after the existing ones. */
function addCategorical(f, name, levels, values) {
  const obs = f.get('obs')
  const g = enc(obs.create_group(name), 'categorical', '0.2.0')
  g.create_attribute('ordered', new Int8Array([0]), [], '<b')
  g.create_dataset({ name: 'categories', data: levels })
  g.create_dataset({ name: 'codes', data: Int8Array.from(values, v => levels.indexOf(v)) })
  const order = obs.attrs['column-order'].value
  obs.delete_attribute('column-order')
  obs.create_attribute('column-order', [...order, name])
}

// ---- Scans built in memory ---------------------------------------------------

/** Cell-major sparse from dense rows. */
function fromDense(rows, nGenes) {
  const indptr = [0], indices = [], data = []
  for (const r of rows) {
    for (let g = 0; g < nGenes; g++) if (r[g] !== 0) { indices.push(g); data.push(r[g]) }
    indptr.push(indices.length)
  }
  return {
    nCells: rows.length, nGenes,
    indptr: Int32Array.from(indptr), indices: Int32Array.from(indices), data: Float32Array.from(data),
  }
}

// Eight cells in three samples. s1 holds both conditions, and its first cell
// is `treat`, so first-cell relabelling puts all four of its cells in treat.
// s2 is all ctrl and s3 all treat.
const COUNTS = [[1, 0], [2, 1], [0, 3], [4, 0], [1, 1], [0, 2], [3, 0], [1, 4]]
const OWN = ['treat', 'ctrl', 'ctrl', 'treat', 'ctrl', 'ctrl', 'treat', 'treat']
const xy = (n, f) => Float32Array.from({ length: 2 * n }, (_v, i) => f(i))

function memScan(over = {}) {
  const n = over.n ?? 8
  const rows = over.rows ?? COUNTS
  return {
    format: 'h5ad', source: 'synthetic.h5ad', nCells: n,
    columns: over.columns ?? [
      categorical('cell_type', ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']),
      categorical('sample', ['s1', 's1', 's1', 's1', 's2', 's2', 's3', 's3']),
      fromFactor('cond', Int32Array.from(OWN, v => (v === 'ctrl' ? 0 : 1)), ['ctrl', 'treat']),
    ],
    matrices: [{ key: 'X', geneSet: 'var', nGenes: 2, nCells: n, kind: 'counts', load: () => fromDense(rows, 2) }],
    embeddings: over.embeddings ?? [{ key: 'umap', nDims: 2, load: () => xy(n, i => i) }],
    geneSets: { var: ['GeneA', 'GeneB'] },
    geneAliases: {},
    provenance: {},
    notes: [],
  }
}
const CHOICES = { cluster: 'cell_type', sample: 'sample', condition: 'cond', embedding: 'umap', label: 'syn' }

/** A bundle built in memory, read back. */
function openBuilt(res) {
  const files = unzipSync(res.zip)
  return { files, meta: JSON.parse(strFromU8(files['meta.json'])) }
}

/** One part of `parent`: the given cells, its matrix cut to them. */
function partOf(parent, key, cells, rows = COUNTS) {
  const mat = fromDense(cells.map(i => rows[i]), 2)
  return shard.shardScan(parent, { key, cells: Int32Array.from(cells), nnz: 0 }, new Map([['X', mat]]), `syn ${key}`)
}

// ===========================================================================
console.log('\nA SAMPLE IN TWO CONDITIONS IS SPLIT, AND EVERY CELL KEEPS ITS OWN')
{
  let b = null
  attempt('a crossed object builds', () => { b = openBuilt(build.buildBundle(memScan(), CHOICES)); return b.meta.nCells }, 8)
  attempt('every cell keeps its own condition', () => conditionsOf(b), OWN)
  attempt('the sample that spans is split, one sample per condition, named sample|condition',
    () => b.meta.samples, [
      { id: 's1|ctrl', condition: 'ctrl' }, { id: 's1|treat', condition: 'treat' },
      { id: 's2', condition: 'ctrl' }, { id: 's3', condition: 'treat' },
    ])
  attempt('the conditions are the object\'s, in its order', () => b.meta.conditions, ['ctrl', 'treat'])
  attempt('pseudobulk sums each split sample apart: no column under the unsplit name',
    () => strFromU8(b.files['pseudobulk.tsv']).split('\n')[0].split('\t').slice(1)
      .map(h => h.split('||')[0]).filter((s, i, a) => a.indexOf(s) === i),
    ['s1|ctrl', 's1|treat', 's2', 's3'])
  attempt('a note names the split sample and the count before and after',
    () => b.meta.notes.filter(s => /\bs1\b/.test(s) && /3 samples became 4/.test(s)).length, 1)

  const scanUnsplit = memScan({ columns: [
    categorical('cell_type', ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']),
    categorical('sample', ['s1', 's1', 's2', 's2', 's3', 's3', 's4', 's4']),
    fromFactor('cond', Int32Array.from([0, 0, 0, 0, 1, 1, 1, 1]), ['ctrl', 'treat']),
  ] })
  attempt('an object whose samples each sit in one condition comes back untouched (the same Scan)',
    () => build.oneConditionPerSample(scanUnsplit, CHOICES) === scanUnsplit, true)

  // A split name that is already a sample's is not reused for another.
  const clash = memScan({ columns: [
    categorical('cell_type', ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']),
    categorical('sample', ['s1', 's1', 's1', 's1', 's1|treat', 's1|treat', 's2', 's2']),
    fromFactor('cond', Int32Array.from([1, 0, 0, 1, 0, 0, 1, 1]), ['ctrl', 'treat']),
  ] })
  let c = null
  attempt('a split name already taken by another sample is not reused',
    () => { c = openBuilt(build.buildBundle(clash, CHOICES)); return c.meta.samples.map(s => s.id) },
    ['s1|ctrl', 's1|treat#2', 's1|treat', 's2'])
  attempt('and every cell still keeps its own condition',
    () => conditionsOf(c), ['treat', 'ctrl', 'ctrl', 'treat', 'ctrl', 'ctrl', 'treat', 'treat'])
}

console.log('\nWITH NO SAMPLE COLUMN, THE CELLS OF EACH CONDITION ARE ONE SAMPLE')
{
  let b = null
  attempt('it builds', () => { b = openBuilt(build.buildBundle(memScan(), { ...CHOICES, sample: null })); return b.meta.nCells }, 8)
  attempt('every cell keeps its own condition', () => conditionsOf(b), OWN)
  attempt('one sample per condition', () => b.meta.samples,
    [{ id: 'all cells|ctrl', condition: 'ctrl' }, { id: 'all cells|treat', condition: 'treat' }])
  attempt('and the note says what a sample is here',
    () => b.meta.notes.filter(s => /no sample column/.test(s) && /each condition/.test(s)).length, 1)
  attempt('with neither column it is one sample in one group, as before',
    () => openBuilt(build.buildBundle(memScan(), { ...CHOICES, sample: null, condition: null })).meta.samples,
    [{ id: 'all cells', condition: 'all cells' }])
}

console.log('\nTHE SPLIT IS DECIDED ON THE WHOLE OBJECT, SO EVERY PART NAMES A SAMPLE THE SAME WAY')
{
  // Part A holds only s1's treat cell; part B holds both of its conditions.
  // Decided part by part, A would call it "s1" and B "s1|ctrl" and "s1|treat"
  // — and the studio, which merges parts by sample id, would count one animal
  // as two replicates.
  let parts = null
  attempt('a parent resolved before the split builds both parts', () => {
    const parent = build.oneConditionPerSample(memScan(), CHOICES)
    parts = [[0, 4, 5], [1, 2, 3, 6, 7]].map((cells, k) => ({
      cells, b: openBuilt(build.buildBundle(partOf(parent, `p${k}`, cells), CHOICES)),
    }))
    return parts.length
  }, 2)
  attempt('each part names its samples as the whole object does',
    () => parts.map(p => p.b.meta.samples.map(s => `${s.id}:${s.condition}`)),
    [['s1|treat:treat', 's2:ctrl'], ['s1|ctrl:ctrl', 's1|treat:treat', 's3:treat']])
  attempt('and every cell in every part keeps its own condition',
    () => parts.every(p => JSON.stringify(conditionsOf(p.b)) === JSON.stringify(p.cells.map(i => OWN[i]))), true)
  attempt('the note is in every part, and counts the whole object',
    () => parts.map(p => p.b.meta.notes.filter(s => /3 samples became 4/.test(s)).length), [1, 1])
}

console.log('\nTHE APP\'S WARNING SAYS WHAT NOW HAPPENS')
{
  const s = memScan()
  attempt('it no longer promises to keep only the first condition',
    () => /only the first/.test(scanMod.crossCheck(s.columns[1], s.columns[2])), false)
  attempt('it says the sample is split, and names the new samples',
    () => { const m = scanMod.crossCheck(s.columns[1], s.columns[2]); return [/appears in both/.test(m), /split/.test(m), m.includes('s1|ctrl')] },
    [true, true, true])
}

console.log('\nTHROUGH THE COMMAND: A CROSSED OBJECT KEEPS EVERY CELL\'S CONDITION')
{
  // A copy of the family fixture with one more obs column, `crossed`, that
  // varies inside s1 and s3 and is constant inside every other sample.
  const sampleOf = obsLevels(FIX.classic, 'sample')
  const at = new Map()
  const crossedOf = sampleOf.map(s => {
    const k = at.get(s) ?? 0
    at.set(s, k + 1)
    if (s === 's1') return k % 2 === 0 ? 'c2' : 'c1'
    if (s === 's2') return 'c2'
    if (s === 's3') return k < 10 ? 'c1' : 'c2'
    return 'c1'
  })
  const crossedFile = copyWith(FIX.classic, 'crossed.h5ad', f => addCategorical(f, 'crossed', ['c1', 'c2'], crossedOf))
  const a = convert(crossedFile, 'crossed.zip', ['--condition', 'crossed', '--label', 'crossed'])
  check('it exits 0', a.r.code, 0)
  const p = a.o?.parts[0]
  attempt('every cell keeps its own condition', () => {
    const got = conditionsOf(p)
    return got.filter((c, i) => c !== crossedOf[i]).length
  }, 0)
  attempt('s1 and s3 are split, and only they are renamed',
    () => p.meta.samples.map(s => `${s.id}:${s.condition}`),
    ['s1|c1:c1', 's1|c2:c2', 's2:c2', 's3|c1:c1', 's3|c2:c2', 's4:c1', 's5:c1', 's6:c1', 's7:c1', 's8:c1'])
  attempt('a note names both split samples and the count',
    () => p.meta.notes.filter(s => /\bs1\b/.test(s) && /\bs3\b/.test(s) && /8 samples became 10/.test(s)).length, 1)
  attempt('and the command says it on stderr', () => /^note: .*\bs1\b.*\bs3\b.*8 samples became 10/m.test(a.r.err), true)
  const b = convert(crossedFile, 'crossed-again.zip', ['--condition', 'crossed', '--label', 'crossed'])
  attempt('converted twice, the same bytes', () => [b.r.code, sha256(readFileSync(b.out)) === sha256(a.o.bytes)], [0, true])
}

console.log('\nTHROUGH THE COMMAND: NO SAMPLE COLUMN, EVERY CELL KEEPS ITS CONDITION')
{
  const a = convert(FIX.classic, 'nosample.zip', ['--sample', 'none', '--condition', 'factor_a', '--label', 'nosample'])
  check('it exits 0', a.r.code, 0)
  const own = obsLevels(FIX.classic, 'factor_a')
  attempt('every cell keeps its own factor_a level',
    () => conditionsOf(a.o.parts[0]).filter((c, i) => c !== own[i]).length, 0)
  attempt('one sample per level', () => a.o.parts[0].meta.samples,
    [{ id: 'all cells|a1', condition: 'a1' }, { id: 'all cells|a2', condition: 'a2' }])
}

// ===========================================================================
// The committed fixture, every role guessed: the bundle the gate diffs.
const G = convert(FIX.classic, 'guessed.zip', ['--label', 'family'])
const g = G.o?.parts[0]

console.log('\nTHE RUN\'S RECORD TRAVELS WITH THE BUNDLE')
{
  // What scripts/make-family-fixture.py writes as uns['scintegrate']['provenance'].
  const PROVENANCE = {
    tool: 'scintegrate', version: '0.0.0+fixture', commit: '0'.repeat(40), state_version: 1,
    argv: ['scintegrate', 'integrate', '--h5ad', 'joint.h5ad'], python: '3.12', job: 'fixture',
    wrapped_versions: { anndata: '0.13.2' }, sees: { none: [], harmony: ['sample'] }, input: 'joint.h5ad',
    batch_key: 'sample', label_key: 'cell_type', seed: 0,
  }
  check('the guessed conversion exits 0', G.r.code, 0)
  attempt('meta.run lists the one record the object has, and where it was',
    () => g.meta.run.records.map(r => r.path), ['uns/scintegrate/provenance'])
  attempt('and the record is the provenance itself: the tool, its version and commit, the job, argv, the seed',
    () => canon(g.meta.run.records[0].record), canon(PROVENANCE))
  attempt('meta.run says where it looked', () => typeof g.meta.run.searched === 'string' && /provenance/.test(g.meta.run.searched) && /harness_run/.test(g.meta.run.searched), true)
  attempt('a note names the record carried', () => g.meta.notes.filter(s => s.includes('uns/scintegrate/provenance')).length, 1)

  const d = convert(FIX.dflt, 'default.zip', ['--label', 'family'])
  attempt("anndata's default encoding gives the same run", () => [d.r.code, d.o.parts[0].meta.run.records.length, canon(d.o.parts[0].meta.run) === canon(g.meta.run)], [0, 1, true])

  const m = convert(FIX.mask, 'norecord.zip', ['--label', 'norecord'])
  attempt('an object whose uns holds no record says so: an empty list', () => [m.r.code, m.o.parts[0].meta.run.records], [0, []])
  attempt('and a note says the bundle cannot tell which run made it',
    () => m.o.parts[0].meta.notes.filter(s => /no run record/.test(s)).length, 1)

  // Every form the rule names, and two it does not, in one copy. uns lists
  // its entries by name, so the records come in that order.
  const many = copyWith(FIX.classic, 'records.h5ad', f => {
    const uns = f.get('uns')
    const an = dict(uns, 'analysis')
    const prov = dict(an, 'provenance')
    str(prov, 'tool', 'synthetic-tool')
    int(prov, 'n', 3)
    // A real h5py boolean (an HDF5 enum), linked from the fixture's own uns.
    prov.create_hard_link('/uns/scintegrate/assessment/summary/needed', 'dirty')
    const ct = dict(uns, 'cell_type_provenance')
    str(ct, 'tree', 'synthetic')
    int(ct, 'n_cells', 240)
    const deep = dict(dict(uns, 'deep'), 'inner')
    str(dict(deep, 'provenance'), 'x', 'two levels down: not a record by the rule')
    str(uns, 'harness_run', JSON.stringify({ stages: [{ stage: 'integration', tool: 'synthetic', run_key: 'rk-1', job: 'j-1', sealed: true }] }))
    const qc = dict(uns, 'qc_decl')
    str(qc, 'schema', 'synthetic/provenance@1')
    str(qc, 'run_key', 'rk-0')
    str(qc, 'commit', 'abc1234')
    const other = dict(uns, 'unrelated')
    str(other, 'schema', 'synthetic/results@1')
  })
  const r = convert(many, 'records.zip', ['--label', 'records'])
  const run = r.o?.parts[0].meta.run
  attempt('every record the rule names is carried, in uns order, and nothing else',
    () => run.records.map(x => x.path),
    ['uns/analysis/provenance', 'uns/cell_type_provenance', 'uns/harness_run', 'uns/qc_decl', 'uns/scintegrate/provenance'])
  const rec = (p) => run.records.find(x => x.path === p)?.record
  attempt('a boolean stays a boolean, an int64 a number', () => [rec('uns/analysis/provenance').dirty, rec('uns/analysis/provenance').n, rec('uns/cell_type_provenance').n_cells], [true, 3, 240])
  attempt("the harness's run record is parsed from its JSON",
    () => rec('uns/harness_run'), { stages: [{ stage: 'integration', tool: 'synthetic', run_key: 'rk-1', job: 'j-1', sealed: true }] })
  attempt('a record that names provenance in its own schema is carried whole',
    () => canon(rec('uns/qc_decl')), canon({ commit: 'abc1234', run_key: 'rk-0', schema: 'synthetic/provenance@1' }))

  const bad = copyWith(FIX.classic, 'bad-run.h5ad', f => { str(f.get('uns'), 'harness_run', 'not json {') })
  const br = convert(bad, 'bad-run.zip', ['--label', 'bad'])
  const badRun = br.o?.parts[0].meta.run
  attempt('a harness_run that is not JSON is carried as its text, not dropped',
    () => badRun.records.find(x => x.path === 'uns/harness_run')?.record, 'not json {')
  attempt('and a note says it is not JSON',
    () => br.o.parts[0].meta.notes.filter(s => /harness_run/.test(s) && /not JSON/i.test(s)).length, 1)
}

console.log('\nA RECORD IS COPIED AS FAR AS JSON CAN HOLD IT, AND THE REST IS NAMED IN PLACE')
{
  // scripts/fixtures/uns-records.h5ad, written by anndata and h5py (see
  // make-family-fixture.py): its integration record holds a structured array
  // with an int64 field, which crashed the whole conversion once ("Do not know
  // how to serialize a BigInt", exit 1, no file); its other record holds
  // leaves a JSON copy cannot hold as they are.
  const u = convert(FIX.records, 'uns-records.zip', ['--label', 'records'])
  check('an object whose provenance holds an int64 rec-array converts', u.r.code, 0)
  const rec = (p) => u.o?.parts[0].meta.run.records.find(x => x.path === p)?.record
  attempt('the rec-array is copied row by row, under its field names, int64 as a number',
    () => rec('uns/scintegrate/provenance').ranked,
    [{ method: 'harmony', rank: 3, score: 0.6 }, { method: 'none', rank: 1, score: 0.5 }])
  attempt('a float16 leaf is named in place, and the rest of its record is still copied',
    () => { const o = rec('uns/odd_provenance'); return [/^\(not copied: .*float16/.test(o.half), o.tool, o.flag] },
    [true, 'synthetic-tool', false])
  attempt('an HDF5 reference is named in place, not copied as bytes',
    () => /^\(not copied: .*HDF5 reference/.test(rec('uns/odd_provenance').ref), true)
  attempt('uint64 past 2**53 is its exact digits, and NaN and the infinities are words',
    () => { const o = rec('uns/odd_provenance'); return [o.big, o.special] },
    ['9223372036854775813', ['NaN', 'Infinity', '-Infinity', 1]])
  attempt('a dataframe is named in place', () => /^\(not copied: .*dataframe/.test(rec('uns/odd_provenance').frame), true)
}

console.log('\nA RECORD PAST ITS BUDGET IS NAMED WITH ITS SIZE, NOT COPIED INTO EVERY PART')
{
  // The run section is written into every part of a split object, so what one
  // record and the whole section may take is bounded, in bytes as meta.json
  // writes them: 262,144 for a record, 1,048,576 for the section.
  const big = copyWith(FIX.classic, 'big-record.h5ad', f => {
    const p = dict(f.get('uns'), 'log_provenance')
    str(p, 'log', 'x'.repeat(400_000))
    str(p, 'tool', 'synthetic-tool')
  })
  const b = convert(big, 'big-record.zip', ['--label', 'big'])
  const run = b.o?.parts[0].meta.run
  attempt('a record over 262,144 bytes as written is named in place, with its size',
    () => { const r = run.records.find(x => x.path === 'uns/log_provenance').record; return [typeof r, /400,0\d\d bytes as written/.test(r), /262,144/.test(r)] },
    ['string', true, true])
  attempt('the record under the budget is still copied',
    () => run.records.find(x => x.path === 'uns/scintegrate/provenance').record.tool, 'scintegrate')
  attempt('a note names the record left out, and why',
    () => b.o.parts[0].meta.notes.filter(s => s.includes('uns/log_provenance') && /budget/.test(s)).length, 1)
  attempt('so meta.json stays small', () => b.o.parts[0].files['meta.json'].length < 20_000, true)

  // Six records of about 200,000 bytes each, every one under the record budget:
  // after the fixture's own record and five of them, the sixth would pass the
  // section's budget, and is left out.
  const many = copyWith(FIX.classic, 'many-records.h5ad', f => {
    const uns = f.get('uns')
    for (let i = 1; i <= 6; i++) str(dict(uns, `z${i}_provenance`), 'log', 'y'.repeat(200_000))
  })
  const m = convert(many, 'many-records.zip', ['--label', 'many'])
  attempt('records are copied until the section would pass 1,048,576 bytes; the rest are left out, in one entry',
    () => { const r = m.o.parts[0].meta.run; return [r.records.map(x => `${x.path}:${typeof x.record === 'string' ? 'named' : 'copied'}`), r.leftOut] },
    [['uns/scintegrate/provenance:copied', 'uns/z1_provenance:copied', 'uns/z2_provenance:copied',
      'uns/z3_provenance:copied', 'uns/z4_provenance:copied', 'uns/z5_provenance:copied'], { count: 1, paths: ['uns/z6_provenance'] }])
  attempt('so meta.json stays under the section budget and a little more', () => m.o.parts[0].files['meta.json'].length < 1_100_000, true)
}

console.log('\nTHE BUDGETS COUNT WHAT meta.json WRITES, AND THE SECTION FILLS FAIRLY')
{
  /** int8 arrays of zeros, in anndata's array encoding. */
  const zeros = (g, name) => enc(g.create_dataset({ name, data: new Int8Array(10_000) }), 'array', '0.2.0')

  // Four records of 13 x 10,000 zeros: about 260,000 characters of compact
  // JSON each, under the budget, but meta.json writes one number a line, and
  // the four of them came to 4.7 MB.
  const flat = copyWith(FIX.classic, 'written.h5ad', f => {
    for (let r = 1; r <= 4; r++) {
      const p = dict(f.get('uns'), `w${r}_provenance`)
      for (let i = 0; i < 13; i++) zeros(p, `z${String(i).padStart(2, '0')}`)
    }
  })
  const w = convert(flat, 'written.zip', ['--label', 'written'])
  attempt('a record under the budget as compact JSON but over it as written is named with its written size',
    () => w.o.parts[0].meta.run.records.filter(x => /^w\d/.test(x.path.slice(4)))
      .map(x => /^\(not copied: [\d,]+ bytes as written, over the 262,144-byte budget/.test(x.record)),
    [true, true, true, true])
  attempt('so meta.json stays small', () => w.o.parts[0].files['meta.json'].length < 20_000, true)

  // A record forty groups deep: read to sixteen, named below that.
  const deep = copyWith(FIX.classic, 'deep.h5ad', f => {
    let g = dict(f.get('uns'), 'deep_provenance')
    for (let i = 0; i < 40; i++) g = dict(g, 'a')
    str(g, 'leaf', 'x')
  })
  const d = convert(deep, 'deep.zip', ['--label', 'deep'])
  attempt('a record is read sixteen groups deep, and what lies deeper is named in place',
    () => {
      let at = d.o.parts[0].meta.run.records.find(x => x.path === 'uns/deep_provenance').record
      for (let i = 0; i < 16; i++) at = at.a
      return [typeof at, /^\(not copied: .*nested more than 16 groups below its record/.test(at.a)]
    }, ['object', true])

  // 400 records of 3,000 bytes: about a sixth do not fit.
  const crowd = copyWith(FIX.classic, 'crowd.h5ad', f => {
    const uns = f.get('uns')
    for (let i = 0; i < 400; i++) str(uns, `r${String(i).padStart(3, '0')}_provenance`, 'x'.repeat(3_000))
  })
  const c = convert(crowd, 'crowd.zip', ['--label', 'crowd'])
  attempt('once the section is full, one entry counts every record left out, and no record is written as "not copied"',
    () => {
      const run = c.o.parts[0].meta.run
      return [run.records.length - 1 + run.leftOut.count, run.leftOut.paths.length,
        run.records.filter(x => typeof x.record === 'string' && x.record.startsWith('(not copied')).length]
    }, [400, 10, 0])
  attempt('and one note says so', () => c.o.parts[0].meta.notes.filter(s => /left out of run/.test(s)).length, 1)

  // Five records of 250,002 bytes, each under the record budget, all before
  // the tool's own record in uns order: filled in uns order, they pushed it out.
  const first = copyWith(FIX.classic, 'first.h5ad', f => {
    for (let i = 1; i <= 5; i++) str(f.get('uns'), `a${i}_provenance`, 'q'.repeat(250_000))
  })
  const o = convert(first, 'first.zip', ['--label', 'first'])
  attempt("a tool's own provenance is written before *_provenance records that come earlier in uns",
    () => { const run = o.o.parts[0].meta.run; return [run.records.find(x => x.path === 'uns/scintegrate/provenance')?.record?.tool, run.leftOut] },
    ['scintegrate', { count: 1, paths: ['uns/a5_provenance'] }])
}

console.log('\nHARD LINKS CANNOT MAKE A SMALL RECORD AN ENDLESS READ')
{
  // Eighteen groups, each holding two hard links to the next: 2^17 paths in a
  // small file, which read for minutes without finishing. 30 s is many times
  // what the fixed reader needs; the old one is killed there.
  const dag = copyWith(FIX.classic, 'dag.h5ad', f => {
    const store = f.create_group('dagstore')
    const levels = Array.from({ length: 18 }, (_v, i) => dict(store, `L${i}`))
    str(levels[17], 'leaf', 'x')
    for (let i = 0; i < 17; i++) for (const c of ['c0', 'c1']) levels[i].create_hard_link(`/dagstore/L${i + 1}`, c)
    f.get('uns').create_hard_link('/dagstore/L0', 'dag_provenance')
  })
  const g = convert(dag, 'dag.zip', ['--label', 'dag'], 30_000)
  attempt('a record whose groups are reached by many hard links converts, in seconds', () => [g.r.code, g.r.seconds < 30], [0, true])
  const rec = () => g.o.parts[0].meta.run.records.find(x => x.path === 'uns/dag_provenance').record
  attempt('each group is copied once; another link to it is named in place, not read again',
    () => [typeof rec().c0, /^\(not copied: uns\/dag_provenance\/c1 is uns\/dag_provenance\/c0 again, reached by another link/.test(rec().c1)],
    ['object', true])
  attempt('and a note says so', () => g.o.parts[0].meta.notes.filter(s => /uns\/dag_provenance/.test(s) && /more than one link/.test(s)).length, 1)

  // A group holding three hard links to itself.
  const loop = copyWith(FIX.classic, 'loop.h5ad', f => {
    const p = dict(f.get('uns'), 'loop_provenance')
    str(p, 'tool', 'synthetic-tool')
    for (const c of ['self0', 'self1', 'self2']) p.create_hard_link('/uns/loop_provenance', c)
  })
  const l = convert(loop, 'loop.zip', ['--label', 'loop'], 30_000)
  attempt('a group that links to itself converts, in seconds', () => [l.r.code, l.r.seconds < 30], [0, true])
  attempt('each link back to it is named in place, and the rest of the record copied',
    () => {
      const r = l.o.parts[0].meta.run.records.find(x => x.path === 'uns/loop_provenance').record
      return [r.tool, ...['self0', 'self1', 'self2'].map(k => String(r[k]).startsWith(`(not copied: uns/loop_provenance/${k} is uns/loop_provenance again`))]
    }, ['synthetic-tool', true, true, true])

  // 10,050 datasets in one record: read to 10,000 groups and datasets, the
  // rest named in one sentence.
  const wide = copyWith(FIX.classic, 'wide.h5ad', f => {
    const p = dict(f.get('uns'), 'wide_provenance')
    for (let i = 0; i < 10_050; i++) p.create_dataset({ name: `d${String(i).padStart(5, '0')}`, data: 'x' })
  })
  const w = convert(wide, 'wide.zip', ['--label', 'wide'], 60_000)
  attempt('a record is read to 10,000 groups and datasets, and the rest of it is named in one sentence',
    () => {
      const r = w.o.parts[0].meta.run.records.find(x => x.path === 'uns/wide_provenance').record
      return [Object.keys(r).length, /^51 more entries of uns\/wide_provenance, past the 10,000 groups and datasets/.test(r['(not copied)'])]
    }, [10_000, true])
  attempt('and a note says so', () => w.o.parts[0].meta.notes.filter(s => /uns\/wide_provenance/.test(s) && /read only to the 10,000 groups and datasets one record may take/.test(s)).length, 1)
}

console.log('\nTHE RUN SECTION IS READ IN ONE WALK: AN OBJECT IS READ ONCE, WHATEVER REACHES IT')
{
  // 400 records, half hard links and half soft links to one group of 9,990
  // datasets. Read record by record, every link read the whole group again,
  // about two seconds each; the old reader is killed at 60 s.
  const shared = copyWith(FIX.classic, 'shared.h5ad', f => {
    const g = dict(f, 'shared')
    for (let i = 0; i < 9_990; i++) g.create_dataset({ name: `k${String(i).padStart(5, '0')}`, data: 'v' })
    const uns = f.get('uns')
    for (let i = 0; i < 400; i++) {
      const name = `r${String(i).padStart(4, '0')}_provenance`
      if (i % 2) uns.create_soft_link('/shared', name)
      else uns.create_hard_link('/shared', name)
    }
  })
  const s = convert(shared, 'shared.zip', ['--label', 'shared'], 60_000)
  attempt('400 records linked to one large group convert, in seconds', () => [s.r.code, s.r.seconds < 30], [0, true])
  // Each record is read on its own, so five copy the group whole and the
  // sixth reaches the section's 50,000 groups and datasets part-way through;
  // the 394 after it are left out, unread.
  attempt('the section stops reading at 50,000 groups and datasets, and leaves the rest out, unread',
    () => {
      const run = s.o.parts[0].meta.run
      const recs = run.records.filter(x => /^uns\/r\d{4}_provenance$/.test(x.path))
      return [recs.length, recs.slice(0, 5).map(x => Object.keys(x.record).length), /past the 50,000 groups and datasets read for the run section$/.test(recs[5].record['(not copied)']), run.leftOut.count]
    }, [6, [9_990, 9_990, 9_990, 9_990, 9_990], true, 394])
  attempt('and the note says 394 were not read', () => s.o.parts[0].meta.notes.some(n => /394 of them were not read/.test(n)), true)

  // The section's own budget, on groups shaped like h5wasm's: 60 records of
  // 1,000 datasets each, every one its own object. After 50,000 groups and
  // datasets, the records not yet read are left out, unread.
  let id = 0
  const ref = () => { const k = ++id; return () => ({ ref_data: Uint8Array.from([k & 255, (k >> 8) & 255, (k >> 16) & 255]) }) }
  const ds = () => ({ shape: [], value: 'v', metadata: { type: 3, size: 1 }, attrs: {}, create_reference: ref() })
  const grp = (children) => ({ attrs: {}, keys: () => Object.keys(children), get: (k) => children[k] ?? null, create_reference: ref() })
  const uns = grp(Object.fromEntries(Array.from({ length: 60 }, (_v, r) =>
    [`r${String(r).padStart(2, '0')}_provenance`, grp(Object.fromEntries(Array.from({ length: 1_000 }, (_w, i) => [`d${String(i).padStart(3, '0')}`, ds()])))])))
  const h5 = await tryImport('../src/lib/h5ad.ts')
  const said = []
  attempt('after 50,000 groups and datasets for the section, the records not yet reached are left out, unread',
    () => { const run = h5.readRun(uns, s => said.push(s)); return [run.records.length, run.leftOut.count, run.leftOut.paths[0]] },
    [50, 10, 'uns/r50_provenance'])
  attempt('and the note says how many were not read', () => said.filter(s => /left out of run/.test(s) && /10 of them were not read/.test(s)).length, 1)
}

console.log('\nAN OBJECT ONE RECORD READ IS NOT TAKEN FROM ANOTHER')
{
  // A holds a link to the tool's group and is over its budget, so its content
  // is never written; B is the tool's group itself. Remembered across records,
  // B became a pointer into A, and the tool's record was lost.
  const toobig = copyWith(FIX.classic, 'toobig.h5ad', f => {
    const t = dict(f, 'toolrec')
    str(t, 'tool', 'B-tool')
    str(t, 'version', '9.9')
    str(t, 'commit', 'deadbeef')
    const a = dict(f.get('uns'), 'a_provenance')
    str(a, 'tool', 'A-tool')
    for (const [name, k] of [['big', 1], ['big2', 2]]) {
      enc(a.create_dataset({ name, data: Float64Array.from({ length: 9_000 }, (_v, i) => Math.sin(i * k + 0.5) * 1e6) }), 'array', '0.2.0')
    }
    a.create_hard_link('/toolrec', 'shared')
    f.get('uns').create_hard_link('/toolrec', 'b_provenance')
  })
  const t = convert(toobig, 'toobig.zip', ['--label', 'toobig'])
  attempt("a record that reaches what an unwritten record reached is still copied: the tool's record is not lost",
    () => {
      const run = t.o.parts[0].meta.run
      return [/^\(not copied: [\d,]+ bytes as written/.test(run.records.find(x => x.path === 'uns/a_provenance').record),
        run.records.find(x => x.path === 'uns/b_provenance').record]
    }, [true, { commit: 'deadbeef', tool: 'B-tool', version: '9.9' }])

  // Another file, a copy of the fixture with one more field in its
  // integration record: the record's group sits at the same address in both
  // files. Linked from this file, it was taken for this file's own.
  const other = copyWith(FIX.classic, 'other.h5ad', f => str(f.get('uns/scintegrate/provenance'), 'origin', 'OTHER-FILE'))
  const ext = copyWith(FIX.classic, 'extcollide.h5ad', f => {
    const uns = f.get('uns')
    uns.create_external_link(other, '/uns/scintegrate/provenance', 'extp_provenance')
    const m = dict(uns, 'mix_provenance')
    m.create_hard_link('/uns/scintegrate/provenance', 'local')
    m.create_external_link(other, '/uns/scintegrate/provenance', 'remote')
  })
  const e = convert(ext, 'extcollide.zip', ['--label', 'extcollide'])
  const rec = (p) => e.o.parts[0].meta.run.records.find(x => x.path === p)?.record
  attempt('a record in another file, at the same address as one here, is copied as its own',
    () => [rec('uns/extp_provenance')?.origin, rec('uns/scintegrate/provenance')?.origin], ['OTHER-FILE', undefined])
  attempt('and inside one record, an object here and one in another file at the same address are both copied',
    () => [rec('uns/mix_provenance')?.local?.tool, rec('uns/mix_provenance')?.remote?.origin], ['scintegrate', 'OTHER-FILE'])
}

console.log("\nWITHIN A KIND, FIRST FIT IN uns ORDER; WHAT IS LEFT OUT IS NAMED LARGEST FIRST")
{
  // The annotator's record (250,000 bytes) comes before four 200,000-byte
  // *_provenance records in uns. Filled smaller first, the four took the
  // section and the annotator's record was left out.
  const smaller = copyWith(FIX.classic, 'smaller.h5ad', f => {
    const uns = f.get('uns')
    for (let i = 0; i < 4; i++) str(uns, `z${i}_provenance`, 'x'.repeat(200_000))
    const a = dict(uns, 'cell_type_provenance')
    str(a, 'tool', 'annotator')
    str(a, 'version', '1.0')
    str(a, 'detail', 'y'.repeat(250_000))
  })
  const s = convert(smaller, 'smaller.zip', ['--label', 'smaller'])
  attempt('a larger record earlier in uns is not left out for smaller ones after it',
    () => { const run = s.o.parts[0].meta.run; return [run.records.find(x => x.path === 'uns/cell_type_provenance')?.record?.tool, run.leftOut] },
    ['annotator', { count: 1, paths: ['uns/z3_provenance'] }])

  // A hundred 12,000-byte `provenance` records earlier in uns than the
  // integration stage's own, made 50,000 bytes here: they fill the section, the
  // tool's record is left out, and it must be named. Listed in fill order,
  // the first ten left out were all the small ones.
  const flood = copyWith(FIX.classic, 'flood.h5ad', f => {
    const uns = f.get('uns')
    for (let i = 0; i < 100; i++) str(dict(uns, `g${String(i).padStart(3, '0')}`), 'provenance', 'x'.repeat(12_000))
    str(f.get('uns/scintegrate/provenance'), 'notes', 'n'.repeat(50_000))
  })
  const fl = convert(flood, 'flood.zip', ['--label', 'flood'])
  attempt("a tool's own record left out is named first among those left out",
    () => { const lo = fl.o.parts[0].meta.run.leftOut; return [lo.count, lo.paths[0]] }, [15, 'uns/scintegrate/provenance'])
  attempt('and the note names it first', () => fl.o.parts[0].meta.notes.some(n => /left out of run, largest first: uns\/scintegrate\/provenance,/.test(n)), true)
}

console.log('\nAN UNS ENTRY THAT CANNOT BE INSPECTED IS SAID, NOT FATAL')
{
  // scripts/fixtures/uns-unreadable.h5ad: a group whose encoding-type
  // attribute is float16, which h5wasm cannot read. Classifying it stopped the
  // conversion with exit 1.
  const u = convert(FIX.unreadable, 'unreadable.zip', ['--label', 'unreadable'])
  check('it converts', u.r.code, 0)
  attempt('a note names the entry and why it could not be inspected',
    () => u.o.parts[0].meta.notes.filter(s => s.includes('uns/plain_group') && /could not be inspected/.test(s) && /Float16/.test(s)).length, 1)
  attempt('and the record beside it is still copied',
    () => u.o.parts[0].meta.run.records.map(x => [x.path, x.record?.tool]), [['uns/x_provenance', 'synthetic-tool']])
}

// ===========================================================================
console.log('\nONLY TWO-COLUMN EMBEDDINGS ARE OFFERED, AND EACH ONLY ONCE')
{
  // The committed fixture's obsm, measured: X_harmony, X_none and X_pca are
  // 240 x 10; X_umap, X_umap_harmony and X_umap_none are 240 x 2. X_umap is a
  // byte copy of X_umap_harmony, and X_none of X_pca.
  attempt('the fixture carries X_umap and X_umap_none, the default first',
    () => g.meta.embeddings, [{ key: 'X_umap', file: 'embed.f32' }, { key: 'X_umap_none', file: 'embed.X_umap_none.f32' }])
  attempt('collection.json lists the same two', () => G.o.index.embeddings, ['X_umap', 'X_umap_none'])
  attempt('no file is written for what is not offered',
    () => Object.keys(g.files).filter(n => n.startsWith('embed.')).sort(), ['embed.X_umap_none.f32', 'embed.f32'])
  const skip = () => g.meta.notes.find(s => /X_harmony/.test(s) && /X_umap_harmony/.test(s)) ?? ''
  attempt('one note names each latent space with its width',
    () => ['X_harmony', 'X_none', 'X_pca'].map(k => new RegExp(`${k} \\(10 columns\\)`).test(skip())), [true, true, true])
  attempt('and the copy, with the embedding it copies', () => /X_umap_harmony[^;]*X_umap\b/.test(skip()), true)

  const h = convert(FIX.classic, 'emb-harmony-umap.zip', ['--label', 'family', '--embedding', 'X_umap_harmony'])
  attempt('choosing the copy keeps it, under its own name, and drops the other',
    () => h.o.parts[0].meta.embeddings, [{ key: 'X_umap_harmony', file: 'embed.f32' }, { key: 'X_umap_none', file: 'embed.X_umap_none.f32' }])
  attempt('and the note names X_umap as its copy',
    () => h.o.parts[0].meta.notes.some(s => /X_umap \(the same coordinates as X_umap_harmony/.test(s)), true)

  const p = convert(FIX.classic, 'emb-pca.zip', ['--label', 'family', '--embedding', 'X_pca'])
  const pm = p.o?.parts[0].meta
  attempt('choosing X_pca draws its first two components, flagged as before',
    () => [pm.embedding, pm.notes.some(s => /first two principal components/.test(s))], ['X_pca', true])
  attempt('and only the two-column ones travel beside it, each once',
    () => pm.embeddings, [
      { key: 'X_pca', file: 'embed.f32' }, { key: 'X_umap', file: 'embed.X_umap.f32' },
      { key: 'X_umap_none', file: 'embed.X_umap_none.f32' },
    ])

  const x = convert(FIX.classic, 'emb-harmony.zip', ['--label', 'family', '--embedding', 'X_harmony'])
  attempt('choosing a latent space that is not a PCA is flagged with its width',
    () => x.o.parts[0].meta.notes.some(s => /^X_harmony has 10 columns/.test(s)), true)
}

console.log('\nCOPIES ARE FOUND ON THE WHOLE OBJECT')
{
  const a = xy(8, i => i * 0.5)
  const b = Float32Array.from(a)
  const embs = (third) => [
    { key: 'umap', nDims: 2, load: () => xy(8, i => i) },
    { key: 'tsne_a', nDims: 2, load: () => a },
    { key: 'tsne_b', nDims: 2, load: () => third },
    { key: 'pca10', nDims: 10, load: () => xy(8, i => -i) },
  ]
  let one = null
  attempt('two copies that are not the default are offered once, under the first name',
    () => { one = openBuilt(build.buildBundle(memScan({ embeddings: embs(b) }), CHOICES)); return one.meta.embeddings.map(e => e.key) },
    ['umap', 'tsne_a'])
  attempt('and the note names both the copy and the latent space',
    () => one.meta.notes.some(s => /tsne_b \(the same coordinates as tsne_a\)/.test(s) && /pca10 \(10 columns\)/.test(s)), true)

  // Differs from tsne_a in the last cell only, which part 0 does not hold.
  const near = Float32Array.from(a)
  near[15] += 1
  const decide = (third) => {
    const parent = memScan({ embeddings: embs(third) })
    return [[0, 1, 2, 3], [4, 5, 6, 7]].map((cells, k) =>
      openBuilt(build.buildBundle(partOf(parent, `p${k}`, cells), { ...CHOICES, condition: null })).meta.embeddings.map(e => e.key))
  }
  attempt('two embeddings equal in one part\'s cells but not the whole object\'s are both kept, in every part',
    () => decide(near), [['umap', 'tsne_a', 'tsne_b'], ['umap', 'tsne_a', 'tsne_b']])
  attempt('a copy over the whole object is dropped from every part',
    () => decide(b), [['umap', 'tsne_a'], ['umap', 'tsne_a']])
}

// ===========================================================================
console.log('\nNOTHING ELSE IN THE FIXTURE\'S BUNDLE CHANGED (against 51a00fa)')
{
  // Measured on 51a00fa's conversion of the committed fixture (every role
  // guessed, --label family), sha256 189fde2d...: every entry this wave does
  // not touch, and every field and note it does not touch, hashed.
  const BEFORE = {
    'cluster.u16': '96c6c3ae1c4a743508b8d62595842b075304c18fcf331cbed2b1e69ea5dcdc62',
    'embed.X_umap_none.f32': '24aab30f2919bd7dcb9ce2fd7293141e478bc5e099bc89f2ec38f54028e4d44b',
    'embed.f32': '31c6f71e5e02390a8180b75a837369797b8bb568cba5a034f94aca73bd0eb590',
    'expr.chunk.bin': '08f496e7207a70be2394e988b0e031264c8f91c5b22e4a2ef66a42c3e7343ddc',
    'expr.chunkptr.i32': '8c06f6566490d1784e3b3172808ce596edc9bcb833cf6d92cc7bdec01069aea8',
    'expr.data.f32': '6934018a109bd9de0540231f0028d70467999b90f5aed9906e2fd39aeedb3723',
    'expr.indices.i32': 'a6f1bff0312d175708ef551f9ab9e33d7aee436d0c801854413b7b1d371f59b0',
    'expr.indptr.i32': 'b8b1549d12b564e43b66c40629dd76a582106a80a89ae942fde7c928effa818a',
    'extra.factor_a_factor_b.u16': '5d6ec6063dbd8476144bc11b5d64cc8032d39ae56dff6a230a03dfeb34c07c8f',
    'extra.factor_b.u16': 'dfeca1c10658b88d11df7ddd3f863a18f8cae46fd157d4177d249524d508100c',
    'gene_alias.txt': '1011067faff7a7ade595b964bba267ad9c8c5671be42159b890681d97b168619',
    'genes.txt': '55c184e89408fd5d4ce5ad334f80bb8f8dd4a8ec35f5c57fb1d972ab91af86e7',
    'pseudobulk.tsv': '41181276d8c6f9f715e855d616600e0dd9fe9c245abd10b56af446d601da7eb1',
    'qc.f32': '17672fcd11e4bc68ace858db263b896da12a5eca4121afc6017dc123aefe1a51',
    'sample.u16': 'f857be3baf7ca82b3d075d3108a08aa7e1b7b3465a719c28cfda5de4c459ba90',
  }
  attempt('every other entry of the part is byte-identical (guard)',
    () => Object.keys(BEFORE).filter(k => !g.files[k] || sha256(g.files[k]) !== BEFORE[k]), [])
  attempt('meta.json is the same in every field but embeddings, notes and run (guard)', () => {
    const m = { ...g.meta }
    delete m.embeddings; delete m.notes; delete m.run
    return sha256(JSON.stringify(m))
  }, '89f9acf2d69cde69c7e6977dc14e2bf45b30050e0b247e6d63e51306aea779cb')
  attempt('and the notes not about embeddings or uns are the same 13, in order (guard)',
    () => sha256(JSON.stringify(g.meta.notes.filter(s => !/embedding/i.test(s) && !/\buns\b/.test(s)))),
    '669678f83b4aefc8cf3835f130fe3e22aee382e190f9d084df9f515d586f892f')
  attempt('collection.json is the same but for embeddings, notes and the part\'s bytes and sha256 (guard)', () => {
    const ix = { ...G.o.index, parts: G.o.index.parts.map(({ bytes: _b, sha256: _s, ...q }) => q) }
    delete ix.embeddings; delete ix.notes
    return sha256(JSON.stringify(ix))
  }, 'a9985a264cf4021d54153cc03fa2bb5ba1d901804355295cfac697d2f3380cc4')
  const again = convert(FIX.classic, 'guessed-again.zip', ['--label', 'family'])
  attempt('two conversions are still byte-identical',
    () => [again.r.code, sha256(readFileSync(again.out)) === sha256(G.o.bytes)], [0, true])
}

rmSync(tmp, { recursive: true, force: true })
console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll integrity tests passed\n')
process.exit(failed ? 1 : 0)
