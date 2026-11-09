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
}
const BIN = here('../bin/scrnaseq-lab.mjs')
const tmp = mkdtempSync(join(tmpdir(), 'lab-integrity-test-'))
const sha256 = (b) => createHash('sha256').update(b).digest('hex')

const build = await tryImport('../src/lib/build.ts')
const scanMod = await tryImport('../src/lib/scan.ts')
const shard = await tryImport('../src/lib/shard.ts')
const { categorical, fromFactor } = scanMod

/** Run the real command. */
function cli(args, env = {}) {
  const r = spawnSync(process.execPath, [BIN, ...args], {
    encoding: 'utf8', env: { ...process.env, ...env }, maxBuffer: 64 << 20,
  })
  let summary = null
  try { summary = JSON.parse(r.stdout) } catch { /* not a summary */ }
  return { code: r.status, out: r.stdout ?? '', err: r.stderr ?? '', summary }
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
function convert(input, name, args) {
  const out = join(tmp, name)
  const r = cli(['convert', input, out, ...args])
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
  attempt('two embeddings equal in one part\'s cells but not the whole object\'s are both kept, in every part (guard)',
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
    const ix = { ...G.o.index, parts: G.o.index.parts.map(({ bytes, sha256: _s, ...q }) => q) }
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
