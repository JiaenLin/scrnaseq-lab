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
import { copyFileSync, mkdirSync, mkdtempSync, readFileSync, rmSync } from 'node:fs'
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
  stored: here('./fixtures/stored-outside.h5ad'),
  storedBin: here('./fixtures/stored-outside.bin'),
  damaged: here('./fixtures/damaged-reads.h5ad'),
  half: here('./fixtures/damaged-half.h5ad'),
  opaque: here('./fixtures/damaged-opaque.h5ad'),
  onlyMatrix: here('./fixtures/damaged-only-matrix.h5ad'),
  indptr: here('./fixtures/damaged-indptr.h5ad'),
  lzfCodes: here('./fixtures/lzf-codes.h5ad'),
  lzfGenes: here('./fixtures/lzf-genes.h5ad'),
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

/** Every file of a written bundle, unzipped, as one text: collection.json and every entry of every part. */
const bundleText = (o) => [strFromU8(unzipSync(o.bytes)['collection.json']),
  ...o.parts.flatMap(p => Object.values(p.files).map(b => Buffer.from(b).toString('latin1')))].join('\n')
const SIX = ['tool', 'version', 'commit', 'stage', 'run_key', 'job']
const recOf = (part, p) => part?.meta.run.records.find(x => x.path === p)?.record
const notesOf = (part, re) => part.meta.notes.filter(s => re.test(s))
/** A new HDF5 file written through h5wasm. */
const write = (path, edit) => { const f = new h5wasm.File(path, 'w'); try { edit(f) } finally { f.flush(); f.close() } }
/** Edit a file in place through h5wasm. */
const edit = (path, fn) => { const f = new h5wasm.File(path, 'a'); try { fn(f) } finally { f.flush(); f.close() } }

console.log('\nTHE RUN\'S IDENTITY TRAVELS WITH THE BUNDLE, NOT ITS RECORD')
{
  check('the guessed conversion exits 0', G.r.code, 0)
  attempt('meta.run lists the one record the object has, and where it was',
    () => g.meta.run.records.map(r => r.path), ['uns/scintegrate/provenance'])
  attempt('and carries its identity as recorded, under the six names, and nothing else',
    () => JSON.stringify(g.meta.run.records[0].record),
    JSON.stringify({ tool: 'scintegrate', version: '0.0.0+fixture', commit: '0'.repeat(40), job: 'fixture' }))
  attempt('meta.run says where it looked', () => typeof g.meta.run.searched === 'string' && /provenance/.test(g.meta.run.searched) && /harness_run/.test(g.meta.run.searched), true)
  attempt('a note says what of the record is carried, and counts the rest without naming it',
    () => notesOf(g, /^uns\/scintegrate\/provenance:/), ['uns/scintegrate/provenance: run carries its tool, version, commit and job; not carried: 9 other entries'])
  attempt("the record's argv, input path and settings are in no file of the bundle",
    () => ['joint.h5ad', 'wrapped_versions', 'label_key', 'argv'].filter(s => bundleText(G.o).includes(s)), [])

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
    const prov = dict(dict(uns, 'analysis'), 'provenance')
    str(prov, 'tool', 'synthetic-tool')
    int(prov, 'n', 3)
    const ct = dict(uns, 'cell_type_provenance')
    str(ct, 'tree', 'synthetic')
    int(ct, 'n_cells', 240)
    str(dict(dict(dict(uns, 'deep'), 'inner'), 'provenance'), 'tool', 'two levels down: not a record by the rule')
    str(uns, 'harness_run', JSON.stringify({ stages: [{ stage: 'integration', tool: 'synthetic', run_key: 'HARNESS-RK-TOKEN', job: 'HARNESS-JOB-TOKEN' }] }))
    const qc = dict(uns, 'qc_decl')
    str(qc, 'schema', 'synthetic/provenance@1')
    str(qc, 'run_key', 'rk-0')
    str(qc, 'commit', 'abc1234')
    str(dict(uns, 'unrelated'), 'schema', 'synthetic/results@1')
  })
  const r = convert(many, 'records.zip', ['--label', 'records'])
  const part = r.o?.parts[0]
  attempt('each record the rule finds that holds one of the six is listed, in uns order, and nothing else is',
    () => part.meta.run.records.map(x => x.path), ['uns/analysis/provenance', 'uns/qc_decl', 'uns/scintegrate/provenance'])
  attempt('a record known by its schema carries its identity fields',
    () => JSON.stringify(recOf(part, 'uns/qc_decl')), JSON.stringify({ commit: 'abc1234', run_key: 'rk-0' }))
  attempt("a record with none of the six, an annotator's column record, is not listed, and its note says so",
    () => notesOf(part, /^uns\/cell_type_provenance holds none of tool, version, commit, stage, run_key and job as a string of 256 characters or fewer, so run does not list it; not carried: 2 other entries$/).length, 1)
  attempt("the harness's record, one JSON string, is named as a string and nothing of it is read",
    () => [notesOf(part, /^uns\/harness_run is a string, not a dict of named fields, so run reads nothing from it$/).length,
      ['HARNESS-RK-TOKEN', 'HARNESS-JOB-TOKEN'].filter(s => bundleText(r.o).includes(s))], [1, []])
}

console.log('\nNO STRING OF A RECORD BUT ITS SIX FIELDS REACHES THE BUNDLE')
{
  // A record shaped as scIntegrate's writer builds it (cli.py: the job as a
  // dict of scheduler, id and host, argv, the input path, nested settings),
  // with a host path and an argument list planted in it.
  const host = copyWith(FIX.classic, 'planted.h5ad', f => {
    const p = dict(f.get('uns'), 'site_provenance')
    str(p, 'tool', 'planted-tool')
    str(p, 'version', '1.2.3')
    str(p, 'commit', 'f'.repeat(40))
    enc(p.create_dataset({ name: 'argv', data: ['tool', 'run', '--h5ad', '/planted/host/dir/input.h5ad', '--token', 'PLANTED-ARGV-TOKEN'] }), 'string-array', '0.2.0')
    str(p, 'input', '/planted/host/dir/input.h5ad')
    const job = dict(p, 'job')
    str(job, 'scheduler', 'pbs')
    str(job, 'id', '4242.PLANTED-HOST')
    str(job, 'host', 'PLANTED-HOST')
    str(dict(dict(p, 'method_settings'), 'harmony'), 'note', 'PLANTED-DEEP-TOKEN')
  })
  const h = convert(host, 'planted.zip', ['--label', 'planted'])
  const all = h.o ? bundleText(h.o) : ''
  attempt('it converts, and carries the tool, version and commit',
    () => [h.r.code, JSON.stringify(recOf(h.o?.parts[0], 'uns/site_provenance'))], [0, JSON.stringify({ tool: 'planted-tool', version: '1.2.3', commit: 'f'.repeat(40) })])
  attempt('no host path, argument, job host or nested value of the record is in any file of the bundle',
    () => ['/planted/host', 'PLANTED-ARGV-TOKEN', 'PLANTED-HOST', 'PLANTED-DEEP-TOKEN', 'method_settings'].filter(s => all.includes(s)), [])
  attempt('a job written as a dict, as scIntegrate writes it, is named as a group and not carried',
    () => notesOf(h.o.parts[0], /^uns\/site_provenance: run carries its tool, version and commit; not carried: job \(a group\) and 3 other entries$/).length, 1)
}

console.log('\nONE FILE: A LINK OUT OF THE INPUT IS NAMED ONCE AND NEVER READ')
{
  // The last re-checker's first reproduction: a record that is an external
  // link by absolute path to a file outside the input's folder, and one to a
  // file beside the input. Both were read as the object's own, named nowhere,
  // and the bundle changed when only the other file did. Links out of the
  // matrix's layers, the embeddings and one of a record's six fields, and a
  // soft link to nothing, are each named too.
  const dir = mkdtempSync(join(tmpdir(), 'lab-onefile-'))
  const away = mkdtempSync(join(tmpdir(), 'lab-away-'))
  const secret = join(away, 'secret.h5')
  const sibling = join(dir, 'other.h5')
  write(secret, f => str(dict(f, 'secret'), 'token', 'OUTSIDE-INPUT-TOKEN'))
  write(sibling, f => { const r = dict(f, 'rec'); str(r, 'tool', 'SIBLING-TOOL'); str(r, 'origin', 'SIBLING-FILE-CONTENT') })
  const input = join(dir, 'ext.h5ad')
  copyFileSync(FIX.classic, input)
  edit(input, f => {
    const uns = f.get('uns')
    uns.create_external_link(secret, '/secret', 'abs_provenance')
    uns.create_external_link('other.h5', '/rec', 'e_provenance')
    const p = dict(uns, 'f_provenance')
    str(p, 'tool', 'F-TOOL')
    p.create_external_link('other.h5', '/rec/tool', 'version')
    uns.create_soft_link('/nowhere', 'o_provenance')
    uns.create_soft_link('/uns/scintegrate/provenance', 'mirror_provenance')
    f.get('obsm').create_external_link('other.h5', '/rec', 'X_away')
    f.get('layers').create_external_link('other.h5', '/rec', 'away')
  })
  const one = convert(input, 'ext1.zip', ['--label', 'ext'])
  const part = one.o?.parts[0]
  const all = one.o ? bundleText(one.o) : ''
  attempt('it converts, and nothing of either other file is in the bundle',
    () => [one.r.code, ['OUTSIDE-INPUT-TOKEN', 'SIBLING-FILE-CONTENT', 'SIBLING-TOOL'].filter(s => all.includes(s))], [0, []])
  attempt('one note names each link not followed, once, by its own place in the file, and carries no target',
    () => {
      const n = notesOf(part, /^the converter reads only the file it is given/)
      const out = 'an external link, which points outside this file'
      return [n.length, ...[`uns/abs_provenance (${out})`, `uns/e_provenance (${out})`, `uns/f_provenance/version (${out})`,
        'uns/o_provenance (a soft link that reaches nothing in this file)', `obsm/X_away (${out})`, `layers/away (${out})`]
        .map(s => (n[0] ?? '').split(s).length - 1), [secret, away, 'other.h5', '/nowhere'].filter(t => all.includes(t))]
    }, [1, 1, 1, 1, 1, 1, 1, []])
  attempt('an ordinary soft link inside the file still reads: the record it reaches carries its identity',
    () => [JSON.stringify(recOf(part, 'uns/mirror_provenance')), JSON.stringify(recOf(part, 'uns/f_provenance'))],
    [JSON.stringify(recOf(g, 'uns/scintegrate/provenance')), JSON.stringify({ tool: 'F-TOOL' })])
  // The same input bytes, whatever the other files hold, and when they are gone.
  write(secret, f => str(dict(f, 'secret'), 'token', 'CHANGED-TOKEN'))
  write(sibling, f => str(dict(f, 'rec'), 'tool', 'OTHER-TOOL-2'))
  const two = convert(input, 'ext2.zip', ['--label', 'ext'])
  rmSync(secret)
  rmSync(sibling)
  const three = convert(input, 'ext3.zip', ['--label', 'ext'])
  attempt('the same input bytes give the same bundle bytes when the linked files change, and when they are gone',
    () => [two.r.code, three.r.code, sha256(two.o.bytes) === sha256(one.o.bytes), sha256(three.o.bytes) === sha256(one.o.bytes)], [0, 0, true, true])
  rmSync(dir, { recursive: true, force: true })
  rmSync(away, { recursive: true, force: true })
}

console.log('\nA DATASET WHOSE VALUES ARE STORED IN ANOTHER FILE IS NAMED, NOT READ')
{
  // scripts/fixtures/stored-outside.h5ad, written by h5py: an embedding whose
  // raw values HDF5 keeps in stored-outside.bin (external storage), which
  // h5wasm cannot write. Under the seal HDF5 could not open that file, and
  // h5wasm handed back whatever its buffer held, with no error: offered as an
  // embedding, it was memory, not the object.
  const dir = mkdtempSync(join(tmpdir(), 'lab-stored-'))
  const input = join(dir, 'stored-outside.h5ad')
  copyFileSync(FIX.stored, input)
  const a = convert(input, 'stored1.zip', ['--label', 'stored'])
  copyFileSync(FIX.storedBin, join(dir, 'stored-outside.bin'))
  const b = convert(input, 'stored2.zip', ['--label', 'stored'])
  attempt('it converts, the embedding is not offered, and one note names it',
    () => [a.r.code, JSON.stringify(a.o.parts[0].meta.embeddings).includes('X_stored_outside'),
      notesOf(a.o.parts[0], /obsm\/X_stored_outside \(a dataset whose values are stored in another file\)/).length], [0, false, 1])
  attempt('and the bundle is the same whether that file sits beside the input or not',
    () => [b.r.code, sha256(a.o.bytes) === sha256(b.o.bytes)], [0, true])
  rmSync(dir, { recursive: true, force: true })
}

console.log('\nTWO FILES OF ONE NAME CANNOT COLLIDE: NEITHER IS READ')
{
  // The last re-checker's third reproduction: x.h5 beside the input and x.h5
  // in sub/, reached by two links from one record. An identity made of the
  // link's file name and an address made the second "the first again": its
  // content was lost, and the bundle said what was false.
  const dir = mkdtempSync(join(tmpdir(), 'lab-collide-'))
  mkdirSync(join(dir, 'sub'))
  for (const [where, marker] of [[dir, 'X-TOP'], [join(dir, 'sub'), 'X-SUB']]) {
    write(join(where, 'x.h5'), f => { const r = dict(f, 'r'); str(r, 'tool', marker); str(r, 'marker', marker) })
  }
  write(join(dir, 'sub', 'y.h5'), f => dict(f, 'g').create_external_link('x.h5', '/r', 'link'))
  const input = join(dir, 'collide.h5ad')
  copyFileSync(FIX.classic, input)
  edit(input, f => {
    const p = dict(f.get('uns'), 'col_provenance')
    str(p, 'tool', 'COL-TOOL')
    p.create_external_link('x.h5', '/r', 'a_top')
    p.create_external_link('sub/y.h5', '/g', 'b_sub')
  })
  const c = convert(input, 'collide.zip', ['--label', 'collide'])
  const part = c.o?.parts[0]
  attempt('the record carries its own identity, and nothing of either file of the one name is in the bundle',
    () => [c.r.code, JSON.stringify(recOf(part, 'uns/col_provenance')), ['X-TOP', 'X-SUB'].filter(s => bundleText(c.o).includes(s))],
    [0, JSON.stringify({ tool: 'COL-TOOL' }), []])
  attempt('and no sentence says that one link reaches what another did',
    () => /again/.test(JSON.stringify(part.meta.run) + part.meta.notes.join('\n')), false)
  rmSync(dir, { recursive: true, force: true })
}

console.log('\nNO NAME FROM THE OBJECT BECOMES A KEY, SO A FIELD NAMED LIKE AN OBJECT PROPERTY CANNOT VANISH')
{
  // The last re-checker's fourth reproduction: a record's field named
  // __proto__, a string in one record and a group in another, vanished from
  // the copy, and nothing said so. No field name comes from the object now:
  // each such entry is counted among those not carried.
  const proto = copyWith(FIX.classic, 'proto.h5ad', f => {
    const uns = f.get('uns')
    const a = dict(uns, 'proto_provenance')
    str(a, 'tool', 'PROTO-TOOL')
    str(a, '__proto__', 'PROTO-VALUE')
    str(a, 'constructor', 'CTOR-VALUE')
    const b = dict(uns, 'proto2_provenance')
    str(b, 'tool', 'PROTO2-TOOL')
    str(dict(b, '__proto__'), 'hidden', 'HIDDEN-VALUE')
  })
  const pr = convert(proto, 'proto.zip', ['--label', 'proto'])
  const part = pr.o?.parts[0]
  attempt("each record carries its tool, under the converter's own key, and nothing else",
    () => [pr.r.code, ...['uns/proto_provenance', 'uns/proto2_provenance'].map(p => JSON.stringify(recOf(part, p)))],
    [0, JSON.stringify({ tool: 'PROTO-TOOL' }), JSON.stringify({ tool: 'PROTO2-TOOL' })])
  attempt('and its note counts the entries named __proto__ and constructor among those not carried',
    () => [notesOf(part, /^uns\/proto_provenance: run carries its tool; not carried: 2 other entries$/).length,
      notesOf(part, /^uns\/proto2_provenance: run carries its tool; not carried: 1 other entry$/).length], [1, 1])
  attempt('every key of every record is one of the six', () => part.meta.run.records.every(x => Object.keys(x.record).every(k => SIX.includes(k))), true)
}

console.log('\nAT MOST 64 RECORDS ARE READ; THE REST ARE COUNTED, AND THE FIRST TEN NAMED')
{
  const crowd = copyWith(FIX.classic, 'crowd.h5ad', f => {
    const uns = f.get('uns')
    for (let i = 0; i < 100; i++) str(dict(uns, `r${String(i).padStart(3, '0')}_provenance`), 'tool', `tool-${i}`)
  })
  const c = convert(crowd, 'crowd.zip', ['--label', 'crowd'])
  const run = c.o?.parts[0].meta.run
  attempt('the first 64 in uns order are read, and the 37 after them counted, the first ten named',
    () => [c.r.code, run.records.length, run.records[63].path, run.leftOut.count, run.leftOut.paths],
    [0, 64, 'uns/r063_provenance', 37, Array.from({ length: 10 }, (_v, i) => `uns/r${String(64 + i).padStart(3, '0')}_provenance`)])
  attempt('and one note says so',
    () => notesOf(c.o.parts[0], /^101 records were found and run reads the first 64: the 37 after them were not read: uns\/r064_provenance, /).length, 1)
}

console.log('\nA FIELD OF THE WRONG KIND IS NAMED WITH ITS KIND, NOT CARRIED')
{
  const kinds = copyWith(FIX.classic, 'kinds.h5ad', f => {
    const uns = f.get('uns')
    const a = dict(uns, 'a_provenance')
    str(a, 'tool', 'kind-tool')
    int(a, 'version', 2)
    enc(a.create_dataset({ name: 'commit', data: Float64Array.from([1.5]), shape: [] }), 'numeric-scalar', '0.2.0')
    enc(a.create_dataset({ name: 'stage', data: ['s1', 's2', 's3'] }), 'string-array', '0.2.0')
    // A real h5py boolean (an HDF5 enum), linked from the fixture's own uns.
    a.create_hard_link('/uns/scintegrate/assessment/summary/needed', 'run_key')
    str(dict(a, 'job'), 'id', '1')
    const b = dict(uns, 'b_provenance')
    str(b, 'tool', 'b-tool')
    str(b, 'commit', 'c'.repeat(256))
    str(b, 'version', '\u{1F600}'.repeat(256))
    const c = dict(uns, 'c_provenance')
    str(c, 'tool', 'c-tool')
    str(c, 'commit', '\u{1F600}'.repeat(257))
    str(dict(uns, 'd_provenance'), 'tool', 'd-tool')
  })
  const k = convert(kinds, 'kinds.zip', ['--label', 'kinds'])
  const part = k.o?.parts[0]
  attempt('the string is carried; a number, an array, a boolean and a group are each named with their kind',
    () => [k.r.code, JSON.stringify(recOf(part, 'uns/a_provenance')), notesOf(part, /^uns\/a_provenance:/)],
    [0, JSON.stringify({ tool: 'kind-tool' }),
      ['uns/a_provenance: run carries its tool; not carried: version (a number), commit (a number), stage (an array of 3 values), run_key (a boolean) and job (a group)']])
  attempt('a string of 256 characters is carried whole, counted as characters, not UTF-16 units',
    () => { const b = recOf(part, 'uns/b_provenance'); return [b?.commit?.length, [...(b?.version ?? '')].length] }, [256, 256])
  attempt('the first string longer than 256 characters is named, ends the reading, and the records after it are counted, unread',
    () => [JSON.stringify(recOf(part, 'uns/c_provenance')), part.meta.run.leftOut,
      notesOf(part, /^uns\/c_provenance\/commit is a string longer than 256 characters, so run read no value after it: 2 records found after it were not read: uns\/d_provenance, uns\/scintegrate\/provenance$/).length],
    [JSON.stringify({ tool: 'c-tool' }), { count: 2, paths: ['uns/d_provenance', 'uns/scintegrate/provenance'] }, 1])
}

console.log("\nTHE READER ON THE SHAPES THE FAMILY'S TOOLS WRITE")
{
  // Small objects shaped as each writer's code builds its record (the root
  // design review lists them): an integration stage's uns/<tool>/provenance,
  // the QC stage's declaration (found by its schema), the annotator's column
  // record, the profiler's uns entry (no record by the rule), the harness's
  // uns/harness_run (one JSON string).
  const fam = copyWith(FIX.classic, 'family-shapes.h5ad', f => {
    const uns = f.get('uns')
    const si = dict(dict(uns, 'integrate'), 'provenance')
    str(si, 'tool', 'scintegrate')
    str(si, 'version', '0.4.1')
    str(si, 'commit', 'b'.repeat(40))
    int(si, 'state_version', 3)
    enc(si.create_dataset({ name: 'argv', data: ['scintegrate', 'integrate', '--h5ad', 'joint.h5ad'] }), 'string-array', '0.2.0')
    str(si, 'python', '3.12.1')
    const job = dict(si, 'job')
    str(job, 'scheduler', 'pbs')
    str(job, 'id', '77.queue')
    str(job, 'host', 'node-1')
    str(dict(si, 'wrapped_versions'), 'harmony', '0.0.10')
    enc(dict(si, 'sees').create_dataset({ name: 'harmony', data: ['sample'] }), 'string-array', '0.2.0')
    str(si, 'input', 'joint.h5ad')
    int(si, 'seed', 0)
    int(dict(dict(si, 'method_settings'), 'harmony'), 'theta', 2)
    const qc = dict(uns, 'scqc')
    for (const [k2, v] of [['schema', 'scqc/provenance@1'], ['tool', 'scQC'], ['version', '0.9.0'], ['commit', 'a'.repeat(40)],
      ['run_key', 'RK-QC-1'], ['sample', 'S1'], ['written', '2026-01-01T00:00:00+00:00'], ['flag_column', 'qc_flag'],
      ['flag_meaning', 'flagged by step 6'], ['flag_digest', 'd'.repeat(64)], ['flag_na_as', 'False'], ['watch_column', ''],
      ['watch_meaning', ''], ['note', 'scQC removed nothing on the flag']]) str(qc, k2, v)
    for (const [k2, v] of [['n_obs', 240], ['n_flagged', 3], ['n_examined', 240], ['removed_on_flag', 0]]) int(qc, k2, v)
    const an = dict(uns, 'cell_type_provenance')
    str(an, 'column', 'cell_type')
    str(an, 'source', 'independent')
    str(an, 'tree', 'tree-1')
    int(an, 'n_cells', 240)
    int(dict(an, 'labels'), 'B', 120)
    str(an, 'replaced', '')
    int(an, 'n_disagree', -1)
    dict(an, 'disagreements')
    const pf = dict(uns, 'scprofile')
    str(pf, 'contract', '1.0')
    str(dict(pf, 'input'), 'path', 'object.h5ad')
    str(dict(dict(pf, 'kernels'), 'de'), 'method', 'pseudobulk')
    str(uns, 'harness_run', JSON.stringify({ stages: [{ stage: 'qc', tool: 'scQC', run_key: 'RK-QC-1' }] }))
  })
  const fr = convert(fam, 'family-shapes.zip', ['--label', 'shapes'])
  const part = fr.o?.parts[0]
  attempt('the records found are those that hold one of the six, in uns order, each with its identity',
    () => [fr.r.code, JSON.stringify(part.meta.run.records)],
    [0, JSON.stringify([
      { path: 'uns/integrate/provenance', record: { tool: 'scintegrate', version: '0.4.1', commit: 'b'.repeat(40) } },
      { path: 'uns/scintegrate/provenance', record: { tool: 'scintegrate', version: '0.0.0+fixture', commit: '0'.repeat(40), job: 'fixture' } },
      { path: 'uns/scqc', record: { tool: 'scQC', version: '0.9.0', commit: 'a'.repeat(40), run_key: 'RK-QC-1' } }])])
  attempt("each record's note: the integration stage's job is a group; the declaration's other 14 entries, the column record's 8, are counted",
    () => [notesOf(part, /^uns\/integrate\/provenance: run carries its tool, version and commit; not carried: job \(a group\) and 8 other entries$/).length,
      notesOf(part, /^uns\/scqc: run carries its tool, version, commit and run_key; not carried: 14 other entries$/).length,
      notesOf(part, /^uns\/cell_type_provenance holds none of .* so run does not list it; not carried: 8 other entries$/).length,
      notesOf(part, /^uns\/harness_run is a string/).length, notesOf(part, /scprofile/).length], [1, 1, 1, 1, 0])
  // The deployed Studio says of every leftOut entry beside a listed record that
  // it was left out over the bundle's budget: no budget left these two out.
  attempt('beside the listed records, the column record and the harness string are in the notes only: nothing is counted in leftOut',
    () => [fr.r.code, part.meta.run.leftOut ?? null], [0, null])
}

console.log('\nRUN.LEFTOUT BY THE STATE OF THE RECORDS: BESIDE A LISTED RECORD, ONLY WHAT THE BUDGET LEFT OUT')
{
  // Fix round 2's truth table: {0, 1, 2 listed} x {no other, one holding none
  // of the six, one not a dict} x {no tail, past the cap, after a long value}.
  // uns iterates by name, so names set the order: the others (a..), the listed
  // (l..), the tail (z..). Past the cap needs 64 records read before it, so
  // there the other kind fills what the listed do not (64 - L of it); with no
  // other kind the 64 are listed, the one row of that column. The long value is
  // the last listed record's version; with none listed it is a tool, and that
  // record holds none of the six. The deployed Studio says, with no record
  // listed, "None of the N run records the converter found in the object was
  // copied into this bundle: ..."; beside one, "N more records the converter
  // found, left out of this bundle over its budget: ...".
  const p2 = (i) => String(i).padStart(2, '0')
  const rows = []
  for (const L of [0, 1, 2]) {
    for (const O of ['none', 'six-less', 'not-a-dict']) {
      for (const T of ['none', 'cap', 'long']) {
        if (O === 'none' && T === 'cap') { if (L === 0) rows.push({ L: 64, O, T }); continue }
        rows.push({ L, O, T })
      }
    }
  }
  for (const { L, O, T } of rows) {
    const others = O === 'none' ? 0 : T === 'cap' ? 64 - L : 1
    const name = `tt-${L}-${O}-${T}`
    const c = convert(copyWith(FIX.mask, `${name}.h5ad`, f => {
      const uns = f.get('uns')
      for (let i = 0; i < others; i++) {
        if (O === 'six-less') int(dict(uns, `a${p2(i)}_provenance`), 'tool', 7)
        else str(uns, `a${p2(i)}_provenance`, '{"tool": "h"}')
      }
      for (let i = 0; i < L; i++) {
        const r = dict(uns, `l${p2(i)}_provenance`)
        str(r, 'tool', `l${p2(i)}-tool`)
        if (T === 'long' && i === L - 1) str(r, 'version', 'L'.repeat(257))
      }
      if (T === 'long' && L === 0) str(dict(uns, 'l00_provenance'), 'tool', 'L'.repeat(257))
      for (let i = 0; i < (T === 'cap' ? 2 : T === 'long' ? 1 : 0); i++) str(dict(uns, `z${i}_provenance`), 'tool', `z${i}-tool`)
    }), `${name}.zip`, ['--label', 'p'])
    const run = c.o?.parts[0].meta.run
    const a = Array.from({ length: others }, (_v, i) => `uns/a${p2(i)}_provenance`)
    const z = Array.from({ length: T === 'cap' ? 2 : T === 'long' ? 1 : 0 }, (_v, i) => `uns/z${i}_provenance`)
    const left = L === 0 ? [...a, ...(T === 'long' ? ['uns/l00_provenance'] : []), ...z] : z
    attempt(`${L} listed, ${O === 'none' ? 'no other' : `${others} ${O}`}, ${{ none: 'no tail', cap: '2 past the cap', long: '1 after a long value' }[T]}: leftOut ${left.length ? `counts ${left.length}` : 'is absent'}${L && others ? ', the others in the notes only' : ''}`,
      () => [c.r.code, run?.records.map(x => x.path), run?.leftOut ?? null,
        notesOf(c.o.parts[0], /^uns\/a\d\d_provenance (holds none of|is a string, not a dict)/).length],
      [0, Array.from({ length: L }, (_v, i) => `uns/l${p2(i)}_provenance`),
        left.length ? { count: left.length, paths: left.slice(0, 10) } : null, others])
  }
}

console.log('\nA RECORD OF ODD SHAPES CONVERTS, AND ONLY ITS IDENTITY IS READ')
{
  // scripts/fixtures/uns-records.h5ad, written by anndata and h5py: a
  // structured array with an int64 field (which once crashed the whole
  // conversion), float16, an HDF5 reference, uint64 past 2**53, NaN and a
  // dataframe. None is one of the six, so none is read; each is counted.
  const u = convert(FIX.records, 'uns-records.zip', ['--label', 'records'])
  check('an object whose provenance holds an int64 rec-array converts', u.r.code, 0)
  attempt('its odd record carries its tool, and its other entries are counted, none read',
    () => [JSON.stringify(recOf(u.o?.parts[0], 'uns/odd_provenance')), notesOf(u.o.parts[0], /^uns\/odd_provenance: run carries its tool; not carried: \d+ other entries$/).length],
    [JSON.stringify({ tool: 'synthetic-tool' }), 1])
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
  attempt('and the record beside it still carries its identity',
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
  // Flagged by its width, never by its name (W064): a key holding "pca" is not a record of a PCA.
  attempt('choosing X_pca draws its first two columns, flagged by its width and not by its name',
    () => [pm.embedding, pm.notes.some(s => /principal components/.test(s)),
      pm.notes.some(s => /^X_pca has 10 columns and only its first two are drawn/.test(s))], ['X_pca', false, true])
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
console.log('\nA READ HDF5 REPORTS AS FAILED IS NEVER USED AS VALUES')
{
  // scripts/fixtures/damaged-*.h5ad, written by h5py (scripts/make-damaged-fixtures.py): 8 bytes
  // overwritten inside one gzip chunk of obsm/X_umap, X/data and obs/sample/codes; three
  // embeddings whose stored type is not numbers; an embedding of 16-bit floats. h5wasm drops
  // the status of HDF5's read, so on 8b4997e each came back as plausible numbers with no note
  // (or, for the 16-bit one, stopped the conversion). Every value carried is compared with the
  // source's: the same object converted from its undamaged original.
  const embOf = (part) => Object.fromEntries((part?.meta.embeddings ?? []).map(e => [e.key, sha256(part.files[e.file])]))
  /** An embedding's values as the source file stores them, read straight from it. */
  const source = (file, key) => {
    const f = new h5wasm.File(file, 'r')
    try { const v = f.get(`obsm/${key}`).value; return sha256(new Uint8Array(v.buffer, v.byteOffset, v.byteLength)) } finally { f.close() }
  }
  const rest = (part) => Object.keys(part?.files ?? {}).filter(k => k !== 'meta.json' && !k.startsWith('embed.')).sort()
    .map(k => `${k}:${sha256(part.files[k])}`)
  const guessed = convert(FIX.damaged, 'dmg-guess.zip', ['--label', 'p'])
  attempt('codes HDF5 cannot read stop the conversion when their column would be the guessed sample, naming the dataset',
    () => [guessed.r.code, /obs\/sample could not be read \(HDF5 could not read obs\/sample\/codes /.test(guessed.r.err)], [3, true])
  const asked = convert(FIX.damaged, 'dmg-asked.zip', ['--label', 'p', '--sample', 'sample'])
  attempt('and when a flag names that column',
    () => [asked.r.code, /--sample sample: obs\/sample could not be read \(HDF5 could not read obs\/sample\/codes /.test(asked.r.err)], [3, true])
  const d = convert(FIX.damaged, 'dmg.zip', ['--label', 'p', '--sample', 'none'])
  const c = convert(FIX.classic, 'dmg-clean.zip', ['--label', 'p', '--sample', 'none'])
  const dp = d.o?.parts[0], cp = c.o?.parts[0]
  attempt('without that column it converts, the column named as unreadable in a note',
    () => [d.r.code, notesOf(dp, /^obs\/sample could not be read \(HDF5 could not read obs\/sample\/codes /).length], [0, 1])
  attempt('an embedding with a damaged chunk is named and never offered; every embedding carried is the source\'s, value for value',
    () => { const de = embOf(dp); return [Object.keys(de).sort(), Object.keys(de).every(k => de[k] === source(FIX.classic, k)),
      notesOf(dp, /^obsm\/X_umap could not be read \(HDF5 could not read obsm\/X_umap /).length] },
    [['X_umap_harmony', 'X_umap_none'], true, 1])
  attempt('an opaque, a compound and an empty embedding are each named by what they are stored as, and never offered',
    () => [Object.keys(embOf(dp)).filter(k => /X_opaque|X_compound|X_null/.test(k)), notesOf(dp, /^obsm\/X_(opaque|compound|null) /)],
    [[], ['obsm/X_compound is stored as a compound type, not numbers; it is not offered as an embedding',
      'obsm/X_null is a dataset with no values (an empty dataspace); it is not offered as an embedding',
      'obsm/X_opaque is stored as opaque bytes, not numbers; it is not offered as an embedding']])
  attempt('a matrix with a damaged chunk is named and not used, the note says which was used, and the expression and every other file are the source\'s',
    () => [notesOf(dp, /^X could not be read \(HDF5 could not read X\/data /).length, notesOf(dp, /^expression taken from /),
      JSON.stringify(rest(dp)) === JSON.stringify(rest(cp)), rest(dp).length],
    [1, ['expression taken from layers["lognorm"], already log-normalized'], true, rest(cp).length])
  const only = convert(FIX.onlyMatrix, 'dmg-only.zip', ['--label', 'p'])
  attempt('when the matrix HDF5 cannot read is the only one, the conversion stops and names it',
    () => [only.r.code, /no matrix in this object can be read, so the conversion stops: X \(HDF5 could not read X\/data /.test(only.r.err)], [3, true])
  // The split reads the first matrix's row lengths whole, after the scan: on
  // f2ee369 a failed read there was an internal error with a trace (exit 1).
  const ip = convert(FIX.indptr, 'dmg-indptr.zip', ['--label', 'p'])
  attempt('a row-pointer array HDF5 cannot read stops the conversion as the input\'s error, naming it: exit 3, no trace',
    () => [ip.r.code, /^scrnaseq-lab: X cannot be read, so the conversion stops: HDF5 could not read X\/indptr /m.test(ip.r.err), /internal error|\n\s+at /.test(ip.r.err)], [3, true, false])
  const h = convert(FIX.half, 'dmg-half.zip', ['--label', 'p'])
  const m = convert(FIX.mask, 'dmg-mask.zip', ['--label', 'p'])
  attempt('an embedding of 16-bit floats, the one a guess takes first, is named and left out; the conversion goes on with the next, the source\'s values',
    () => (h.o ? [h.r.code, notesOf(h.o.parts[0], /^obsm\/X_umap is stored as 16-bit floats, which this reader does not read; it is not offered/).length,
      Object.keys(embOf(h.o.parts[0])), embOf(h.o.parts[0]).X_tsne === source(FIX.half, 'X_tsne'), embOf(m.o?.parts[0]).X_umap === source(FIX.half, 'X_tsne')] : [h.r.code]),
    [0, 1, ['X_tsne'], true, true])
  const op = convert(FIX.opaque, 'dmg-opaque.zip', ['--label', 'p'])
  attempt('an embedding of opaque bytes, the one a guess takes first, is named and never offered: its bytes are never read as coordinates',
    () => (op.o ? [op.r.code, notesOf(op.o.parts[0], /^obsm\/X_umap is stored as opaque bytes, not numbers; it is not offered/).length,
      Object.keys(embOf(op.o.parts[0])), embOf(op.o.parts[0]).X_tsne === source(FIX.opaque, 'X_tsne')] : [op.r.code]),
    [0, 1, ['X_tsne'], true])

  // The browser's worker reads through h5wasm's other build, with the File mounted as the worker mounts it.
  const esm = await tryImport('../node_modules/h5wasm/dist/esm/hdf5_hl.js')
  const h5ad = await tryImport('../src/lib/h5ad.ts')
  let w = null
  try {
    globalThis.FileReaderSync ??= class { readAsArrayBuffer(b) { return b.buf.buffer.slice(b.buf.byteOffset, b.buf.byteOffset + b.buf.byteLength) } }
    const bytes = readFileSync(FIX.damaged)
    const r = await h5ad.scanH5adFile(esm.default ?? esm, { name: 'damaged-reads.h5ad', size: bytes.length, slice: (a, b) => ({ buf: bytes.subarray(a, b) }) })
    w = r.scan
    // Loaded while the file is open: on 8b4997e an embedding was read only when loaded.
    w.loaded = Object.fromEntries(w.embeddings.map(e => { const v = e.load(); return [e.key, sha256(new Uint8Array(v.buffer, v.byteOffset, v.byteLength))] }))
    r.handle.close()
  } catch (e) { w = { error: e.message } }
  attempt('in the browser build too: the same parts are named, none is offered, and an embedding carried is the source\'s values',
    () => [w.embeddings.map(e => e.key), w.matrices.map(x => x.key), (w.unreadableColumns ?? []).map(u => u.name),
      w.notes.filter(n => /^(X|obsm\/X_umap|obs\/sample) could not be read \(HDF5 could not read /.test(n)).length,
      ['X_umap', 'X_umap_harmony', 'X_umap_none'].filter(k => k in w.loaded).every(k => w.loaded[k] === source(FIX.classic, k))], [['X_harmony', 'X_none', 'X_pca', 'X_umap_harmony', 'X_umap_none'], ['layers["counts"]', 'layers["lognorm"]'], ['sample'], 3, true])
}

console.log('\nA RECORD FOUND AND NOT LISTED IS COUNTED WHERE A READER OF THE BUNDLE LOOKS')
{
  const runOf = (x) => x.o?.parts[0].meta.run
  const bare = convert(FIX.mask, 'run-bare.zip', ['--label', 'p'])
  const none = convert(copyWith(FIX.mask, 'none-of-six.h5ad', f => int(dict(f.get('uns'), 'x_provenance'), 'tool', 7)), 'none-of-six.zip', ['--label', 'p'])
  const odd = convert(copyWith(FIX.mask, 'not-a-dict.h5ad', f => {
    str(dict(f.get('uns'), 'a_provenance'), 'tool', 'a-tool')
    str(f.get('uns'), 'harness_run', '{"tool": "h"}')
  }), 'not-a-dict.zip', ['--label', 'p'])
  attempt('no record found: none listed, nothing counted',
    () => [bare.r.code, runOf(bare)?.records, runOf(bare)?.leftOut ?? null], [0, [], null])
  attempt('the only record found holds none of the six: none listed, and it is counted and named in leftOut',
    () => [none.r.code, runOf(none)?.records, runOf(none)?.leftOut], [0, [], { count: 1, paths: ['uns/x_provenance'] }])
  attempt('a record that is not a dict, beside one that is listed, is named in the notes and not counted in leftOut',
    () => [odd.r.code, runOf(odd)?.records, runOf(odd)?.leftOut ?? null, notesOf(odd.o?.parts[0], /^uns\/harness_run is a string, not a dict/).length],
    [0, [{ path: 'uns/a_provenance', record: { tool: 'a-tool' } }], null, 1])
  const cut = convert(copyWith(FIX.mask, 'cut-record.h5ad', f => {
    const r = dict(f.get('uns'), 'a_provenance')
    str(r, 'tool', 'CUT-TOOL'); str(r, 'version', 'L'.repeat(257)); str(r, 'commit', 'CUT-COMMIT'); str(r, 'job', 'CUT-JOB')
  }), 'cut-record.zip', ['--label', 'p'])
  attempt('a record cut by a long value keeps the fields before it, and its note names the fields it holds that were not read',
    () => [cut.r.code, runOf(cut)?.records, notesOf(cut.o?.parts[0], /^uns\/a_provenance:/), ['CUT-COMMIT', 'CUT-JOB'].filter(t => bundleText(cut.o).includes(t))],
    [0, [{ path: 'uns/a_provenance', record: { tool: 'CUT-TOOL' } }],
      ['uns/a_provenance: run carries its tool; not carried: version (a string longer than 256 characters), commit (not read: it comes after that string) and job (not read: it comes after that string)'], []])
}

console.log('\nA DATASET COMPRESSED WITH A FILTER THIS READER LACKS IS SAID AS THAT, WITH THE WAY FORWARD')
{
  // scripts/fixtures/lzf-*.h5ad, LZF-compressed by h5py: HDF5 looks for the
  // filter in its own plugin directory, and the seal refuses that open. On
  // f2ee369 each was said to be values stored in another file, false, with no
  // way forward to keep the part.
  const FWD = 'write the file again with gzip compression'
  const g = convert(FIX.lzfGenes, 'lzf-genes.zip', ['--label', 'p'])
  const c = convert(FIX.lzfCodes, 'lzf-codes.zip', ['--label', 'p'])
  const n = convert(FIX.lzfCodes, 'lzf-codes-none.zip', ['--label', 'p', '--sample', 'none'])
  attempt('gene names compressed with LZF stop the conversion, said as a filter this reader does not have, with the way forward',
    () => [g.r.code, g.r.err.includes(`var/_index is a dataset compressed with a filter this reader does not have, LZF for one: ${FWD}`)], [3, true])
  attempt('sample codes compressed with LZF stop a guessed conversion, said as that, never as stored in another file',
    () => [c.r.code, c.r.err.includes(`obs/sample/codes is a dataset compressed with a filter this reader does not have, LZF for one: ${FWD}`), /stored in another file/.test(c.r.err)], [3, true, false])
  attempt("with --sample none it converts, and the bundle's note says the same",
    () => [n.r.code, notesOf(n.o?.parts[0], new RegExp(`^obs/sample/codes is a dataset compressed with a filter this reader does not have, LZF for one, so it is not read: ${FWD}$`)).length,
      notesOf(n.o?.parts[0], /stored in another file/).length], [0, 1, 0])
}

console.log('\nA NAME FROM THE INPUT IS CUT IN A NOTE OR A LIST, SO META.JSON STAYS SMALL')
{
  // 70 records whose names are 20,013 characters each: on 8b4997e, meta.json was 3.2 MB.
  const big = convert(copyWith(FIX.mask, 'long-names.h5ad', f => {
    for (let i = 0; i < 70; i++) str(dict(f.get('uns'), `${String(i).padStart(2, '0')}${'n'.repeat(20000)}_provenance`), 'tool', `t${i}`)
  }), 'long-names.zip', ['--label', 'p'])
  const part = big.o?.parts[0]
  const cutName = `uns/00${'n'.repeat(198)}… (cut: a name of 20,013 characters)`
  attempt('meta.json and collection.json each stay under 100,000 bytes',
    () => [big.r.code, part?.files['meta.json'].length < 100_000, strFromU8(unzipSync(big.o.bytes)['collection.json']).length < 100_000], [0, true, true])
  attempt('each name is cut at 200 characters, and the cut is marked: in the records, in leftOut (past the 64th) and in the notes',
    () => [part.meta.run.records.length, part.meta.run.records[0].path, part.meta.run.leftOut.count, part.meta.run.leftOut.paths[0],
      notesOf(part, /^uns\/00n+… \(cut: a name of 20,013 characters\): run carries its tool$/).length],
    [64, cutName, 6, cutName.replace('uns/00', 'uns/64'), 1])
}

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
  attempt('meta.json is the same in every field but embeddings, notes, run and facts (guard)', () => {
    const m = { ...g.meta }
    delete m.embeddings; delete m.notes; delete m.run; delete m.facts
    return sha256(JSON.stringify(m))
  }, '89f9acf2d69cde69c7e6977dc14e2bf45b30050e0b247e6d63e51306aea779cb')
  // W064 restated one of the 13, the sentinel note ("They stay drawn" was not
  // measured); its coordinates are now said in a note of their own, measured.
  attempt('and the notes not about embeddings or uns are the same 13, in order, the sentinel note restated (guard)',
    () => sha256(JSON.stringify(g.meta.notes.filter(s => !/embedding/i.test(s) && !/\buns\b/.test(s)))),
    'a1070dbd1905511b3d75cfd3322cb8e5a47262a24d8bdd5b89ae03f705be4470')
  attempt('collection.json is the same but for embeddings, notes and the part\'s bytes and sha256 (guard)', () => {
    const ix = { ...G.o.index, parts: G.o.index.parts.map(({ bytes: _b, sha256: _s, ...q }) => q) }
    delete ix.embeddings; delete ix.notes
    return sha256(JSON.stringify(ix))
  }, 'a9985a264cf4021d54153cc03fa2bb5ba1d901804355295cfac697d2f3380cc4')
  const again = convert(FIX.classic, 'guessed-again.zip', ['--label', 'family'])
  // W064: a55a9f06... before; only meta.json (facts added, the sentinel note
  // restated, one measured note) and collection.json's notes moved.
  attempt('the bundle is the one predicted: the closed reader, then the facts', () => sha256(G.o.bytes),
    'd596404e79864ac31351d498001b48b34b620d950684d69f8e8818dd56b44d06')
  attempt('two conversions are still byte-identical',
    () => [again.r.code, sha256(readFileSync(again.out)) === sha256(G.o.bytes)], [0, true])
}

rmSync(tmp, { recursive: true, force: true })
console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll integrity tests passed\n')
process.exit(failed ? 1 : 0)
