// Bundle-building regressions.
//
// This is the third independent writer of scrnaseq-studio/bundle@1 (after the
// Python and R exporters), so the format is checked here by reading the zip
// back and comparing it against a dense reference computed the obvious way.
// A transpose that is subtly wrong still produces a bundle that opens.

import { unzipSync, strFromU8 } from 'fflate'
import { buildBundle, chooseMatrices, BuildError } from '../src/lib/build.ts'
import { categorical, numeric, fromFactor } from '../src/lib/scan.ts'
import { cellQC, lognormalize, sumCells, toGeneMajor, dropZeros } from '../src/lib/matrix.ts'
import { CHUNK_GENES, readGenes } from '../src/lib/chunked.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const near = (name, got, want, tol = 1e-5) => {
  const ok = Math.abs(got - want) < tol
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${got}\n        want ${want}`}`)
}
const rejects = (name, fn, fragment) => {
  try {
    fn()
    failed++
    console.log(`  FAIL ${name}\n        it was accepted`)
  } catch (e) {
    const ok = String(e.message).toLowerCase().includes(fragment.toLowerCase())
    if (!ok) failed++
    console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        message was: ${e.message}`}`)
  }
}

/** Cell-major sparse from a dense array of rows. */
function fromDense(rows, nGenes) {
  const indptr = [0], indices = [], data = []
  for (const r of rows) {
    for (let g = 0; g < nGenes; g++) if (r[g] !== 0) { indices.push(g); data.push(r[g]) }
    indptr.push(indices.length)
  }
  return {
    nCells: rows.length, nGenes,
    indptr: Int32Array.from(indptr),
    indices: Int32Array.from(indices),
    data: Float32Array.from(data),
  }
}

// Six cells, four genes, three cell types, two samples in two conditions.
const COUNTS = [
  [5, 0, 1, 0],
  [7, 0, 2, 0],
  [0, 9, 0, 1],
  [0, 8, 0, 0],
  [1, 1, 6, 0],
  [0, 2, 4, 3],
]
const GENES = ['Aaa', 'Bbb', 'Ccc', 'Ddd']
const TYPES = ['T cell', 'B cell', 'Myeloid']

const scanOf = (over = {}) => ({
  format: 'seurat',
  source: 'test.rds',
  nCells: 6,
  columns: [
    categorical('cell_type', ['T cell', 'T cell', 'B cell', 'B cell', 'Myeloid', 'Myeloid']),
    categorical('sample', ['s1', 's1', 's1', 's2', 's2', 's2']),
    categorical('group', ['ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat']),
    numeric('percent.mt', Float32Array.from([1, 2, 3, 4, 5, 6])),
    ...(over.columns ?? []),
  ],
  matrices: over.matrices ?? [
    { key: 'RNA@counts', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'counts',
      load: () => fromDense(COUNTS, 4) },
  ],
  embeddings: [{ key: 'umap', nDims: 2, load: () => Float32Array.from([0,0, 1,1, 2,2, 3,3, 4,4, 5,5]) }],
  geneSets: { RNA: GENES },
  provenance: { normalization: null },
  notes: ['a note from the scan'],
  ...over.scan,
})

const CHOICES = {
  cluster: 'cell_type', sample: 'sample', condition: 'group',
  embedding: 'umap', label: 'Test object',
}

const open = res => {
  const f = unzipSync(res.zip)
  const view = (name, Ctor) => new Ctor(f[name].slice().buffer)
  return {
    meta: JSON.parse(strFromU8(f['meta.json'])),
    genes: strFromU8(f['genes.txt']).split('\n'),
    indptr: view('expr.indptr.i32', Int32Array),
    indices: view('expr.indices.i32', Int32Array),
    data: view('expr.data.f32', Float32Array),
    cluster: view('cluster.u16', Uint16Array),
    sample: view('sample.u16', Uint16Array),
    embed: view('embed.f32', Float32Array),
    qc: view('qc.f32', Float32Array),
    chunkptr: view('expr.chunkptr.i32', Int32Array),
    pseudobulk: f['pseudobulk.tsv'] ? strFromU8(f['pseudobulk.tsv']) : null,
    names: Object.keys(f),
  }
}

console.log('\nTHE TRANSPOSE AGREES WITH A DENSE REFERENCE')
// Gene-major CSC is what every view reads. Getting it wrong shows up as one
// gene's values appearing under another's name, which renders perfectly.
{
  const m = fromDense(COUNTS, 4)
  const g = toGeneMajor(m)
  const dense = []
  for (let gi = 0; gi < 4; gi++) {
    const col = new Array(6).fill(0)
    for (let k = g.indptr[gi]; k < g.indptr[gi + 1]; k++) col[g.indices[k]] = g.data[k]
    dense.push(col)
  }
  check('every value lands under its own gene',
    dense, GENES.map((_x, gi) => COUNTS.map(r => r[gi])))
  check('cell indices come out ascending per gene', (() => {
    for (let gi = 0; gi < 4; gi++) {
      for (let k = g.indptr[gi] + 1; k < g.indptr[gi + 1]; k++) {
        if (g.indices[k] <= g.indices[k - 1]) return false
      }
    }
    return true
  })(), true)
  check('indptr ends at the number of stored values', g.indptr[4], g.indices.length)
  check('nothing is lost', g.data.length, m.data.length)
}

console.log('\nLOG-NORMALIZATION IS log1p(CP10K)')
{
  const m = lognormalize(fromDense([[5, 0, 5, 0]], 4))
  near('a gene at half the library', m.data[0], Math.log1p(5 * 1e4 / 10))
  const zero = lognormalize(fromDense([[0, 0, 0, 0]], 4))
  check('an empty cell does not divide by zero', zero.data.length, 0)
}

console.log('\nQC AND PSEUDOBULK ARE SUMS OF THE RIGHT CELLS')
{
  const m = fromDense(COUNTS, 4)
  const q = cellQC(m)
  check('total counts per cell', Array.from(q.total), COUNTS.map(r => r.reduce((a, b) => a + b, 0)))
  check('detected genes per cell', Array.from(q.genes), COUNTS.map(r => r.filter(v => v !== 0).length))
  check('a pseudobulk column sums only its cells',
    Array.from(sumCells(m, Int32Array.from([0, 1]))),
    [12, 0, 3, 0])
  check('dropZeros removes stored zeros', dropZeros(fromDense([[0, 1]], 2)).data.length, 1)
}

console.log('\nA BUILT BUNDLE IS INTERNALLY CONSISTENT')
{
  const res = buildBundle(scanOf(), CHOICES)
  const b = open(res)
  check('every required entry is present',
    ['meta.json', 'genes.txt', 'cluster.u16', 'sample.u16', 'embed.f32', 'qc.f32',
     'expr.indptr.i32', 'expr.indices.i32', 'expr.data.f32']
      .every(n => b.names.includes(n)), true)
  check('the schema is the one the studio reads', b.meta.schema, 'scrnaseq-studio/bundle@1')
  check('genes', b.genes, GENES)
  check('cell types keep file order, not alphabetical', b.meta.clusters, TYPES)
  check('cluster codes', Array.from(b.cluster), [0, 0, 1, 1, 2, 2])
  check('samples carry their condition', b.meta.samples,
    [{ id: 's1', condition: 'ctrl' }, { id: 's2', condition: 'treat' }])
  check('conditions', b.meta.conditions, ['ctrl', 'treat'])
  check('nnz matches the stored values', b.meta.nnz, b.data.length)
  check('indptr has one entry per gene plus one', b.indptr.length, GENES.length + 1)
  check('the embedding is interleaved x,y', Array.from(b.embed.slice(0, 4)), [0, 0, 1, 1])
  check('QC is interleaved counts,genes,mito',
    Array.from(b.qc.slice(0, 3)), [6, 2, 1])
  check('raw counts were found', b.meta.hasRawCounts, true)
  check('scan notes travel into the bundle',
    b.meta.notes.includes('a note from the scan'), true)
  check('the chosen cluster column is recorded', b.meta.provenance.clustering, 'cell_type')

  // Expression was built from counts, so it must be log1p(CP10K) of them.
  const cell0 = COUNTS[0]
  const tot = cell0.reduce((a, x) => a + x, 0)
  let got = 0
  for (let k = b.indptr[0]; k < b.indptr[1]; k++) if (b.indices[k] === 0) got = b.data[k]
  near('expression is log1p(CP10K) of the counts', got, Math.log1p(cell0[0] * 1e4 / tot), 1e-4)

  const pb = b.pseudobulk.split('\n')
  check('pseudobulk has a column per sample × cluster present',
    pb[0].split('\t').length - 1, res.pseudobulkColumns)
  check('pseudobulk columns are named sample||cluster||n',
    pb[0].split('\t')[1], 's1||T cell||2')
  check('a pseudobulk row is the summed raw counts',
    pb[1].split('\t').slice(0, 3), ['Aaa', '12', '0'])
  check('every gene has a pseudobulk row', pb.length - 1, GENES.length)
}

console.log('\nOMITTING SAMPLE AND CONDITION IS HONEST, NOT INVENTED')
{
  const b = open(buildBundle(scanOf(), { ...CHOICES, sample: null, condition: null }))
  check('one sample', b.meta.samples, [{ id: 'all cells', condition: 'all cells' }])
  check('one condition', b.meta.conditions, ['all cells'])
  check('and it says so', b.meta.notes.some(n => /no condition column/.test(n)), true)
}

console.log('\nMATRIX CHOICE FOLLOWS THE NUMBERS, NOT THE SLOT NAME')
{
  const notes = []
  const scaled = { key: 'X', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'scaled',
    load: () => fromDense(COUNTS.map(r => r.map(v => v - 3)), 4) }
  const logCounts = { key: 'raw.X', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'log-counts',
    load: () => fromDense(COUNTS.map(r => r.map(v => v ? Math.log1p(v) : 0)), 4) }
  const out = chooseMatrices(scanOf({ matrices: [scaled, logCounts] }), s => notes.push(s))
  check('scaled X is passed over for the log-counts in raw', out.counts != null, true)
  check('counts are recovered by exponentiating',
    Array.from(out.counts.data.slice(0, 2)), [5, 1])
  check('and it says where they came from',
    notes.some(n => /exponentiated back to integer counts/.test(n)), true)
}
{
  const lognorm = { key: 'RNA@data', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'lognorm',
    load: () => fromDense(COUNTS.map(r => r.map(v => v ? Math.log1p(v * 100) : 0)), 4) }
  const counts = { key: 'RNA@counts', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'counts',
    load: () => fromDense(COUNTS, 4) }
  const notes = []
  const out = chooseMatrices(scanOf({ matrices: [counts, lognorm] }), s => notes.push(s))
  check('an already-normalized matrix is used as it is, not renormalized',
    notes.some(n => /already log-normalized/.test(n)), true)
  check('and the counts are still kept for pseudobulk', out.counts != null, true)
}
rejects('an object with only scaled data is refused with a fix',
  () => buildBundle(scanOf({
    matrices: [{ key: 'X', geneSet: 'RNA', nGenes: 4, nCells: 6, kind: 'scaled',
      load: () => fromDense(COUNTS, 4) }],
  }), CHOICES),
  'cannot plot z-scores')

console.log('\nMISMATCHES ARE REFUSED RATHER THAN WRITTEN')
rejects('a gene list that does not match the matrix',
  () => buildBundle(scanOf({ scan: { geneSets: { RNA: ['only', 'two'] } } }), CHOICES),
  'gene list has 2')
rejects('an embedding that is not there',
  () => buildBundle(scanOf(), { ...CHOICES, embedding: 'tsne' }), 'no embedding')
rejects('a cluster column that is not a grouping',
  () => buildBundle(scanOf(), { ...CHOICES, cluster: 'percent.mt' }), 'not a usable cell annotation')
check('BuildError is the error type', new BuildError('x') instanceof Error, true)

console.log('\nUNANNOTATED CELLS BECOME THEIR OWN LEVEL')
// Folding NA into level 0 would quietly file every unannotated cell under the
// first cell type, which is the kind of wrongness that renders.
{
  const withNA = fromFactor('anno', Int32Array.from([0, 1, -1, 1, 0, -1]), ['A', 'B'])
  const b = open(buildBundle(scanOf({ columns: [withNA] }), { ...CHOICES, cluster: 'anno' }))
  check('NA is its own cluster', b.meta.clusters, ['A', 'B', 'NA'])
  check('and the codes point at it', Array.from(b.cluster), [0, 1, 2, 1, 0, 2])
}

console.log('\nMITOCHONDRIAL FRACTION IS PUT ON ONE SCALE')
{
  const asProportion = numeric('percent.mt', Float32Array.from([0.01, 0.02, 0.03, 0.04, 0.05, 0.06]))
  const b = open(buildBundle(scanOf({ scan: { columns: [
    categorical('cell_type', ['T cell', 'T cell', 'B cell', 'B cell', 'Myeloid', 'Myeloid']),
    asProportion,
  ] } }), { ...CHOICES, sample: null, condition: null }))
  near('a proportion is converted to a percentage', b.qc[2], 1)
  check('and it says so', b.meta.notes.some(n => /converted to a percentage/.test(n)), true)
}

console.log('\nTHE CHUNKED MATRIX IS STORED, AND SAYS THE SAME THING AS THE FLAT ONE')
// The whole lazy reader rests on two properties of the written zip: that
// expr.chunk.bin is method 0, so chunk k is a contiguous byte range the studio
// can slice out of a Blob; and that what comes back out of it is the same
// matrix as expr.indices/expr.data, which are still written unchanged. A
// deflated entry here would still unzip perfectly and still pass every other
// test in this file — and would make every gene lookup read the whole matrix.
{
  const res = buildBundle(scanOf(), CHOICES)
  const zip = res.zip
  const b = open(res)
  check('both copies of the matrix are in the bundle',
    ['expr.indices.i32', 'expr.data.f32', 'expr.chunk.bin', 'expr.chunkptr.i32']
      .every(n => b.names.includes(n)), true)
  check('meta carries the chunk size a reader must use', b.meta.chunkGenes, CHUNK_GENES)

  // Walk the central directory from the tail — exactly what the studio does.
  const dv = new DataView(zip.buffer, zip.byteOffset, zip.byteLength)
  let eocd = -1
  for (let i = zip.length - 22; i >= 0; i--) if (dv.getUint32(i, true) === 0x06054b50) { eocd = i; break }
  const count = dv.getUint16(eocd + 10, true)
  const found = {}
  let p = dv.getUint32(eocd + 16, true)
  for (let i = 0; i < count; i++) {
    const nameLen = dv.getUint16(p + 28, true)
    const name = strFromU8(zip.subarray(p + 46, p + 46 + nameLen))
    found[name] = { method: dv.getUint16(p + 10, true), size: dv.getUint32(p + 20, true),
      local: dv.getUint32(p + 42, true) }
    p += 46 + nameLen + dv.getUint16(p + 30, true) + dv.getUint16(p + 32, true)
  }
  check('expr.chunk.bin is stored, not deflated', found['expr.chunk.bin'].method, 0)
  check('the flat entries are still compressed', found['expr.data.f32'].method, 8)

  const e = found['expr.chunk.bin']
  const base = e.local + 30 + dv.getUint16(e.local + 26, true) + dv.getUint16(e.local + 28, true)
  check('the stored payload is exactly as long as the writer made it',
    e.size, b.chunkptr[b.chunkptr.length - 1])

  // Read every gene through the byte-range path and compare against the flat
  // entries, value for value.
  const getBytes = async (from, to) => zip.subarray(base + from, base + to)
  const all = GENES.map((_g, i) => i)
  const vecs = await readGenes(getBytes, b.chunkptr, b.indptr, b.meta.chunkGenes, all)
  const flat = all.map(g => ({
    gene: g,
    cells: Array.from(b.indices.subarray(b.indptr[g], b.indptr[g + 1])),
    values: Array.from(b.data.subarray(b.indptr[g], b.indptr[g + 1])),
  }))
  check('every gene reads back identically',
    vecs.map(v => ({ gene: v.gene, cells: Array.from(v.cells), values: Array.from(v.values) })),
    flat)
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll build tests passed\n')
process.exit(failed ? 1 : 0)
