// Bundle-building regressions.
//
// This is the third independent writer of scrnaseq-studio/bundle@1 (after the
// Python and R exporters), so the format is checked here by reading the zip
// back and comparing it against a dense reference computed the obvious way.
// A transpose that is subtly wrong still produces a bundle that opens.

import { unzipSync, strFromU8 } from 'fflate'
import { buildBundle, chooseMatrices, streamableMatrix, BuildError } from '../src/lib/build.ts'
import { categorical, numeric, fromFactor, geneNameKind, pickGeneAlias } from '../src/lib/scan.ts'
import { cellQC, lognormalize, sumCells, toGeneMajor, dropZeros } from '../src/lib/matrix.ts'
import { shardScan } from '../src/lib/shard.ts'
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


console.log(String.fromCharCode(10) + 'THE SPLIT PATH PICKS THE SAME MATRIX AS THE WHOLE PATH')
// It used to take the first streamable matrix, which on the commonest atlas
// layout - a scaled X beside raw counts in .raw - is the scaled one. That
// either ships z-scores as expression or dies at the end of a twenty-minute
// conversion saying every matrix is scaled, on a file that plainly is not.
{
  const mk = (key, kind, streamable) => ({
    key, geneSet: 'RNA', nGenes: 4, nCells: 6, kind,
    load: () => fromDense(COUNTS, 4),
    ...(streamable ? { loadGroups: () => [fromDense(COUNTS, 4)] } : {}),
  })
  const pick = ms => { const m = streamableMatrix(scanOf({ matrices: ms })); return m && m.key }
  check('a scaled X is never chosen over raw counts',
    pick([mk('X', 'scaled', true), mk('raw.X', 'counts', true)]), 'raw.X')
  check('an already log-normalized matrix wins, as it does unsplit',
    pick([mk('X', 'counts', true), mk('RNA@data', 'lognorm', true)]), 'RNA@data')
  check('a matrix that cannot be streamed is not chosen',
    pick([mk('X', 'lognorm', false), mk('raw.X', 'counts', true)]), 'raw.X')
  check('nothing streamable but scaled data returns null, so the caller can say why',
    pick([mk('X', 'scaled', true)]), null)
}


console.log(String.fromCharCode(10) + 'A MATRIX OVER A CORNER OF THE GENE LIST IS NOT THE EXPRESSION')
// Measured on the real atlas: layers["expected"] is a velocity model's output.
// Its values are small positive non-integers, which classifies as lognorm and
// therefore wins the preference order — and it holds values for exactly 2000 of
// 31053 genes, so the studio would report 94% of the gene list as expressed
// nowhere. Value classification cannot see this; coverage can.
{
  const mk = (key, kind, coverage) => ({
    key, geneSet: 'RNA', nGenes: 4, nCells: 6, kind,
    load: () => fromDense(COUNTS, 4),
    loadGroups: () => [fromDense(COUNTS, 4)],
    ...(coverage === undefined ? {} : { coverage: () => coverage }),
  })
  const notes = []
  const pick = ms => {
    notes.length = 0
    const m = streamableMatrix(scanOf({ matrices: ms }), s => notes.push(s))
    return m && m.key
  }
  check('a 6%-coverage lognorm loses to counts that cover the list',
    pick([mk('layers["expected"]', 'lognorm', 0.064), mk('X', 'counts', 0.9)]), 'X')
  check('and it says why, in the notes that travel into the bundle',
    notes.some(n => n.includes('layers["expected"]') && n.includes('6% of the genes')), true)
  check('a lognorm with real coverage is still preferred',
    pick([mk('RNA@data', 'lognorm', 0.8), mk('X', 'counts', 0.9)]), 'RNA@data')
  check('a sparse matrix with nothing to lose against is kept',
    pick([mk('X', 'lognorm', 0.05)]), 'X')
  check('two equally thin matrices are both kept — this chooses, it never refuses',
    pick([mk('X', 'lognorm', 0.05), mk('raw.X', 'counts', 0.06)]), 'X')
  check('a source that cannot report coverage is never vetoed',
    pick([mk('layers["expected"]', 'lognorm', undefined), mk('X', 'counts', 0.9)]),
    'layers["expected"]')
  // The unsplit path has to make the same decision, or an object converts one
  // way whole and another way split.
  const res = buildBundle(scanOf({ matrices: [
    mk('layers["expected"]', 'lognorm', 0.064), mk('X', 'counts', 0.9),
  ] }), CHOICES)
  check('the whole-object path agrees with the split path',
    res.meta.notes.some(n => n.includes('expression computed here as log1p(CP10K) from the counts')), true)
}

console.log(String.fromCharCode(10) + 'EVERY 2D EMBEDDING IS CARRIED, NOT JUST THE CHOSEN ONE')
// An object routinely holds a UMAP and a t-SNE of the same cells. The lab used
// to keep one and discard the rest, which is a decision the user never made and
// could not undo without converting again.
{
  const umap = Float32Array.from([0,0, 1,1, 2,2, 3,3, 4,4, 5,5])
  const tsne = Float32Array.from([9,9, 8,8, 7,7, 6,6, 5,5, 4,4])
  const pca = Float32Array.from([1,2, 3,4, 5,6, 7,8, 9,10, 11,12])
  const res = buildBundle(scanOf({ scan: { embeddings: [
    { key: 'X_umap', nDims: 2, load: () => umap },
    { key: 'X_tSNE', nDims: 2, load: () => tsne },
    { key: 'X_pca', nDims: 50, load: () => pca },
  ] } }), { ...CHOICES, embedding: 'X_umap' })
  const b = open(res)
  check('the chosen one is still embed.f32, byte for byte', Array.from(b.embed), Array.from(umap))
  check('meta names it as the default', b.meta.embedding, 'X_umap')
  check('and lists them all, the default first', b.meta.embeddings, [
    { key: 'X_umap', file: 'embed.f32' },
    { key: 'X_tSNE', file: 'embed.X_tSNE.f32' },
    { key: 'X_pca', file: 'embed.X_pca.f32' },
  ])
  const f = unzipSync(res.zip)
  const read = name => Array.from(new Float32Array(f[name].slice().buffer))
  check('the t-SNE is really in the file', read('embed.X_tSNE.f32'), Array.from(tsne))
  check('and so is the PCA', read('embed.X_pca.f32'), Array.from(pca))

  // Choosing a different default must move the bytes, not just the label: a
  // reader that knows nothing of meta.embeddings reads embed.f32 alone.
  const b2 = open(buildBundle(scanOf({ scan: { embeddings: [
    { key: 'X_umap', nDims: 2, load: () => umap },
    { key: 'X_tSNE', nDims: 2, load: () => tsne },
  ] } }), { ...CHOICES, embedding: 'X_tSNE' }))
  check('picking the t-SNE puts the t-SNE in embed.f32', Array.from(b2.embed), Array.from(tsne))
  check('and the UMAP moves to a named entry', b2.meta.embeddings, [
    { key: 'X_tSNE', file: 'embed.f32' },
    { key: 'X_umap', file: 'embed.X_umap.f32' },
  ])
}

console.log(String.fromCharCode(10) + 'ONE EMBEDDING STILL WRITES EXACTLY WHAT IT USED TO')
{
  const b = open(buildBundle(scanOf(), CHOICES))
  check('no extra embedding entries', b.names.filter(n => n.startsWith('embed.')), ['embed.f32'])
  check('and the list has the one', b.meta.embeddings, [{ key: 'umap', file: 'embed.f32' }])
}

console.log(String.fromCharCode(10) + 'AN UNREADABLE EXTRA EMBEDDING DOES NOT LOSE THE CONVERSION')
// Twenty minutes of reading must not be thrown away over a decoration.
{
  const b = open(buildBundle(scanOf({ scan: { embeddings: [
    { key: 'umap', nDims: 2, load: () => Float32Array.from([0,0, 1,1, 2,2, 3,3, 4,4, 5,5]) },
    { key: 'broken', nDims: 2, load: () => { throw new Error('bad chunk') } },
    { key: 'short', nDims: 2, load: () => Float32Array.from([0, 0]) },
  ] } }), CHOICES))
  check('the bundle is written anyway', b.meta.nCells, 6)
  check('only the readable one is listed', b.meta.embeddings, [{ key: 'umap', file: 'embed.f32' }])
  check('and it says why', b.meta.notes.some(n => n.includes('broken') && n.includes('bad chunk')), true)
  check('a wrong-length one is refused too',
    b.meta.notes.some(n => n.includes('short') && n.includes('1 points for 6 cells')), true)
}

console.log(String.fromCharCode(10) + 'GENE NAMES ARE CLASSIFIED BY WHAT THEY LOOK LIKE')
{
  check('Ensembl mouse', geneNameKind(['ENSMUSG00000038751', 'ENSMUSG00000030762']), 'accession')
  check('Ensembl human', geneNameKind(['ENSG00000141510', 'ENSG00000012048']), 'accession')
  check('versioned Ensembl', geneNameKind(['ENSG00000141510.17', 'ENSG00000012048.23']), 'accession')
  check('symbols', geneNameKind(['Sox2', 'Actb', 'Gapdh', 'Mki67']), 'symbol')
  // A symbol that starts with an accession's letters is still a symbol; the
  // digits are what make an identifier.
  check('ENSA is a real gene symbol', geneNameKind(['ENSA', 'Sox2', 'Actb']), 'symbol')
  check('bare Entrez ids', geneNameKind(['100038431', '17869', '11461']), 'accession')
  check('two references concatenated',
    geneNameKind(['ENSG00000141510', 'ENSG00000012048', 'Sox2', 'Actb']), 'mixed')
}

console.log(String.fromCharCode(10) + 'THE SYMBOL COLUMN IS FOUND BY CONTENT, NOT BY ITS NAME')
{
  const acc = ['ENSMUSG00000000001', 'ENSMUSG00000000002', 'ENSMUSG00000000003', 'ENSMUSG00000000004']
  const pick = cols => pickGeneAlias(acc, cols)
  check('a column called Gene holding symbols',
    pick([{ name: 'Gene', values: ['Sox2', 'Actb', 'Gapdh', 'Mki67'] }]),
    { kind: 'symbol', column: 'Gene', names: ['Sox2', 'Actb', 'Gapdh', 'Mki67'] })
  check('a chromosome column is not a naming — too few distinct values',
    pick([{ name: 'Chromosome', values: ['1', '1', '2', '2'] }]), null)
  check('nor is a biotype', pick([{ name: 'feature_type',
    values: ['Gene Expression', 'Gene Expression', 'Gene Expression', 'Gene Expression'] }]), null)
  check('a column called Gene holding accessions is not the symbols',
    pick([{ name: 'Gene', values: acc }]), null)
  check('the conventional name wins when two columns both qualify', pick([
    { name: 'other_name', values: ['aa', 'bb', 'cc', 'dd'] },
    { name: 'gene_symbol', values: ['Sox2', 'Actb', 'Gapdh', 'Mki67'] },
  ]).column, 'gene_symbol')
  // The reverse layout: rows are symbols already and the accessions are the
  // column. Scanpy writes exactly this.
  check('symbols indexed, accessions in var',
    pickGeneAlias(['Sox2', 'Actb', 'Gapdh', 'Mki67'], [{ name: 'gene_ids', values: acc }]),
    { kind: 'accession', column: 'gene_ids', names: acc })
  // A mixed index — two references concatenated — takes symbols if it can get
  // them, because that is what makes gene sets match.
  check('a mixed index prefers the symbols', pickGeneAlias(
    ['ENSG00000141510', 'ENSG00000012048', 'Sox2', 'Actb'], [
      { name: 'ensembl_id', values: ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000181449', 'ENSG00000075624'] },
      { name: 'gene_symbol', values: ['TP53', 'BRCA1', 'SOX2', 'ACTB'] },
    ]).kind, 'symbol')
  check('and falls back to the accessions when there are no symbols', pickGeneAlias(
    ['ENSG00000141510', 'ENSG00000012048', 'Sox2', 'Actb'],
    [{ name: 'ensembl_id', values: ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000181449', 'ENSG00000075624'] }]
  ).kind, 'accession')
}

console.log(String.fromCharCode(10) + 'BOTH NAMINGS TRAVEL IN THE BUNDLE')
{
  const ACC = ['ENSMUSG00000000001', 'ENSMUSG00000000002', 'ENSMUSG00000000003', 'ENSMUSG00000000004']
  const res = buildBundle(scanOf({ scan: {
    geneSets: { RNA: ACC },
    geneAliases: { RNA: { kind: 'symbol', column: 'Gene', names: ['Sox2', 'Actb', 'Sox2', ''] } },
  } }), CHOICES)
  const b = open(res)
  const f = unzipSync(res.zip)
  check('genes.txt is still what the matrix rows are indexed by', b.genes, ACC)
  check('meta says which naming that is', b.meta.geneIdKind, 'accession')
  // Row 4 had no symbol, so it keeps its accession — every line is usable.
  check('the alias is a file beside it', strFromU8(f['gene_alias.txt']).split('\n'),
    ['Sox2', 'Actb', 'Sox2', 'ENSMUSG00000000004'])
  check('and meta describes it', b.meta.geneAlias,
    { kind: 'symbol', column: 'Gene', file: 'gene_alias.txt', missing: 1, duplicated: 2 })
  check('rows sharing a symbol are not merged', b.meta.nGenes, 4)
  check('the empty one is reported',
    b.meta.notes.some(n => n.includes('1 of 4 genes have no symbol')), true)
  check('and so is the collision',
    b.meta.notes.some(n => n.includes('2 genes share a symbol')), true)
}

console.log(String.fromCharCode(10) + 'WITHOUT AN ALIAS THE BUNDLE IS WHAT IT ALWAYS WAS')
{
  const b = open(buildBundle(scanOf(), CHOICES))
  check('no alias file', b.names.includes('gene_alias.txt'), false)
  check('meta says there is none', b.meta.geneAlias, null)
  check('and still classifies what it does have', b.meta.geneIdKind, 'symbol')
}

console.log(String.fromCharCode(10) + 'AN ACCESSION-ONLY OBJECT IS TOLD IT WILL NOT MATCH GENE SETS')
// The reason the studio's built-in sets matched nothing on the real atlas.
{
  const ACC = ['ENSMUSG00000000001', 'ENSMUSG00000000002', 'ENSMUSG00000000003', 'ENSMUSG00000000004']
  const b = open(buildBundle(scanOf({ scan: { geneSets: { RNA: ACC } } }), CHOICES))
  check('it says the sets will not match',
    b.meta.notes.some(n => n.includes('written as symbols will not match')), true)
}

console.log(String.fromCharCode(10) + 'AN EXTRA COLUMN TRAVELS WITH THE CELLS')
// Whatever else the object knows about a cell. Not a fourth role: it is carried
// under its own name, and nothing downstream is told which one it is.
const REGION = categorical('dissection',
  ['Forebrain', 'Forebrain', 'Midbrain', 'Midbrain', 'Hindbrain', 'Hindbrain'])
/** One extra column's codes, read back out of the zip. */
const extraOf = (res, file) =>
  [...new Uint16Array(unzipSync(res.zip)[file].slice().buffer)]
{
  const res = buildBundle(scanOf({ columns: [REGION] }), { ...CHOICES, extras: ['dissection'] })
  const b = open(res)
  check('meta names it, under the object\'s own name',
    b.meta.extras.map(e => e.key), ['dissection'])
  check('with its levels in the object\'s order',
    b.meta.extras[0].levels, ['Forebrain', 'Midbrain', 'Hindbrain'])
  check('and an entry named by meta, one level index per cell',
    extraOf(res, b.meta.extras[0].file), [0, 0, 1, 1, 2, 2])
  check('the three roles are untouched',
    [b.meta.clusters.length, b.meta.samples.length, b.meta.conditions.length], [3, 2, 2])
  check('and which column each role came from is on the record',
    [b.meta.provenance.clustering, b.meta.provenance.condition, b.meta.provenance.sample],
    ['cell_type', 'group', 'sample'])
}

console.log(String.fromCharCode(10) + 'A PART DROPS THE LEVELS IT HAS NO CELLS FOR')
// Already fixed once for cluster and sample. A shard that inherited the
// parent's level list would make its own "region 2" a different region from
// every other shard's, and every figure would still draw.
{
  const parent = scanOf({ columns: [REGION] })
  const cells = Int32Array.from([2, 3, 4, 5])
  const sub = shardScan(parent, { key: 'later half', cells, nnz: 0 },
    new Map([['RNA@counts', fromDense(COUNTS.slice(2), 4)]]), 'test.rds · shard')
  const res = buildBundle(sub, { ...CHOICES, extras: ['dissection'] })
  const b = open(res)
  check('the level this part does not hold is gone',
    b.meta.extras[0].levels, ['Midbrain', 'Hindbrain'])
  check('and the codes were renumbered onto what is left',
    extraOf(res, b.meta.extras[0].file), [0, 0, 1, 1])
  check('exactly as the cluster column is', b.meta.clusters, ['B cell', 'Myeloid'])
  check('the note quotes the whole object, so every part says the same',
    b.meta.notes.some(n => n.includes('dissection travels with the cells') && n.includes('3 levels')),
    true)
}

console.log(String.fromCharCode(10) + 'WITHOUT ONE, THE BUNDLE IS WHAT IT ALWAYS WAS')
{
  const b = open(buildBundle(scanOf(), CHOICES))
  check('no extra entry', b.names.some(n => n.startsWith('extra.')), false)
  check('and nothing claims there is one', b.meta.extras, [])
}

console.log(String.fromCharCode(10) + 'A COLUMN ALREADY SPOKEN FOR IS NOT CARRIED TWICE')
{
  const b = open(buildBundle(scanOf(), { ...CHOICES, extras: ['group', 'percent.mt'] }))
  check('the condition column is not duplicated', b.meta.extras, [])
  check('and a measurement is refused with a reason',
    b.meta.notes.some(n => n.includes('percent.mt') && n.includes('not a categorical')), true)
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll build tests passed\n')
process.exit(failed ? 1 : 0)
