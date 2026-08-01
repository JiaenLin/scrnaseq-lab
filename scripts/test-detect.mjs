// Detection regressions: what the lab infers before it asks.
//
// Everything here decides what the user is shown, so a wrong answer is not a
// crash — it is a default that quietly produces the wrong bundle.

import { classify, stridedSample } from '../src/lib/matrix.ts'
import {
  candidates, categorical, crossCheck, duplicateOf, fromFactor, guessCluster,
  guessEmbedding, guessRole, internGeneSet, maybeDiscrete, numeric,
  replicatesPerCondition, samePartition, unusable,
} from '../src/lib/scan.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const matches = (name, got, re) => {
  const ok = re.test(String(got))
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${got}`}`)
}

console.log('\nA MATRIX IS CLASSIFIED BY ITS NUMBERS')
// The whole point: slot names lie. An .h5ad keeps scaled values in X and
// log1p(counts) in .raw about as often as not.
check('integers are counts', classify([0, 1, 5, 12, 300]), 'counts')
check('any negative means scaled', classify([1.2, -0.5, 3.1]), 'scaled')
check('log1p of integers is log-counts',
  classify([1, 2, 5, 9].map(Math.log1p)), 'log-counts')
check('small positive non-integers are lognorm', classify([0.31, 1.27, 2.9, 4.4]), 'lognorm')
check('large positive non-integers are linear', classify([12.5, 880.2, 3100.7]), 'linear')
check('an all-zero matrix does not crash', classify([0, 0, 0]), 'lognorm')
check('an empty matrix does not crash', classify([]), 'lognorm')
check('stored zeros do not make a lognorm matrix look scaled',
  classify([0, 0, 1.27, 0, 2.9]), 'lognorm')
check('non-finite values are ignored', classify([NaN, 1, 2, Infinity]), 'counts')

console.log('\nSAMPLING SPREADS ACROSS THE WHOLE VECTOR')
// Reading the first 20 000 values of a 200-million-value array and calling it
// representative is a guess; striding costs nothing and is not one.
{
  const n = 1000
  const s = stridedSample(n, 10, i => i)
  check('ten values spread over a thousand', Array.from(s), [0, 100, 200, 300, 400, 500, 600, 700, 800, 900])
  check('asking for more than exists returns what exists',
    stridedSample(3, 10, i => i).length, 3)
  check('the last element is reachable', stridedSample(10, 10, i => i)[9], 9)
  check('a scaled tail is still found',
    classify(stridedSample(1000, 20, i => (i > 900 ? -1 : 1))), 'scaled')
}

console.log('\nNUMERIC COLUMNS THAT ARE REALLY GROUPINGS')
// Cluster ids arrive as 0..8 integers more often than as strings.
{
  const ids = maybeDiscrete(numeric('seurat_clusters', Float32Array.from([0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2])), 12)
  check('small integer columns become categorical', ids.kind, 'categorical')
  check('levels are sorted numerically, not as text', ids.levels, ['0', '1', '2'])
  const mt = maybeDiscrete(numeric('percent.mt', Float32Array.from([1.2, 3.4, 5.6])), 3)
  check('a continuous column stays numeric', mt.kind, 'numeric')
  const counts = maybeDiscrete(numeric('nCount', Float32Array.from(Array.from({ length: 100 }, (_v, i) => i))), 100)
  check('an integer column with one value per cell stays numeric', counts.kind, 'numeric')
}

console.log('\nLEVEL ORDER IS FILE ORDER, NEVER SORTED')
// Sorting would reorder a time course into 0 h, 24 h, 6 h.
check('first appearance wins',
  categorical('time', ['0 h', '0 h', '6 h', '24 h', '6 h']).levels, ['0 h', '6 h', '24 h'])
check('codes point at those levels',
  Array.from(categorical('t', ['b', 'a', 'b']).codes), [0, 1, 0])
check('NA is counted', categorical('x', ['a', 'NA', 'b']).missing, 1)

console.log('\nTHE SAME PARTITION UNDER DIFFERENT NAMES')
// seurat_clusters (0..2) and seurat_annotations (B, NK, T) are routinely the
// same grouping twice. Offering both as "condition" invites a bundle whose
// comparison has nothing to compare.
{
  const anno = categorical('seurat_annotations', ['B', 'NK', 'T', 'B', 'NK', 'T'])
  const nums = categorical('seurat_clusters', ['2', '0', '1', '2', '0', '1'])
  const other = categorical('group', ['ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat'])
  check('a relabelled partition is recognised', samePartition(anno, nums), true)
  check('a genuinely different one is not', samePartition(anno, other), false)
  check('duplicates point back to the first',
    [...duplicateOf([anno, nums, other]).entries()], [['seurat_clusters', 'seurat_annotations']])
  check('different lengths are never the same partition',
    samePartition(anno, categorical('x', ['a', 'b'])), false)
}

console.log('\nCANDIDATES AND EXCLUSIONS')
{
  const scan = {
    nCells: 6,
    columns: [
      numeric('nCount_RNA', Float32Array.from([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])),
      categorical('seurat_clusters', ['0', '1', '2', '0', '1', '2']),
      categorical('cell_type', ['B', 'NK', 'T', 'B', 'NK', 'T']),
      categorical('orig.ident', ['s1', 's1', 's1', 's2', 's2', 's2']),
      categorical('treatment', ['ctrl', 'ctrl', 'ctrl', 'drug', 'drug', 'drug']),
      categorical('barcode', ['c1', 'c2', 'c3', 'c4', 'c5', 'c6']),
    ],
    embeddings: [{ key: 'X_pca', nDims: 50 }, { key: 'X_umap', nDims: 2 }],
  }
  check('a named annotation outranks a numbered one',
    candidates(scan, 'cluster')[0].name, 'cell_type')
  check('a per-cell identifier is not an annotation',
    candidates(scan, 'cluster').some(c => c.name === 'barcode'), false)
  matches('and it says why',
    unusable(scan.columns[5], 'cluster', 6), /nearly one per cell/)
  matches('a measurement is not a grouping',
    unusable(scan.columns[0], 'condition', 6), /not a grouping/)
  check('sample candidates are led by the conventional name',
    candidates(scan, 'sample')[0].name, 'orig.ident')
  check('condition candidates are led by the conventional name',
    candidates(scan, 'condition')[0].name, 'treatment')

  // A column cannot be two roles at once, and neither can its relabelling.
  const afterCluster = candidates(scan, 'condition', ['cell_type'])
  check('the chosen cluster is not offered as a condition',
    afterCluster.some(c => c.name === 'cell_type'), false)
  check('nor is the same partition under another name',
    afterCluster.some(c => c.name === 'seurat_clusters'), false)
  check('unrelated columns are still offered',
    afterCluster.map(c => c.name), ['treatment', 'orig.ident'])

  check('guessing prefers a named annotation', guessCluster(scan), 'cell_type')
  check('sample is guessed', guessRole(scan, 'sample'), 'orig.ident')
  check('condition is guessed', guessRole(scan, 'condition'), 'treatment')
  check('UMAP beats PCA', guessEmbedding(scan), 'X_umap')
  check('PCA is used only when it is all there is',
    guessEmbedding({ ...scan, embeddings: [{ key: 'X_pca', nDims: 50 }] }), 'X_pca')
  check('no embedding at all is reported, not invented',
    guessEmbedding({ ...scan, embeddings: [] }), null)
  check('a one-dimensional reduction is not an embedding',
    guessEmbedding({ ...scan, embeddings: [{ key: 'X_odd', nDims: 1 }] }), null)
}
check('an object with no usable annotation guesses nothing',
  guessCluster({ nCells: 2, columns: [numeric('x', Float32Array.from([1.5, 2.5]))], embeddings: [] }), null)

console.log('\nSAMPLES MUST NEST INSIDE CONDITIONS')
// If they do not, the pairing is two independent labels, and the bundle's
// one-condition-per-sample model would keep whichever came first.
{
  const sample = categorical('sample', ['s1', 's1', 's2', 's2'])
  const good = categorical('group', ['ctrl', 'ctrl', 'drug', 'drug'])
  const bad = categorical('group', ['ctrl', 'drug', 'ctrl', 'drug'])
  check('a nested design passes', crossCheck(sample, good), null)
  matches('a crossed one is caught', crossCheck(sample, bad), /appears in both/)
  check('with nothing to check it passes', crossCheck(undefined, good), null)
}

console.log('\nREPLICATES PER GROUP DECIDE WHETHER PSEUDOBULK IS OFFERED')
{
  const s = categorical('s', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
  const c = categorical('c', ['x', 'x', 'x', 'x', 'y', 'y', 'y', 'y'])
  check('4 v 4', replicatesPerCondition(s, c), { x: 4, y: 4 })
  check('without a sample column every group has one',
    replicatesPerCondition(undefined, c), { x: 1, y: 1 })
  check('without a condition there is nothing to count', replicatesPerCondition(s, undefined), {})
}

console.log('\nGENE SETS ARE INTERNED BY CONTENT')
// Keying by slot name made @counts and @data look like different gene sets,
// and pseudobulk was disabled to avoid a mismatch that did not exist.
{
  const sets = {}
  const a = internGeneSet(sets, ['G1', 'G2'], 'counts')
  const b = internGeneSet(sets, ['G1', 'G2'], 'data')
  const c = internGeneSet(sets, ['G1', 'G3'], 'other')
  check('identical lists share one key', a === b, true)
  check('a different list gets its own', a === c, false)
  check('only the distinct lists are stored', Object.keys(sets).length, 2)
  check('same length, different names, is different',
    internGeneSet(sets, ['G9', 'G8'], 'x') === a, false)
}

console.log('\nFACTORS WITH MISSING LEVELS')
{
  const f = fromFactor('anno', Int32Array.from([0, -1, 1]), ['A', 'B'])
  check('NA codes are counted as missing', f.missing, 1)
  check('levels are kept as given', f.levels, ['A', 'B'])
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll detection tests passed\n')
process.exit(failed ? 1 : 0)
