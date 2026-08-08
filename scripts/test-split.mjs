// The automatic split decision.
//
// The user is never asked how to split, so this function is the only thing
// standing between them and a file that will not open. It has to be right about
// the small case (leave it alone) as firmly as about the large one.

import { categorical, numeric } from '../src/lib/scan.ts'
import { COMFORTABLE_NNZ, chooseSplit, planSplit, splitCandidates } from '../src/lib/shard.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}

/** n cells, `per` nonzeros each, labelled by a repeating set of level names. */
const object = (levels, per) => {
  const donor = categorical('donor_id', levels)
  const nnz = Int32Array.from(levels.map(() => per))
  return { columns: [donor], nnz, nCells: levels.length }
}

const rep = (n, f) => Array.from({ length: n }, (_x, i) => f(i))

console.log('\nA SMALL OBJECT IS LEFT WHOLE')
// Splitting something that already fits only makes it harder to look at, and
// the studio would then have to stitch four files back together for nothing.
{
  const o = object(rep(1000, i => `d${i % 4}`), 1000)  // 1 M nonzeros
  check('one part', chooseSplit(o.columns, o.nnz, o.nCells), null)
  check('and exactly at the limit, still one part',
    chooseSplit(o.columns, Int32Array.from(o.nnz).fill(COMFORTABLE_NNZ / 1000), o.nCells), null)
}

console.log('\nA LARGE OBJECT IS SPLIT WITHOUT BEING ASKED')
{
  // 1000 cells at 500 k nonzeros each: 500 M in all, 31 M per donor once cut.
  const o = object(rep(1000, i => `d${i % 16}`), COMFORTABLE_NNZ / 100)
  const d = chooseSplit(o.columns, o.nnz, o.nCells)
  check('it decided to split', d !== null, true)
  check('along the only column there is', d.column.name, 'donor_id')
  check('into one part per level', d.plan.shards.length, 16)
  check('every part fits', d.plan.oversized.length, 0)
  check('and it says why, with the number and the column in it',
    /500 M stored values.*16 parts.*donor_id/.test(d.reason), true)
  console.log(`        "${d.reason}"`)
}

console.log('\nA SOURCE THAT CANNOT BE COSTED IS NEVER SPLIT')
// .rds has no row offsets; it is decompressed whole anyway, so there is nothing
// to measure and nothing to gain. Guessing here would split objects that fit.
{
  const o = object(rep(1000, i => `d${i % 8}`), 1)
  check('no row offsets, no split', chooseSplit(o.columns, null, o.nCells), null)
}

console.log('\nA COLUMN THAT DOES NOT DIVIDE ANYTHING IS NOT USED')
{
  // One level for every cell, and one level for the whole object: neither is a
  // split. splitCandidates already drops both, so chooseSplit has nothing left.
  const cols = [
    categorical('barcode', rep(200, i => `cell${i}`)),
    categorical('assay', rep(200, () => '10x')),
    numeric('n_counts', Float32Array.from(rep(200, () => 5))),
  ]
  const nnz = Int32Array.from(rep(200, () => COMFORTABLE_NNZ / 100))
  check('nothing splits it, so it stays whole', chooseSplit(cols, nnz, 200), null)
}

console.log('\nTHE BEST CUT WINS, NOT THE FIRST ONE')
{
  // `tissue` makes two pieces, one of which is still far too big; `donor` makes
  // ten that all fit. The ranking has to prefer donor even though tissue makes
  // fewer files.
  const n = 1000
  const tissue = categorical('tissue', rep(n, i => (i < 900 ? 'brain' : 'cord')))
  const donor = categorical('donor', rep(n, i => `d${i % 10}`))
  const nnz = Int32Array.from(rep(n, () => COMFORTABLE_NNZ / 100))
  const d = chooseSplit([tissue, donor], nnz, n)
  check('it picked the column that leaves nothing oversized', d.column.name, 'donor')
  check('tissue really would have left a part too big',
    planSplit(tissue, nnz, n).oversized.length, 1)
  check('and the ranking agrees',
    splitCandidates([tissue, donor], nnz, n).map(c => c.column.name), ['donor', 'tissue'])
}

console.log('\nAN IMPOSSIBLE OBJECT IS STILL CONVERTED, NOT REFUSED')
{
  // Two levels, both enormous. There is no cut that fits, but the studio reads
  // genes lazily now, so a large part is slow to write rather than unopenable —
  // and half a loaf beats handing the file back.
  const n = 100
  const col = categorical('donor', rep(n, i => (i < 50 ? 'a' : 'b')))
  const nnz = Int32Array.from(rep(n, () => COMFORTABLE_NNZ / 10))  // 500 M total
  const d = chooseSplit([col], nnz, n)
  check('it still split', d.plan.shards.length, 2)
  check('and does not pretend the parts are small', d.plan.oversized.length, 2)
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll split tests passed\n')
process.exit(failed ? 1 : 0)
