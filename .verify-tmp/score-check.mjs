import fs from 'node:fs'
import { readCollectionIndex } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/collection.ts'
import { openCollection } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/collection-source.ts'
import { parseBundle } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/bundle.ts'
import { bundleSource } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/source.ts'
import { moduleScore, moduleScoreAsync } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/score.ts'

const V = 'C:/Users/Lin/AppData/Local/Temp/verify'
const flat = bundleSource(parseBundle(new Uint8Array(fs.readFileSync(`${V}/cond_plain.zip`)).buffer))
const blob = await fs.openAsBlob(`${V}/cond-split/pbmc_cond.zip`)
const coll = await openCollection(blob, await readCollectionIndex(blob))

// The set the Gene sets tab shows first (G2/M-looking, from the tab text).
const SET = ['UBE2C', 'CCNB1', 'CDK1', 'PLK1', 'AURKA', 'CDC20']
const a = moduleScore(flat, SET)
const b = await moduleScoreAsync(coll, SET)
console.log('used   ', JSON.stringify(a.used), '|', JSON.stringify(b.used))
console.log('control genes drawn:', a.control.length, 'vs', b.control.length)
const sameCtrl = a.control.length === b.control.length && a.control.every((g, i) => g === b.control[i])
console.log('control sets identical?', sameCtrl)
if (!sameCtrl) {
  const sa = new Set(a.control), sb = new Set(b.control)
  console.log('  in unsplit only:', a.control.filter(g => !sb.has(g)).length,
              '| in split only:', b.control.filter(g => !sa.has(g)).length)
}
// Scores are per cell in DATASET order, and the two datasets order cells
// differently, so compare by cell id.
const key = c => `${c.x.toFixed(6)},${c.y.toFixed(6)}`
const idA = new Map(flat.d.cells.map((c, i) => [key(c), i]))
let maxAbs = 0, n = 0, sumA = 0, sumB = 0
coll.d.cells.forEach((c, i) => {
  const j = idA.get(key(c))
  if (j === undefined) return
  n++; sumA += a.scores[j]; sumB += b.scores[i]
  maxAbs = Math.max(maxAbs, Math.abs(a.scores[j] - b.scores[i]))
})
console.log(`matched ${n} cells; mean score unsplit ${(sumA / n).toFixed(5)} split ${(sumB / n).toFixed(5)}`)
console.log('max |difference| per cell:', maxAbs.toExponential(3))
const rng = arr => { let lo = Infinity, hi = -Infinity; for (const v of arr) { if (v < lo) lo = v; if (v > hi) hi = v } return [lo, hi] }
console.log('score range unsplit', rng(a.scores).map(v => v.toFixed(4)))
console.log('score range split  ', rng(b.scores).map(v => v.toFixed(4)))
