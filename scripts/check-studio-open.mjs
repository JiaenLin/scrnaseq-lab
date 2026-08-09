// Developer script — NOT part of `npm test`.
//
// The last mile: hand a finished collection to the studio's own openCollection
// and see it come back as one dataset. The other checks read parts one at a
// time; this is the merge, which is what the studio actually does on open, and
// it is the path that would notice if the new entries confused the reader.
//
//   node scripts/check-studio-open.mjs <collection.zip>

import { statSync } from 'node:fs'
import { open } from 'node:fs/promises'
import { openCollection } from '../../scrnaseq-studio/src/lib/collection-source.ts'
import { readCollectionIndex } from '../src/lib/collection.ts'

const path = process.argv[2]
const size = statSync(path).size
const fh = await open(path, 'r')

/** Just enough Blob: a size, and ranges read at the offset asked for. */
const blob = (from = 0, to = size) => ({
  size: to - from,
  slice: (a, b = to - from) => blob(from + a, from + Math.min(b, to - from)),
  arrayBuffer: async () => {
    const buf = Buffer.alloc(to - from)
    await fh.read(buf, 0, to - from, from)
    return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)
  },
})

const t = Date.now()
const file = blob()
const index = await readCollectionIndex(file)
if (!index) throw new Error('this is not a collection')
let lastPhase = ''
const src = await openCollection(file, index, (phase, done, total) => {
  if (phase !== lastPhase) { lastPhase = phase; console.log(`  ${phase} (${total})`) }
  void done
})
const d = src.d
console.log(`opened in ${((Date.now() - t) / 1000).toFixed(1)} s`)
console.log(`  ${d.label} — ${d.nCells.toLocaleString()} cells, ${src.genes.length.toLocaleString()} genes`)
console.log(`  clusters ${d.cells.reduce((m, c) => Math.max(m, c.t), 0) + 1}, samples ${d.samples.length}, conditions ${d.conds.length}`)
console.log(`  first cell: cluster=${d.cells[0].t} xy=${d.cells[0].x.toFixed(2)},${d.cells[0].y.toFixed(2)} counts=${d.cells[0].counts}`)
console.log(`  meta.embedding=${src.meta.embedding}`)
console.log(`  notes: ${src.meta.notes.length}`)
await fh.close()
