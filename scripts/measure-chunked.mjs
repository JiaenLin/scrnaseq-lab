// How big should a chunk be?
//
// Not a matter of taste: chunk small and deflate has no window to work with, so
// the entry grows; chunk large and reading one gene inflates megabytes for the
// few thousand values it wanted. This measures both ends on a real shard of the
// developing-mouse atlas rather than on synthetic data, because the answer
// depends entirely on how expression values and cell ids actually look.
//
//   node scripts/measure-chunked.mjs [targetNnzMillions]
//
// Needs E:/developing_mouse_nervous_system.h5ad. Not part of `npm test`.

import h5wasm from 'h5wasm/node'
import { deflateSync, inflateSync } from 'fflate'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { planSplit } from '../src/lib/shard.ts'
import { column } from '../src/lib/scan.ts'
import { toGeneMajor } from '../src/lib/matrix.ts'
import { writeChunked, readGenes, makeChunkCache } from '../src/lib/chunked.ts'

const TARGET = (Number(process.argv[2]) || 25) * 1e6
const ms = t => `${((Date.now() - t) / 1000).toFixed(1)} s`

await h5wasm.ready
let t = Date.now()
const f = new h5wasm.File('E:/developing_mouse_nervous_system.h5ad', 'r')
const scan = scanOpenH5ad(f, 'developing_mouse_nervous_system.h5ad')
const M = scan.matrices[0]
const nnzPer = M.nnzPerCell()
console.log(`scanned ${scan.nCells.toLocaleString()} cells in ${ms(t)}`)

const plan = planSplit(column(scan, 'dissection'), nnzPer, scan.nCells)
// The shard nearest the target: representative of what the studio will open.
const shard = plan.shards.reduce((a, b) =>
  Math.abs(b.nnz - TARGET) < Math.abs(a.nnz - TARGET) ? b : a)
console.log(`shard "${shard.key}": ${shard.cells.length.toLocaleString()} cells, ` +
  `${(shard.nnz / 1e6).toFixed(1)} M nonzeros`)

t = Date.now()
const [cm] = M.loadGroups([shard.cells], fr => process.stdout.write(`\r  reading ${(fr * 100).toFixed(1)}%   `))
console.log(`\n  one forward pass: ${ms(t)}`)
const csc = toGeneMajor(cm)
const nGenes = csc.indptr.length - 1
const nnz = csc.data.length
console.log(`gene-major: ${nGenes.toLocaleString()} genes, ${nnz.toLocaleString()} values, ` +
  `${(nnz * 8 / 1e6).toFixed(0)} MB raw\n`)

// The baseline: what the bundle already pays for these two arrays today.
const flat = deflateSync(new Uint8Array(csc.indices.buffer, 0, nnz * 4), { level: 6 }).length +
  deflateSync(new Uint8Array(csc.data.buffer, 0, nnz * 4), { level: 6 }).length
console.log(`flat deflated expr.indices + expr.data: ${(flat / 1e6).toFixed(1)} MB\n`)

console.log('chunkGenes   entry MB   vs flat   mean chunk KB   max chunk KB   write s   inflate 1 chunk ms')
const rows = []
for (const cg of [16, 32, 64, 128, 256, 512]) {
  const t0 = Date.now()
  const { bin, ptr } = writeChunked(csc.indptr, csc.indices, csc.data, cg)
  const write = (Date.now() - t0) / 1000
  const sizes = []
  for (let k = 1; k < ptr.length; k++) sizes.push(ptr[k] - ptr[k - 1])
  const mean = sizes.reduce((a, b) => a + b, 0) / sizes.length
  const max = Math.max(...sizes)

  // Inflating a mid-sized chunk is what a gene lookup costs.
  const order = sizes.map((s, k) => [s, k]).sort((a, b) => a[0] - b[0])
  const mid = order[Math.floor(order.length / 2)][1]
  const raw = bin.subarray(ptr[mid], ptr[mid + 1])
  let t1 = Date.now(), reps = 0
  while (Date.now() - t1 < 300) { inflateSync(raw); reps++ }
  const inf = (Date.now() - t1) / reps

  rows.push({ cg, bytes: bin.length, mean, max, write, inf })
  console.log(`${String(cg).padStart(9)}   ${(bin.length / 1e6).toFixed(1).padStart(8)}   ` +
    `${((bin.length / flat - 1) * 100).toFixed(1).padStart(6)}%   ${(mean / 1e3).toFixed(0).padStart(13)}   ` +
    `${(max / 1e3).toFixed(0).padStart(12)}   ${write.toFixed(1).padStart(7)}   ${inf.toFixed(1).padStart(18)}`)
}

// What the gap encoding is worth: undo it and store the cell ids as they are.
{
  const cg = 64
  const plain = new Int32Array(nnz)
  for (let g = 0; g < nGenes; g++) {
    // A monotonic run with the same value distribution, so this measures the
    // encoding rather than the data.
    for (let k = csc.indptr[g]; k < csc.indptr[g + 1]; k++) plain[k] = csc.indices[k]
  }
  // writeChunked gaps whatever it is given, so feed it a pre-summed array to
  // cancel that out: the ids come out stored verbatim.
  const undone = new Int32Array(nnz)
  for (let g = 0; g < nGenes; g++) {
    let acc = 0
    for (let k = csc.indptr[g]; k < csc.indptr[g + 1]; k++) { acc += plain[k]; undone[k] = acc }
  }
  const withGaps = writeChunked(csc.indptr, csc.indices, csc.data, cg).bin.length
  const noGaps = writeChunked(csc.indptr, undone, csc.data, cg).bin.length
  console.log(`\ncell-id gaps at chunkGenes=${cg}: ${(withGaps / 1e6).toFixed(1)} MB ` +
    `vs ${(noGaps / 1e6).toFixed(1)} MB storing ids verbatim ` +
    `(${((withGaps / noGaps - 1) * 100).toFixed(1)}%)`)
}

// What one gene costs through the real reader.
{
  const cg = 64
  const { bin, ptr } = writeChunked(csc.indptr, csc.indices, csc.data, cg)
  const get = async (from, to) => bin.subarray(from, to)
  const busiest = (() => {
    let best = 0
    for (let g = 0; g < nGenes; g++) if (csc.indptr[g + 1] - csc.indptr[g] > csc.indptr[best + 1] - csc.indptr[best]) best = g
    return best
  })()

  const cold = async g => { const c = makeChunkCache(); const t0 = performance.now(); await readGenes(get, ptr, csc.indptr, cg, [g], c); return performance.now() - t0 }
  const cache = makeChunkCache()
  await readGenes(get, ptr, csc.indptr, cg, [busiest], cache)
  const t0 = performance.now()
  await readGenes(get, ptr, csc.indptr, cg, [busiest], cache)
  const warm = performance.now() - t0

  // Real data, not a generator: 200 genes checked value for value.
  let wrong = 0
  const sample = Array.from({ length: 200 }, (_, i) => (i * 149) % nGenes)
  for (const v of await readGenes(get, ptr, csc.indptr, cg, sample)) {
    const s = csc.indptr[v.gene], e = csc.indptr[v.gene + 1]
    if (v.cells.length !== e - s) { wrong++; continue }
    for (let i = 0; i < v.cells.length; i++) {
      if (v.cells[i] !== csc.indices[s + i] || v.values[i] !== csc.data[s + i]) { wrong++; break }
    }
  }
  console.log(`\n${sample.length} real genes checked value for value: ${wrong} wrong`)

  const colds = []
  for (let i = 0; i < 20; i++) colds.push(await cold((i * 1499) % nGenes))
  colds.sort((a, b) => a - b)
  console.log(`\nread one gene at chunkGenes=${cg}: cold median ${colds[10].toFixed(1)} ms ` +
    `(worst ${colds[19].toFixed(1)} ms), warm ${warm.toFixed(2)} ms`)
  console.log(`busiest gene has ${(csc.indptr[busiest + 1] - csc.indptr[busiest]).toLocaleString()} cells`)
}

f.close()
