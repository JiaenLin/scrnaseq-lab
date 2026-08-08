// Build a real collection from the 9.3 GB atlas, using the lab's own modules.
// Identical code to the browser except that the bytes go to a file instead of
// through the save picker, which cannot be driven from a test.
import { createWriteStream } from 'node:fs'
import { writeFile } from 'node:fs/promises'
import h5wasm from 'h5wasm/node'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { chooseSplit, shardScan } from '../src/lib/shard.ts'
import { buildBundle } from '../src/lib/build.ts'
import { collectionPieces, COLLECTION_SCHEMA } from '../src/lib/collection.ts'
import { guessCluster, guessEmbedding, guessRole } from '../src/lib/scan.ts'

const OUT = 'C:/Users/Lin/AppData/Local/Temp/loop/atlas_collection.zip'
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)

await h5wasm.ready
const f = new h5wasm.File('E:/developing_mouse_nervous_system.h5ad', 'r')
let t = Date.now()
const scan = scanOpenH5ad(f, 'developing_mouse_nervous_system.h5ad')
const M = scan.matrices.find(m => m.nnzPerCell)
const nnzPerCell = M.nnzPerCell()
let total = 0; for (const v of nnzPerCell) total += v
log(`scanned in ${Date.now() - t} ms · ${scan.nCells.toLocaleString()} cells · ${total.toLocaleString()} nnz`)

const plan = chooseSplit(scan.columns, nnzPerCell, scan.nCells)
log(`auto-split: ${plan ? `${plan.column.name} -> ${plan.plan.shards.length} parts, ${plan.plan.passes.length} passes` : 'none needed'}`)

const choices = {
  cluster: guessCluster(scan), sample: guessRole(scan, 'sample'),
  condition: guessRole(scan, 'condition'), embedding: guessEmbedding(scan),
  label: 'developing_mouse_nervous_system',
}
log(`choices: cluster=${choices.cluster} sample=${choices.sample} condition=${choices.condition}`)

const parts = []
const order = []
const src = scan.matrices.find(m => m.loadGroups)
t = Date.now()
const SH = plan.plan.shards, PA = plan.plan.passes
for (let pi = 0; pi < PA.length; pi++) {
  const group = PA[pi]
  const mats = src.loadGroups(group.map(s => Int32Array.from(s.cells)))
  for (let gi = 0; gi < group.length; gi++) {
    const s = group[gi]
    const sub = shardScan(scan, { ...s, cells: Int32Array.from(s.cells) },
      new Map([[src.key, mats[gi]]]), `${scan.source} · ${choices.label}`)
    const res = buildBundle(sub, { ...choices, label: `${choices.label} ${s.key}` })
    parts.push({ file: `parts/${s.key.replace(/[^\w.-]+/g, '_')}.zip`, bytes: res.zip })
    order.push(s)
  }
  log(`  pass ${pi + 1}/${PA.length} · ${parts.length}/${SH.length} parts`)
}
log(`converted in ${((Date.now() - t) / 60000).toFixed(1)} min`)

const meta = {
  schema: COLLECTION_SCHEMA, label: choices.label,
  source: `${scan.source} (AnnData, converted in scRNA-seq Lab)`,
  splitBy: plan.column.name, reason: `${(total / 1e6).toFixed(0)} M stored values is more than one piece can hold`,
  nCells: scan.nCells, nGenes: scan.geneSets[scan.matrices[0].geneSet].length,
  parts: parts.map((p, i) => ({
    key: order[i]?.key ?? p.file, file: p.file,
    nCells: order[i]?.cells.length ?? 0, nnz: order[i]?.nnz ?? 0, bytes: p.bytes.length,
  })),
  notes: scan.notes,
}
const out = createWriteStream(OUT)
let bytes = 0
// Hand it over in 64 MB writes, exactly as the app does. A single write of a
// whole 400 MB part fails on Windows with ERR_SYSTEM_ERROR "writev failed" and
// silently truncates the file — which then reads back as "not a collection".
const STEP = 64 << 20
for (const piece of collectionPieces(meta, parts)) {
  for (let at = 0; at < piece.length; at += STEP) {
    const slice = piece.subarray(at, Math.min(piece.length, at + STEP))
    bytes += slice.length
    if (!out.write(Buffer.from(slice))) await new Promise(r => out.once('drain', r))
  }
}
await new Promise(r => out.end(r))
log(`wrote ${OUT} — ${(bytes / 1e9).toFixed(2)} GB, ${parts.length} parts`)
await writeFile(OUT + '.meta.json', JSON.stringify(meta, null, 1))
