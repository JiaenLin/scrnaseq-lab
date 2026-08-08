import h5wasm from 'h5wasm/node'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { splitCandidates, planSplit, shardScan, shardFileName } from '../src/lib/shard.ts'
import { buildBundle } from '../src/lib/build.ts'
import { column } from '../src/lib/scan.ts'
import { parseBundle, bundleDataset } from '../../scrnaseq-studio/src/lib/bundle.ts'

await h5wasm.ready
const t0 = Date.now()
const f = new h5wasm.File('E:/developing_mouse_nervous_system.h5ad', 'r')
const scan = scanOpenH5ad(f, 'developing_mouse_nervous_system.h5ad')
console.log(`scanned 9.3 GB in ${Date.now() - t0} ms`)
console.log(`${scan.nCells.toLocaleString()} cells - ${scan.columns.length} columns - matrices:`,
  scan.matrices.map(m => `${m.key}=${m.kind}`).join(' '), '- embeddings:', scan.embeddings.map(e => e.key).join(' '))

const M = scan.matrices[0]
console.log('nnzPerCell available:', !!M.nnzPerCell, '| loadGroups available:', !!M.loadGroups)
const t1 = Date.now()
const nnzPer = M.nnzPerCell()
console.log(`row offsets read in ${Date.now() - t1} ms; total nnz =`,
  nnzPer.reduce((a, b) => a + b, 0).toLocaleString())

console.log('\nSPLIT CANDIDATES (exact, no values read):')
for (const { column: c, plan } of splitCandidates(scan.columns, nnzPer, scan.nCells).slice(0, 6)) {
  console.log(`  ${c.name.padEnd(34)} ${String(plan.shards.length).padStart(3)} shards, ` +
    `largest ${(plan.shards[0].nnz / 1e6).toFixed(1)}M nnz, ` +
    `${plan.passes.length} pass(es), ${plan.oversized.length} too big`)
}

// Convert the two smallest shards of a chosen split, to prove the pipeline.
const col = column(scan, 'dissection')
const plan = planSplit(col, nnzPer, scan.nCells)
const pick = plan.shards.slice(-2)
console.log(`\nCONVERTING ${pick.length} shards of "${plan.column}":`,
  pick.map(s => `${s.key} (${s.cells.length} cells, ${(s.nnz/1e6).toFixed(1)}M nnz)`).join(', '))

const t2 = Date.now()
const mats = M.loadGroups(pick.map(s => s.cells), fr => process.stdout.write(`\r  reading ${(fr*100).toFixed(1)}%   `))
console.log(`\n  one forward pass: ${((Date.now() - t2)/1000).toFixed(1)} s`)

for (let i = 0; i < pick.length; i++) {
  const s = pick[i]
  const sc = shardScan(scan, s, new Map([[M.key, mats[i]]]), `atlas / ${plan.column} = ${s.key}`)
  const res = buildBundle(sc, {
    cluster: 'Subclass', sample: 'sample_id', condition: 'Age',
    embedding: scan.embeddings[0].key, label: `atlas ${s.key}`,
  })
  const b = parseBundle(res.zip.buffer.slice(res.zip.byteOffset, res.zip.byteOffset + res.zip.byteLength))
  const d = bundleDataset(b)
  console.log(`  ${shardFileName('atlas', s.key).padEnd(42)} ${(res.zip.length/1e6).toFixed(1)} MB  ` +
    `-> STUDIO READS: ${b.meta.nCells} cells, ${b.genes.length} genes, ${b.meta.nnz.toLocaleString()} nnz, ` +
    `${b.meta.clusters.length} clusters, ${d.conds.length} conditions`)
}
