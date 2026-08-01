// Developer script — NOT part of `npm test`.
//
// Converts the real PBMC 3k objects and checks the result opens in the studio's
// own reader. Needs the scrnaseq-studio repo checked out alongside this one and
// the test objects fetched (see that repo's testdata/fetch.sh), so it cannot run
// in CI — the self-contained suites are test-rds / test-build / test-detect.
import { readFileSync, writeFileSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import h5wasm from 'h5wasm/node'
import { RdsReader, decompressRds } from '../src/lib/rds.ts'
import { scanSeurat } from '../src/lib/seurat.ts'
import { scanH5ad } from '../src/lib/h5ad.ts'
import { buildBundle } from '../src/lib/build.ts'
import { guessCluster, guessRole, guessEmbedding, candidates } from '../src/lib/scan.ts'
// The contract test: the studio's own reader, unmodified.
import { parseBundle, bundleDataset } from '../../scrnaseq-studio/src/lib/bundle.ts'

const T = 'C:/Users/Lin/scrnaseq-studio/testdata'
const show = (scan) => {
  console.log(`  ${scan.nCells} cells · ${scan.columns.length} columns · ${scan.matrices.length} matrices · ${scan.embeddings.length} embeddings`)
  console.log('  matrices:', scan.matrices.map(m => `${m.key}=${m.kind}[${m.nGenes}]`).join(' '))
  console.log('  embeddings:', scan.embeddings.map(e => `${e.key}(${e.nDims})`).join(' '))
  console.log('  cluster candidates:', candidates(scan, 'cluster').map(c => `${c.name}[${c.levels.length}]`).join(' '))
  console.log('  guesses: cluster=%s sample=%s condition=%s embedding=%s',
    guessCluster(scan), guessRole(scan, 'sample'), guessRole(scan, 'condition'), guessEmbedding(scan))
}
const verify = (res, name) => {
  writeFileSync(`C:/Users/Lin/AppData/Local/Temp/${name}.zip`, res.zip)
  const b = parseBundle(res.zip.buffer.slice(res.zip.byteOffset, res.zip.byteOffset + res.zip.byteLength))
  const d = bundleDataset(b)
  console.log(`  STUDIO READS IT: ${b.meta.nCells} cells, ${b.genes.length} genes, ${b.meta.nnz} nnz, clusters=${b.meta.clusters.length}, pb=${b.pseudobulk ? b.pseudobulk.columns.length + ' cols' : 'none'}, ${(res.zip.length/1e6).toFixed(1)} MB`)
  console.log(`  clusters: ${b.meta.clusters.join(', ')}`)
  console.log(`  first cell: cluster=${d.cells[0].t} counts=${d.cells[0].counts} genes=${d.cells[0].genes} xy=${d.cells[0].x.toFixed(2)},${d.cells[0].y.toFixed(2)}`)
  return b
}

console.log('\n=== H5AD ===')
{
  const t = Date.now()
  const scan = await scanH5ad(h5wasm, readFileSync(`${T}/pbmc3k_processed.h5ad`), 'pbmc3k_processed.h5ad')
  console.log(`scanned in ${Date.now()-t} ms`); show(scan)
  const t2 = Date.now()
  const res = buildBundle(scan, { cluster: guessCluster(scan), sample: null, condition: null, embedding: guessEmbedding(scan), label: 'PBMC 3k (Scanpy)' }, s => process.stdout.write(`    ${s}\n`))
  console.log(`built in ${Date.now()-t2} ms`)
  verify(res, 'lab_h5ad')
  console.log('  notes:'); res.meta.notes.forEach(n => console.log('   ·', n))
}

console.log('\n=== SEURAT RDS ===')
{
  const t = Date.now()
  const raw = decompressRds(readFileSync(`${T}/pbmc3k_final.rds`), gunzipSync)
  const r = new RdsReader(raw)
  const scan = scanSeurat(r, r.read(), 'pbmc3k_final.rds')
  console.log(`scanned in ${Date.now()-t} ms`); show(scan)
  const t2 = Date.now()
  const res = buildBundle(scan, { cluster: guessCluster(scan), sample: guessRole(scan,'sample'), condition: null, embedding: guessEmbedding(scan), label: 'PBMC 3k (Seurat)' }, s => process.stdout.write(`    ${s}\n`))
  console.log(`built in ${Date.now()-t2} ms`)
  verify(res, 'lab_rds')
  console.log('  notes:'); res.meta.notes.forEach(n => console.log('   ·', n))
}
