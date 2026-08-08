import h5wasm from 'h5wasm/node'
import { readFileSync } from 'node:fs'
import { scanH5ad } from '../src/lib/h5ad.ts'
await h5wasm.ready
const f = process.argv[2]
try {
  const r = await scanH5ad(h5wasm, readFileSync(f), f.split(/[\/]/).pop()); const scan = r.scan ?? r
  console.log('OK. cells', scan.nCells, 'matrices', scan.matrices.map(m => `${m.key}:${m.kind}`))
  console.log('columns', scan.columns.map(c => `${c.name}(${c.kind}${c.levels?.length ? ',' + c.levels.length : ''})`).join(' '))
  console.log('embeddings', scan.embeddings.map(e => e.key))
} catch (e) {
  console.log('THREW:', e.message)
  console.log(e.stack.split('\n').slice(0, 12).join('\n'))
}
