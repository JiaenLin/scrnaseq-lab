import { readFileSync } from 'node:fs'
import { parseBundle, bundleDataset } from '../../scrnaseq-studio/src/lib/bundle.ts'
const f = process.argv[2]
const buf = readFileSync(f)
const b = parseBundle(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength))
const d = bundleDataset(b)
console.log('STUDIO READS IT')
console.log('  label      :', b.meta.label)
console.log('  cells      :', b.meta.nCells.toLocaleString())
console.log('  genes      :', b.genes.length.toLocaleString())
console.log('  nonzeros   :', b.meta.nnz.toLocaleString())
console.log('  clusters   :', b.meta.clusters.length, '—', b.meta.clusters.slice(0, 6).join(', '))
console.log('  conditions :', b.meta.conditions.slice(0, 8).join(', '))
console.log('  samples    :', b.meta.samples.length)
console.log('  embedding  :', b.meta.embedding, '| expression:', b.meta.expression)
console.log('  pseudobulk :', b.pseudobulk ? b.pseudobulk.columns.length + ' columns' : 'none')
console.log('  first cell :', `cluster=${d.cells[0].t} counts=${d.cells[0].counts} xy=${d.cells[0].x.toFixed(2)},${d.cells[0].y.toFixed(2)}`)
console.log('  notes      :'); b.meta.notes.slice(0, 6).forEach(n => console.log('     ·', n))
