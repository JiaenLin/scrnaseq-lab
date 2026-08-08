// Open the real atlas collection through the studio's own reader, in node,
// to get a stack for the browser's "Maximum call stack size exceeded".
import { fileBlob } from '../scripts/check-collection.mjs'
import { readCollectionIndex } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/collection.ts'
import { openCollection } from 'file:///C:/Users/Lin/scrnaseq-studio/src/lib/collection-source.ts'

const path = 'C:/Users/Lin/AppData/Local/Temp/verify/atlas_collection.zip'
const blob = await fileBlob(path)
console.log('size', blob.size)
const index = await readCollectionIndex(blob)
console.log('is collection:', !!index, 'parts:', index?.meta.parts.length)
const t0 = Date.now()
try {
  const src = await openCollection(blob, index, p => {
    if (p && p.done % 10 === 0) console.log('  open progress', JSON.stringify(p).slice(0, 120))
  })
  console.log(`opened in ${((Date.now() - t0) / 1000).toFixed(1)}s`)
  console.log('nCells', src.d.nCells, 'genes', src.genes.length, 'clusters', src.clusters.length, 'lazy', src.lazy, 'nParts', src.nParts)
  console.log('first 8 clusters:', src.clusters.slice(0, 8))
  console.log('d.cells length', src.d.cells?.length)
  // The exact call the Cells tab makes.
  const { embedExtent } = await import('file:///C:/Users/Lin/scrnaseq-studio/src/lib/chart.ts')
  try {
    const e = embedExtent(src.d)
    console.log('embedExtent OK', e)
  } catch (err) {
    console.log('embedExtent THREW:', err.message)
  }
} catch (e) {
  console.log('openCollection THREW:', e.message)
  console.log(e.stack.split('\n').slice(0, 14).join('\n'))
}
