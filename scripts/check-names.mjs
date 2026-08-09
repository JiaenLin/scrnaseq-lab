// Developer script — NOT part of `npm test`.
//
// The two format additions, against real objects rather than fixtures: does the
// scan find every embedding and the second naming of the genes, and does a
// bundle built from it carry both? Needs the real files, so it cannot run in CI.
//
// h5wasm/node mounts the real filesystem at the emscripten root, so a path on
// any drive opens without being copied — the same lazy read the app does in a
// worker, with none of the 2 GB ceiling of the copy-it-all path.

import h5wasm from 'h5wasm/node'
import { unzipSync, strFromU8 } from 'fflate'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { buildBundle } from '../src/lib/build.ts'
import { geneNameKind, guessCluster, guessEmbedding, guessRole } from '../src/lib/scan.ts'

const files = process.argv.slice(2)
if (!files.length) files.push('C:/Users/Lin/scrnaseq-studio/testdata/pbmc3k_processed.h5ad')
await h5wasm.ready

for (const path of files) {
  console.log(`\n=== ${path}`)
  const f = new h5wasm.File(path, 'r')
  const scan = scanOpenH5ad(f, path.split(/[\\/]/).pop())
  console.log(`  ${scan.nCells} cells`)
  console.log(`  embeddings: ${scan.embeddings.map(e => `${e.key}(${e.nDims}D)`).join(', ') || 'none'}`)
  for (const [set, genes] of Object.entries(scan.geneSets)) {
    const a = scan.geneAliases[set]
    console.log(`  gene set "${set}": ${genes.length} names, ${geneNameKind(genes)} (${genes[0]})`)
    console.log(`    alias: ${a ? `${a.kind}s from var/${a.column} — ${a.names.slice(0, 3).join(', ')}` : 'none found'}`)
  }

  // Only small objects get built here; an atlas needs the split path.
  if (scan.nCells <= 50000) {
    const res = buildBundle(scan, {
      cluster: guessCluster(scan), sample: guessRole(scan, 'sample'),
      condition: guessRole(scan, 'condition'), embedding: guessEmbedding(scan),
      label: 'names check',
    })
    const z = unzipSync(res.zip)
    console.log(`  bundle entries: ${Object.keys(z).filter(k => k.startsWith('embed') || k.includes('gene')).join(' ')}`)
    console.log(`  meta.embeddings: ${JSON.stringify(res.meta.embeddings)}`)
    console.log(`  meta.geneIdKind: ${res.meta.geneIdKind}  geneAlias: ${JSON.stringify(res.meta.geneAlias)}`)
    if (z['gene_alias.txt']) {
      const alias = strFromU8(z['gene_alias.txt']).split('\n')
      const genes = strFromU8(z['genes.txt']).split('\n')
      console.log(`  alias lines: ${alias.length} (genes ${genes.length}) — ${genes[0]} → ${alias[0]}`)
    }
  }
  f.close()
}
