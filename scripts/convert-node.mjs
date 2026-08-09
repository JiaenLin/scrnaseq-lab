// Developer script — NOT part of `npm test`.
//
// The whole conversion, headless: the same modules the app runs in its worker,
// driven from node instead of from a page. It exists because the browser path
// ends at `showSaveFilePicker`, which no automation can answer, and because a
// quarter-hour conversion of a 9.3 GB atlas is worth being able to repeat from
// a shell.
//
// Nothing here reimplements the format. scan → chooseSplit → buildBundle →
// collectionPieces are the app's own functions in the app's own order; if this
// script and the app ever disagree, this script is wrong.
//
//   node --max-old-space-size=8192 scripts/convert-node.mjs <in.h5ad> <out.zip> [label]
//                                                          [--extra <obs column>]…

import h5wasm from 'h5wasm/node'
import { createWriteStream, mkdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs'
import { stat } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { scanOpenH5ad } from '../src/lib/h5ad.ts'
import { buildBundle, streamableMatrix } from '../src/lib/build.ts'
import { chooseSplit, shardScan } from '../src/lib/shard.ts'
import { guessCluster, guessEmbedding, guessRole } from '../src/lib/scan.ts'
import { COLLECTION_SCHEMA, INDEX_NAME, collectionPieces } from '../src/lib/collection.ts'

const argv = process.argv.slice(2)
// The columns to carry beyond the three roles, by name. Same as the lab's
// Extras card, which is empty by default for the same reason: an object's other
// annotations are worth offering and not worth assuming.
const extras = argv.flatMap((a, i) => (a === '--extra' && argv[i + 1] ? [argv[i + 1]] : []))
const positional = argv.filter((a, i) => a !== '--extra' && argv[i - 1] !== '--extra')
const [input, output, labelArg] = positional
if (!input || !output) {
  console.error('usage: node scripts/convert-node.mjs <in.h5ad> <out.zip> [label] [--extra <col>]…')
  process.exit(2)
}
const label = labelArg ?? input.split(/[\\/]/).pop().replace(/\.(h5ad|h5)$/i, '')
const t0 = Date.now()
const log = (...a) => console.log(`[${((Date.now() - t0) / 1000).toFixed(0).padStart(4)}s]`, ...a)
const safe = s => (s || 'part').replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 80)

await h5wasm.ready
log('opening', input)
const f = new h5wasm.File(input, 'r')
const scan = scanOpenH5ad(f, input.split(/[\\/]/).pop())
log(`${scan.nCells.toLocaleString()} cells · ${scan.columns.length} columns · ${scan.matrices.length} matrices`)
log('embeddings:', scan.embeddings.map(e => `${e.key}(${e.nDims}D)`).join(', '))
for (const [set, genes] of Object.entries(scan.geneSets)) {
  const a = scan.geneAliases[set]
  log(`gene set ${set}: ${genes.length} names (${genes[0]}) · alias ${a ? `${a.kind}s from var/${a.column}` : 'none'}`)
}

const choices = {
  cluster: guessCluster(scan),
  sample: guessRole(scan, 'sample'),
  condition: guessRole(scan, 'condition'),
  embedding: guessEmbedding(scan),
  label,
  extras,
}
log('choices:', JSON.stringify(choices))
for (const name of extras) {
  const c = scan.columns.find(x => x.name === name)
  if (!c) throw new Error(`this object has no obs column called "${name}"`)
  log(`carrying ${name} as an extra column: ${c.levels.length} levels`)
}

// The app measures the matrix before deciding; so does this.
const measurable = scan.matrices.find(m => m.nnzPerCell)
const nnzPerCell = measurable ? measurable.nnzPerCell() : null
const split = chooseSplit(scan.columns, nnzPerCell, scan.nCells)
log(split ? `split along ${split.column.name}: ${split.plan.shards.length} parts in ${split.plan.passes.length} passes` : 'no split needed')

// Parts are parked on disk, not kept in memory.
//
// The app holds every part at once — it must, because the browser has nowhere
// else to put them — and on an atlas that is 5.8 GB of live typed arrays. A
// desktop tab with the machine to itself manages it; a node process sharing a
// 15 GB machine with half a dozen other things does not, and what it does
// instead is page to disk until everything on the machine crawls. So each part
// goes to a scratch file as it is written and comes back one at a time when the
// container is assembled. The bytes are identical either way.
const scratch = `${tmpdir()}/lab-convert-${process.pid}`
mkdirSync(scratch, { recursive: true })
const parked = []
const partInfo = []
const used = new Set()
const keep = (key, res) => {
  let file = `parts/${safe(key)}.zip`
  for (let i = 2; used.has(file); i++) file = `parts/${safe(key)}-${i}.zip`
  used.add(file)
  const tmp = `${scratch}/${parked.length}.zip`
  writeFileSync(tmp, res.zip)
  parked.push({ file, tmp })
  partInfo.push({ key, file, nCells: res.meta.nCells, nnz: res.meta.nnz, bytes: res.zip.length })
  log(`  wrote ${file} — ${res.meta.nCells.toLocaleString()} cells, ${(res.zip.length / 1e6).toFixed(1)} MB`)
  return res.meta
}

/**
 * The part list `collectionPieces` expects, reading one part at a time.
 *
 * It touches `bytes` several times per part — for the CRC, the two size fields
 * and the payload — and always finishes with one before starting the next, so
 * a cache of exactly one entry turns those into a single read.
 */
function lazyParts(list) {
  let atIndex = -1
  let held = null
  return list.map((p, i) => ({
    file: p.file,
    get bytes() {
      if (atIndex !== i) { held = readFileSync(p.tmp); atIndex = i }
      return held
    },
  }))
}

let firstMeta
if (!split) {
  firstMeta = keep(label, buildBundle(scan, choices, s => log('   ' + s)))
} else {
  const src = streamableMatrix(scan, s => { scan.notes.push(s); log('   ' + s) })
  if (!src) throw new Error('no matrix in this object can be read in parts')
  log(`streaming ${src.key}`)
  for (let pi = 0; pi < split.plan.passes.length; pi++) {
    const group = split.plan.passes[pi]
    log(`pass ${pi + 1} of ${split.plan.passes.length} — ${group.length} parts`)
    let last = -1
    const mats = src.loadGroups(group.map(s => Int32Array.from(s.cells)), frac => {
      const pct = Math.floor(frac * 10)
      if (pct !== last) { last = pct; log(`   reading ${(frac * 100).toFixed(0)}%`) }
    })
    for (let gi = 0; gi < group.length; gi++) {
      const s = group[gi]
      const sub = shardScan(scan, { ...s, cells: Int32Array.from(s.cells) },
        new Map([[src.key, mats[gi]]]), `${scan.source} · ${label}`)
      const meta = keep(s.key, buildBundle(sub, { ...choices, label: `${label} ${s.key}` }))
      firstMeta ??= meta
      // The shard's matrix is dead once its bundle exists; a pass can hold a
      // gigabyte across its shards, so drop each as it is used.
      mats[gi] = null
    }
  }
}
f.close()

const meta = {
  schema: COLLECTION_SCHEMA,
  label,
  source: firstMeta.source,
  splitBy: split ? split.column.name : null,
  reason: split ? split.reason : null,
  nCells: partInfo.reduce((n, p) => n + p.nCells, 0),
  nGenes: firstMeta.nGenes,
  clusterOrder: scan.columns.find(c => c.name === choices.cluster)?.levels,
  // The whole object's order for every categorical column a part can only carry
  // part of. Without it the studio reconstructs what the parts imply and falls
  // back to collation for the rest, which compares the digit run after the dot
  // as a number and puts e16.5 before e16.25.
  condOrder: scan.columns.find(c => c.name === choices.condition)?.levels,
  extraOrder: Object.fromEntries(extras.flatMap(name => {
    const levels = scan.columns.find(c => c.name === name)?.levels
    return levels ? [[name, levels]] : []
  })),
  embeddings: firstMeta.embeddings?.map(e => e.key),
  geneIdKind: firstMeta.geneIdKind,
  geneAlias: firstMeta.geneAlias
    ? { kind: firstMeta.geneAlias.kind, column: firstMeta.geneAlias.column }
    : null,
  parts: partInfo,
  notes: firstMeta.notes,
}

log(`writing ${output}`)
const out = createWriteStream(output)
try {
  for (const piece of collectionPieces(meta, lazyParts(parked))) {
    if (!out.write(piece)) await new Promise(r => out.once('drain', r))
  }
  await new Promise((res, rej) => { out.on('error', rej); out.end(res) })
} finally {
  rmSync(scratch, { recursive: true, force: true })
}
log(`done — ${((await stat(output)).size / 1e9).toFixed(2)} GB in ${partInfo.length} parts`)
console.log(JSON.stringify({ ...meta, parts: `${partInfo.length} parts` }, null, 1))
void INDEX_NAME
