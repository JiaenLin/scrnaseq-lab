// Seurat objects written by R itself, read by the lab, compared with what R
// says of them.
//
// The fixtures and seurat-expected.json come from scripts/make-seurat-fixtures.R
// (saveRDS of objects Seurat 5 builds; the names and sums are R's own). The
// other reader tests write R's format by hand, so they can only agree with the
// reader's author. One real Seurat 5 object showed two things that had never
// met the reader:
//   - a function among a logged command's parameters, whose namespace R writes
//     as a bare string vector. The walk came out one item short, the object's
//     `tools` and `class` were never read, and the page said "not an S4 object";
//   - an Assay5, whose layers carry no names (the assay's two maps do), so no
//     matrix "had gene names".
// Each object is read twice: whole, and as the page reads it, through a stream
// fed 64 compressed bytes at a time.

import { readFileSync } from 'fs'
import { gunzipSync } from 'zlib'
import { fileURLToPath } from 'url'
import { RdsReader, classOf, decompressRds, namesOf, slot } from '../src/lib/rds.ts'
import { GzWindow } from '../src/lib/gzwindow.ts'
import { scanSeurat } from '../src/lib/seurat.ts'

const here = p => fileURLToPath(new URL(p, import.meta.url))
const EXPECTED = JSON.parse(readFileSync(here('./fixtures/seurat-expected.json'), 'utf8'))

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  const cut = v => { const s = JSON.stringify(v); return s.length > 300 ? s.slice(0, 300) + '…' : s }
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${cut(got)}\n        want ${cut(want)}`}`)
}

/** The two ways the lab opens an .rds. */
const open = {
  whole: name => new RdsReader(decompressRds(readFileSync(here(`./fixtures/${name}`)), gunzipSync)),
  streamed: name => {
    const raw = readFileSync(here(`./fixtures/${name}`))
    const mk = () => new GzWindow(raw.length, (a, b) => raw.subarray(a, b), () => {}, 64)
    return new RdsReader(mk(), mk)
  },
}

const scanOf = (name, how) => {
  const r = open[how](name)
  let obj
  // A reader that refuses the file while walking it is a result too, not a crash.
  try { obj = r.read() } catch (e) { return { r, obj: undefined, error: `read: ${e.message}` } }
  try { return { r, obj, scan: scanSeurat(r, obj, name) } }
  catch (e) { return { r, obj, error: String(e.message) } }
}

const sums = m => {
  const col = new Array(m.nCells).fill(0), row = new Array(m.nGenes).fill(0)
  for (let c = 0; c < m.nCells; c++) {
    for (let k = m.indptr[c]; k < m.indptr[c + 1]; k++) { col[c] += m.data[k]; row[m.indices[k]] += m.data[k] }
  }
  return { col, row }
}
// The bundle stores float32, so a normalised value agrees with R's double to
// about seven digits; counts are whole numbers and agree exactly.
const near = (a, b) => a.length === b.length && a.every((v, i) => Math.abs(v - b[i]) <= 1e-4 * Math.max(1, Math.abs(b[i])))

/** A scanned matrix against the same layer as R reads it. */
const matchesR = (label, scan, key, want) => {
  const m = scan.matrices.find(x => x.key === key)
  check(`${label}: ${key} is offered`, !!m, true)
  if (!m) return
  check(`${label}: ${key} has R's gene names, in R's order`, scan.geneSets[m.geneSet], want.genes)
  check(`${label}: ${key} has R's shape`, [m.nGenes, m.nCells], want.dim)
  const got = sums(m.load())
  check(`${label}: ${key} per-cell sums are R's`, near(got.col, want.col_sums), true)
  check(`${label}: ${key} per-gene sums are R's`, near(got.row, want.row_sums), true)
}

for (const how of ['whole', 'streamed']) {
  console.log(`\nREAD ${how.toUpperCase()}`)

  for (const name of ['seurat-v5-plain.rds', 'seurat-v5-functions.rds', 'seurat-v4-functions.rds']) {
    const want = EXPECTED[name]
    const { obj, scan, error } = scanOf(name, how)
    const label = name.replace(/^seurat-|\.rds$/g, '')
    check(`${label}: the object's class is read`, classOf(obj), [want.class])
    check(`${label}: its commands keep their own names`, namesOf(slot(obj, 'commands')), want.commands)
    check(`${label}: it scans`, error ?? 'ok', 'ok')
    if (!scan) continue
    check(`${label}: the cells are the object's`, scan.nCells, want.cells.length)
    const v5 = want.assay_class === 'Assay5'
    for (const layer of ['counts', 'data']) {
      matchesR(label, scan, v5 ? `RNA@layers$${layer}` : `RNA@${layer}`, want.layers[layer])
    }
    // A layer holding three of the forty genes, handed to R in another order:
    // its rows are named from the assay's map, in the assay's order.
    if (want.layers.sub) matchesR(label, scan, 'RNA@layers$sub', want.layers.sub)
    check(`${label}: scale.data is not offered`, scan.matrices.some(m => /scale\.data/.test(m.key)), false)
    check(`${label}: the reductions are found`, scan.embeddings.map(e => e.key).sort(), [...want.reductions].sort())
    for (const c of ['sample', 'cell_type']) {
      check(`${label}: meta.data column ${c} is read`, scan.columns.some(x => x.name === c), true)
    }
  }

  // Layers split per sample: no layer holds the object's cells, and the page
  // says so with the two lines of R that join them.
  {
    const { obj, error } = scanOf('seurat-v5-split.rds', how)
    check('v5-split: the object\'s class is read', classOf(obj), ['Seurat'])
    check('v5-split: refused, naming the split', /layers are split, one per sample/.test(error ?? ''), true)
    check('v5-split: refused with the fix', /JoinLayers\(obj\[\["RNA"\]\]\)/.test(error ?? ''), true)
    check('v5-split: a layer\'s own count is said', /holds 12 of the object's 24 cells/.test(error ?? ''), true)
  }

  // meta.data in another order than the layers' cells: a layer is then not the
  // object's matrix, and reading it as one would give every cell another's row.
  {
    const { error } = scanOf('seurat-v5-reordered.rds', how)
    check('v5-reordered: refused, naming the order', /in another order than meta\.data/.test(error ?? ''), true)
  }
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll Seurat tests passed\n')
process.exit(failed ? 1 : 0)
