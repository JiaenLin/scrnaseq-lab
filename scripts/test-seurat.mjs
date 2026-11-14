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
    check('v5-split: refused, naming the split', /the assay's layers are split/.test(error ?? ''), true)
    check('v5-split: refused with the fix', /obj\[\["RNA"\]\] <- JoinLayers\(obj\[\["RNA"\]\]\)/.test(error ?? ''), true)
    check('v5-split: a layer\'s own count is said', /holds 12 of its assay's 24 cells/.test(error ?? ''), true)
  }

  // Cell names made from integers (as.character(1:24)): R saves them as the
  // numbers and a promise to print them. Read as numbers, the object "had no
  // cell names", the question of whose cells a layer holds was skipped, and a
  // split layer of 12 cells was converted as an object of 12: every cell
  // carried another cell's row of meta.data.
  {
    const { scan, error } = scanOf('seurat-v5-intnames-split.rds', how)
    check('v5-intnames-split: not converted as a part of itself', scan ? `scanned as ${scan.nCells} cells` : 'refused', 'refused')
    check('v5-intnames-split: refused, naming the split', /holds 12 of its assay's 24 cells \(the assay's layers are split\)/.test(error ?? ''), true)
  }
  {
    const want = EXPECTED['seurat-v5-intnames.rds']
    const { scan, error } = scanOf('seurat-v5-intnames.rds', how)
    check('v5-intnames: it scans', error ?? 'ok', 'ok')
    if (scan) {
      check('v5-intnames: all the object\'s cells', scan.nCells, want.cells.length)
      for (const layer of ['counts', 'data']) matchesR('v5-intnames', scan, `RNA@layers$${layer}`, want.layers[layer])
    }
  }

  // A second assay on half the cells, with one layer: it is left out, and what
  // is said of it is that it does not cover the cells, never that it is split.
  {
    const want = EXPECTED['seurat-v5-partial-assay.rds']
    const { scan, error } = scanOf('seurat-v5-partial-assay.rds', how)
    check('v5-partial-assay: it scans', error ?? 'ok', 'ok')
    if (scan) {
      check('v5-partial-assay: all the object\'s cells', scan.nCells, want.cells.length)
      matchesR('v5-partial-assay', scan, 'RNA@layers$counts', want.layers.counts)
      check('v5-partial-assay: the other assay is not offered', scan.matrices.some(m => m.key.startsWith('ADT')), false)
      const said = scan.notes.filter(n => n.startsWith('ADT'))
      check('v5-partial-assay: and is said not to cover the cells',
        said.some(n => /ADT@layers\$counts holds 12 of the object's 24 cells: its assay does not cover the others/.test(n)), true)
      check('v5-partial-assay: nothing is called split', scan.notes.some(n => /split/.test(n)), false)
    }
  }

  // The lines of R in a refusal name the assay that is split. Here it is called
  // Spatial: a line that says RNA is an error in R ("'RNA' not found").
  {
    const { error } = scanOf('seurat-v5-spatial-split.rds', how)
    check('v5-spatial-split: the fix names the assay that is split',
      /obj\[\["Spatial"\]\] <- JoinLayers\(obj\[\["Spatial"\]\]\); saveRDS\(obj, "joined\.rds"\)/.test(error ?? ''), true)
    check('v5-spatial-split: and no assay the object does not have', /RNA/.test(error ?? ''), false)
    // The object after exactly those lines, run in R (make-seurat-fixtures.R).
    const want = EXPECTED['seurat-v5-spatial-joined.rds']
    const j = scanOf('seurat-v5-spatial-joined.rds', how)
    check('v5-spatial-joined: after those lines, it scans', j.error ?? 'ok', 'ok')
    if (j.scan) matchesR('v5-spatial-joined', j.scan, 'Spatial@layers$counts', want.layers.counts)
  }

  // Deferred strings by themselves: every number printed as R prints it.
  {
    const want = EXPECTED['r-deferred-strings.rds']
    const r = open[how]('r-deferred-strings.rds')
    const o = r.read()
    for (const k of ['int', 'dbl', 'rnd', 'hard']) {
      const v = slot(o, k)
      check(`deferred strings: ${k} is read as text`, v?.t, 'str')
      if (v?.t !== 'str') continue
      const bad = []
      v.v.forEach((s, i) => { if (s !== want[k][i]) bad.push(`#${i}: lab ${JSON.stringify(s)}, R ${JSON.stringify(want[k][i])}`) })
      check(`deferred strings: ${k}, ${want[k].length} values, each as R prints it`, [v.v.length, bad.slice(0, 6)], [want[k].length, []])
    }
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
