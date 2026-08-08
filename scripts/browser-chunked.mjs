// The chunked reader, in a real browser, against a real Blob.
//
// Node proves the arithmetic; it does not prove the thing the studio actually
// does, which is slice bytes out of a File the user dropped on the page and
// inflate them there. Blob.slice().arrayBuffer() is the whole point of storing
// the chunk entry uncompressed, so it is worth seeing it work in chromium
// rather than assuming Node's Blob behaves the same.
//
//   npx vite --port 4241 --strictPort &
//   node scripts/browser-chunked.mjs

import { chromium } from 'playwright-core'

const CHROME = 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe'
const URL = 'http://localhost:4241/'

// Prove the dev server is really compiling this file before trusting anything
// the page says: a 404 that resolved to the SPA fallback would import nothing
// and fail in a way that looks like a hang.
const res = await fetch(`${URL}src/lib/chunked.ts`)
const src = await res.text()
// A string grepped straight out of src/lib/chunked.ts.
const real = res.status === 200 && src.includes('expr.chunkptr is not ascending at chunk')
console.log(`GET /src/lib/chunked.ts -> ${res.status}, ${src.length} bytes, is the real module: ${real}`)
if (!real) { console.log(src.slice(0, 400)); process.exit(1) }

const b = await chromium.launch({ executablePath: CHROME })
const p = await b.newPage()
p.on('console', m => { if (m.type() === 'error') console.log('[console]', m.text().slice(0, 300)) })
p.on('pageerror', e => console.log('[pageerror]', e.message.slice(0, 300)))

// A bare page on the dev server's origin. The lab's own index would drag in the
// whole React app, whose reloads race this check; the module under test is
// still fetched from the real server, transformed exactly as it will be.
const BLANK = `${URL}__chunked-check`
await p.route(BLANK, r => r.fulfill({
  contentType: 'text/html', body: '<!doctype html><title>chunked check</title><h1>chunked check</h1>',
}))
await p.goto(BLANK, { waitUntil: 'domcontentloaded' })
console.log('page text:', (await p.locator('h1').innerText()))

const out = await p.evaluate(async () => {
  const { writeChunked, readGenes, makeChunkCache, CHUNK_GENES } =
    await import('/src/lib/chunked.ts')

  let seed = 999
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
  const nGenes = 500, nCells = 20000
  const indptr = new Int32Array(nGenes + 1)
  const cells = [], values = []
  for (let g = 0; g < nGenes; g++) {
    const want = g % 7 === 0 ? 0 : Math.round(nCells * 0.03 * (0.2 + rnd() * 1.8))
    const picked = new Set()
    while (picked.size < want) picked.add(Math.floor(rnd() * nCells))
    for (const c of [...picked].sort((a, b) => a - b)) {
      cells.push(c); values.push(Math.round(Math.log1p(rnd() * 4000) * 1000) / 1000)
    }
    indptr[g + 1] = cells.length
  }
  const indices = Int32Array.from(cells), data = Float32Array.from(values)
  const { bin, ptr } = writeChunked(indptr, indices, data, CHUNK_GENES)

  // Exactly what the studio will do: the entry sits inside a larger file, and
  // getBytes slices it out. The padding stands in for the rest of the zip.
  const pad = new Uint8Array(1_000_000)
  const file = new Blob([pad, bin, pad])
  const base = pad.length
  const asked = []
  const getBytes = async (from, to) => {
    asked.push(to - from)
    return new Uint8Array(await file.slice(base + from, base + to).arrayBuffer())
  }

  const cache = makeChunkCache()
  const t0 = performance.now()
  const one = await readGenes(getBytes, ptr, indptr, CHUNK_GENES, [321], cache)
  const cold = performance.now() - t0
  const t1 = performance.now()
  await readGenes(getBytes, ptr, indptr, CHUNK_GENES, [321], cache)
  const warm = performance.now() - t1

  let wrong = 0
  const many = Array.from({ length: 120 }, (_, i) => (i * 37) % nGenes)
  for (const v of await readGenes(getBytes, ptr, indptr, CHUNK_GENES, many, cache)) {
    const s = indptr[v.gene], e = indptr[v.gene + 1]
    if (v.cells.length !== e - s) { wrong++; continue }
    for (let i = 0; i < v.cells.length; i++) {
      if (v.cells[i] !== indices[s + i] || v.values[i] !== data[s + i]) { wrong++; break }
    }
  }

  return {
    fileMB: +(file.size / 1e6).toFixed(1),
    entryMB: +(bin.length / 1e6).toFixed(2),
    chunks: ptr.length - 1,
    firstReadBytes: asked[0],
    geneCells: one[0].cells.length,
    coldMs: +cold.toFixed(1),
    warmMs: +warm.toFixed(3),
    genesChecked: many.length,
    wrong,
    inflations: cache.inflations,
  }
})

console.log(out)
const ok = out.wrong === 0 && out.firstReadBytes < out.entryMB * 1e6 / 4 && out.coldMs < 200
console.log(ok
  ? `\nPASS — in chromium: one gene out of a ${out.fileMB} MB Blob read ${out.firstReadBytes} bytes, ${out.coldMs} ms cold, ${out.warmMs} ms warm; ${out.genesChecked} genes exact`
  : '\nFAIL — see the numbers above')
await b.close()
process.exit(ok ? 0 : 1)
