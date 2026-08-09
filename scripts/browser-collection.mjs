// Developer script — NOT part of `npm test`.
//
// Drives the whole lab in a real browser and then takes the downloaded file
// apart: one collection.zip out, every part inside readable by the studio's own
// parseBundle, and the chunked entry agreeing with the flat one gene for gene.
//
//   npm run build && npx vite preview --port 4242 --strictPort
//   node scripts/browser-collection.mjs [file] [timeout-minutes] [extra,columns]
//
// The container is read through a file handle, never loaded: a real atlas
// collection is gigabytes, and the point of the format is that a reader touches
// only the ranges it needs. Reading it whole here would prove nothing.
//
// Not node:fs openAsBlob either — it reports `size` modulo 2^32, so a 5.83 GB
// file comes back claiming 1.53 GB and every offset past 4 GB lands in the
// wrong place. That cost a run: the container was perfect and the reader said
// "this is not a collection".
//
// It asserts only on strings grepped from the source, and prints the page text
// whenever a selector misses — a silent false negative here is worse than a
// crash, because it looks exactly like a hang.
import { createWriteStream, mkdirSync, statSync } from 'node:fs'
import { createServer } from 'node:http'
import { chromium } from 'playwright-core'
import { checkCollection } from './check-collection.mjs'

const EXE = 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe'
const OUT = 'C:/Users/Lin/AppData/Local/Temp/labshots'
const FILE = process.argv[2] ?? 'C:/Users/Lin/scrnaseq-studio/testdata/pbmc3k_processed.h5ad'
const MINUTES = Number(process.argv[3] ?? 10)
/** Columns to tick on the Extras card, by name. Nothing is ticked by default. */
const EXTRAS = (process.argv[4] ?? '').split(',').map(s => s.trim()).filter(Boolean)
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`  ${ok ? 'ok  ' : 'FAIL'} ${name}${ok ? '' : `\n         got  ${JSON.stringify(got)}\n         want ${JSON.stringify(want)}`}`)
}

mkdirSync(OUT, { recursive: true })

// A stand-in for the operating system's save dialog.
//
// The app writes multi-gigabyte results through showSaveFilePicker, because a
// blob: URL download is cancelled by Chromium somewhere past 256 MB. Playwright
// cannot drive a native picker, so the picker — and only the picker — is
// replaced here: the blob is still read by the app's own streaming code, and
// the bytes still land in a real file, they just arrive over a socket. Sending
// them through CDP instead would take longer than the conversion.
const SINK = 4343
const SINK_FILE = `${OUT}/sink.zip`
let landed
const sinkPath = new Promise(resolve => { landed = resolve })
let out = null
const cors = {
  'access-control-allow-origin': '*',
  'access-control-allow-methods': 'POST, OPTIONS',
  'access-control-allow-headers': '*',
}
const server = createServer((req, res) => {
  if (req.method === 'OPTIONS') { res.writeHead(204, cors); res.end(); return }
  if (req.url === '/close') {
    out?.end(() => landed(SINK_FILE))
    res.writeHead(200, cors); res.end('ok')
    return
  }
  // One slice per request, appended in the order the app wrote them — which is
  // what a real writable stream does with the same calls.
  out ??= createWriteStream(SINK_FILE)
  req.pipe(out, { end: false })
  req.on('end', () => { res.writeHead(200, cors); res.end('ok') })
}).listen(SINK)
// A slice can still take a while over a socket, and Node's default five-minute
// request timeout aborting would look exactly like the app failing to write.
server.requestTimeout = 0
server.headersTimeout = 0

const b = await chromium.launch({ executablePath: EXE })
const page = await b.newPage({ viewport: { width: 1200, height: 1400 } })
await page.addInitScript(port => {
  window.showSaveFilePicker = async () => ({
    createWritable: async () => ({
      async write(data) {
        const r = await fetch(`http://localhost:${port}/w`, { method: 'POST', body: data })
        if (!r.ok) throw new Error(`the sink refused the write: ${r.status}`)
      },
      async close() {
        await fetch(`http://localhost:${port}/close`, { method: 'POST' })
      },
    }),
  })
}, SINK)
const errs = []
page.on('pageerror', e => errs.push('PAGEERROR ' + e.message))
page.on('console', m => { if (m.type() === 'error') errs.push('CONSOLE ' + m.text()) })

const say = async () => (await page.locator('body').innerText()).replace(/\n+/g, ' | ').slice(0, 1500)

await page.goto('http://localhost:4242/', { waitUntil: 'networkidle' })
log('handing over', FILE, `(${(statSync(FILE).size / 1e9).toFixed(2)} GB)`)

const t0 = Date.now()
await page.setInputFiles('input[type=file]', FILE)

// Whichever comes first: the questions, or a refusal.
const scanned = await Promise.race([
  page.waitForSelector('text=Which column is the cell type annotation?', { timeout: MINUTES * 60_000 })
    .then(() => 'ok'),
  page.waitForSelector('.note:not(.note-info) b', { timeout: MINUTES * 60_000 })
    .then(() => 'error'),
]).catch(() => 'neither')
log(`scan finished in ${((Date.now() - t0) / 1000).toFixed(1)} s -> ${scanned}`)
if (scanned !== 'ok') {
  console.log('PAGE TEXT:', await say())
  await b.close()
  process.exit(1)
}

// The split UI is gone. Nothing on this page may offer a choice about it.
const text = await page.locator('body').innerText()
check('no split card', /Split it into several|This object is too large for one bundle/.test(text), false)
check('no "do not split" option', /Do not split/.test(text), false)

// The extra columns, ticked the way a person would. Scoped to that card: the
// role pickers draw the same button for the same column, and clicking one of
// those would change the answer to a different question.
const EXTRA_CARD = 'section.card:has-text("Anything else worth breaking a figure down by")'
for (const name of EXTRAS) {
  const btn = page.locator(`${EXTRA_CARD} button:has(span.mono:text-is("${name}"))`)
  check(`${name} is offered as an extra column`, await btn.count(), 1)
  await btn.first().click()
}
await page.screenshot({ path: `${OUT}/collection-review.png`, fullPage: true })

log('converting...')
const t1 = Date.now()
await page.click('button.btn-primary:has-text("Convert")')

// Say what it is doing while it does it, so a 40-minute run is not a blank wall.
const tick = setInterval(async () => {
  try {
    const s = await page.locator('.btn-primary').first().innerText()
    log(`  ${((Date.now() - t1) / 60_000).toFixed(1)} min · ${s.replace(/\n/g, ' ')}`)
  } catch { /* the button is gone: it finished */ }
}, 60_000)

const done = await Promise.race([
  page.waitForSelector("text=the lab's job is done", { timeout: MINUTES * 60_000 }).then(() => 'ok'),
  page.waitForSelector('text=That did not convert.', { timeout: MINUTES * 60_000 }).then(() => 'error'),
]).catch(() => 'neither')
clearInterval(tick)
log(`convert finished in ${((Date.now() - t1) / 60_000).toFixed(1)} min -> ${done}`)
if (done !== 'ok') {
  console.log('PAGE TEXT:', await say())
  await b.close()
  process.exit(1)
}
console.log('FINISH SCREEN:', await say())
await page.screenshot({ path: `${OUT}/collection-done.png`, fullPage: true })

// Exactly one download button, and it is the only way out of here.
const buttons = await page.locator('.btn:has-text("Download")').allInnerTexts()
check('exactly one download button', buttons.length, 1)

const t2 = Date.now()
await page.locator('.btn:has-text("Download")').first().click()
// The app calls showSaveFilePicker, which the stub above points at the sink;
// "Saved ✓" only appears once the writable has been closed.
await page.waitForSelector('text=Saved ✓', { timeout: MINUTES * 60_000 })
const path = await sinkPath
log(`written in ${((Date.now() - t2) / 1000).toFixed(1)} s — ${(statSync(path).size / 1e6).toFixed(1)} MB`)
await b.close()
log(errs.length ? 'CONSOLE ERRORS: ' + errs.slice(0, 4).join(' | ') : 'no console errors')

failed += await checkCollection(path)

console.log(failed ? `\n${failed} check(s) FAILED\n` : '\nEverything checked out\n')
process.exit(failed ? 1 : 0)
