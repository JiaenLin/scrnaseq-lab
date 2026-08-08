// Drive the lab in a real browser: open a file, convert, save the ONE result.
// Usage: node lab-run.mjs <h5ad path> <out dir> <label>
import { chromium } from 'playwright-core'
import { mkdirSync } from 'node:fs'

const [input, outDir, label] = process.argv.slice(2)
mkdirSync(outDir, { recursive: true })
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)

const b = await chromium.launch({
  executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe',
})
const p = await b.newPage({ viewport: { width: 1280, height: 1600 }, acceptDownloads: true })
p.on('pageerror', e => log('[pageerror]', e.message.slice(0, 300)))
p.on('console', m => { if (m.type() === 'error') log('[console.error]', m.text().slice(0, 200)) })

// Force the anchor fallback: Playwright cannot drive a native save picker.
// Everything before the write is identical either way.
await p.addInitScript(() => { delete window.showSaveFilePicker })

await p.goto('http://localhost:' + (process.env.PORT || 4244) + '/', { waitUntil: 'networkidle' })
await p.setInputFiles('input[type=file]', input)
try {
  await p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 300000 })
} catch (e) {
  // A missed selector and a real hang look identical unless the page is dumped.
  console.log('---- SCAN DID NOT REACH THE REVIEW SCREEN. PAGE TEXT: ----')
  console.log(await p.locator('body').innerText())
  console.log('---------------------------------------------------------')
  await p.screenshot({ path: `${outDir}/${label}-scanfail.png`, fullPage: true })
  await b.close()
  process.exit(1)
}
log('scanned')

// AUDIT: is there any split control on the review screen?
const reviewText = await p.locator('body').innerText()
const splitHits = reviewText.split('\n').filter(l => /split|shard|part\b|parts\b|piece/i.test(l))
log('review-screen lines mentioning split/part:', JSON.stringify(splitHits))
const controls = await p.locator('button, select, input').evaluateAll(
  els => els.map(e => `${e.tagName}[${e.type || ''}]:${(e.innerText || e.value || e.name || '').slice(0, 60)}`))
log('review controls:', JSON.stringify(controls, null, 1))

const convertBtn = await p.locator('button.btn-primary').last().innerText()
log('convert button text:', JSON.stringify(convertBtn))
const t0 = Date.now()
await p.locator('button.btn-primary').last().click()

const done = p.waitForSelector('text=Every decision the converter made', { timeout: 3600000 })
  .then(() => 'DONE').catch(e => 'FAILED: ' + e.message.slice(0, 200))
const poll = (async () => {
  let last = ''
  for (;;) {
    await p.waitForTimeout(10000)
    const btn = await p.locator('button[disabled]').first().innerText().catch(() => '')
    if (btn && btn !== last) { log('  progress:', btn.replace(/\n/g, ' ').slice(0, 90)); last = btn }
    const err = await p.locator('.note').first().innerText().catch(() => '')
    if (/did not convert/i.test(err)) return 'ERROR: ' + err.slice(0, 300)
  }
})()
const r = await Promise.race([done, poll])
log(`convert -> ${r} in ${((Date.now() - t0) / 1000).toFixed(0)}s`)
if (r !== 'DONE') { await p.screenshot({ path: `${outDir}/${label}-fail.png` }); await b.close(); process.exit(1) }

// AUDIT the finish screen.
const finish = await p.locator('body').innerText()
console.log('---- FINISH SCREEN TEXT ----\n' + finish + '\n----------------------------')
const dlButtons = await p.locator('button').evaluateAll(
  els => els.map(e => e.innerText.replace(/\n/g, ' ')).filter(Boolean))
log('all buttons on finish screen:', JSON.stringify(dlButtons))
log('tables on finish screen:', await p.locator('table').count())

const dl = p.waitForEvent('download', { timeout: 1800000 })
await p.locator('button:has-text("Download")').first().click()
const d = await dl
const dest = `${outDir}/${d.suggestedFilename()}`
await d.saveAs(dest)
log('saved ->', dest)
await p.screenshot({ path: `${outDir}/${label}-done.png`, fullPage: true })
await b.close()
