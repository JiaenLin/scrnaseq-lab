import { chromium } from 'playwright-core'
import { mkdirSync } from 'node:fs'
const OUT = 'C:/Users/Lin/AppData/Local/Temp/loop'
mkdirSync(OUT, { recursive: true })
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })

// ---- LAB: one file out, no split control anywhere ------------------------
const p = await b.newPage({ viewport: { width: 1200, height: 1500 } })
const labErr = []
p.on('pageerror', e => labErr.push(e.message.slice(0, 160)))
p.on('console', m => { if (m.type() === 'error') labErr.push('console: ' + m.text().slice(0, 160)) })
await p.goto('http://localhost:4251/', { waitUntil: 'networkidle' })
// Force the Blob path so the save is automatable.
await p.evaluate(() => { delete window.showSaveFilePicker })
await p.setInputFiles('input[type=file]', 'C:/Users/Lin/scrnaseq-studio/testdata/pbmc3k_processed.h5ad')
await p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 300000 })
const body = await p.locator('body').innerText()
log('lab scanned. split control present?',
  /split it into several|how to split|choose a split/i.test(body) ? 'YES — DEFECT' : 'no')
await p.click('button:has-text("Convert")')
const dl = p.waitForEvent('download', { timeout: 300000 })
await p.waitForSelector('text=Every decision the converter made', { timeout: 300000 }).catch(() => {})
const btn = p.locator('button:has-text("Download"), a:has-text("Download")').first()
if (await btn.count()) await btn.click()
const d = await dl
const file = `${OUT}/${d.suggestedFilename()}`
await d.saveAs(file)
log('lab produced ONE file:', d.suggestedFilename())
const done = await p.locator('body').innerText()
log('finish screen mentions parts?', /part/i.test(done) ? 'yes (as information)' : 'no')
await p.screenshot({ path: `${OUT}/lab-done.png`, fullPage: true })
await p.close()

// ---- STUDIO: opens it as one dataset, every tab works ---------------------
const q = await b.newPage({ viewport: { width: 1400, height: 1200 } })
const stErr = []
q.on('pageerror', e => stErr.push(e.message.slice(0, 160)))
q.on('console', m => { if (m.type() === 'error') stErr.push('console: ' + m.text().slice(0, 160)) })
await q.goto('http://localhost:4252/', { waitUntil: 'networkidle' })
await q.setInputFiles('input[type=file]', file)
await q.waitForSelector('text=Overview', { timeout: 300000 })
log('studio opened it')
const sbody = await q.locator('body').innerText()
log('part switcher present?', /viewing part|switch part|part \d+ of/i.test(sbody) ? 'YES — DEFECT' : 'no')
const TABS = ['Overview', 'Cells', 'Composition', 'Markers', 'DEG table', 'Volcano',
              'Enrichment', 'Gene expression', 'Gene sets', 'Methods']
for (const t of TABS) {
  const before = stErr.length
  await q.click(`button:has-text("${t}")`).catch(() => {})
  await q.waitForTimeout(1800)
  const txt = (await q.locator('main, body').first().innerText()).replace(/\s+/g, ' ')
  log(`  ${t.padEnd(16)} ${stErr.length > before ? 'ERROR: ' + stErr.slice(before).join(' | ') : 'ok'} · ${txt.slice(0, 70)}`)
}
await q.screenshot({ path: `${OUT}/studio-final.png`, fullPage: false })
await q.close()
await b.close()
log('lab console errors:', labErr.length ? labErr.slice(0, 3).join(' | ') : 'none')
log('studio console errors:', stErr.length ? stErr.slice(0, 3).join(' | ') : 'none')
