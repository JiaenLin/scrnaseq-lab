import { chromium } from 'playwright-core'
const OUT = 'C:/Users/Lin/AppData/Local/Temp/atlas-bundles'
import { mkdirSync } from 'node:fs'
mkdirSync(OUT, { recursive: true })
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })
const p = await b.newPage({ viewport: { width: 1200, height: 1400 } })
p.on('pageerror', e => log('[pageerror]', e.message.slice(0, 200)))
await p.goto('http://localhost:4222/', { waitUntil: 'networkidle' })
log('opening the atlas')
await p.setInputFiles('input[type=file]', 'E:/developing_mouse_nervous_system.h5ad')
await p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 300000 })
log('scanned; converting with the recommended split')
const t0 = Date.now()
await p.click('button:has-text("Convert to bundle.zip")')
// Report progress as it goes.
let last = ''
const done = p.waitForSelector('text=Every decision the converter made', { timeout: 3600000 })
    .then(() => 'DONE').catch(e => 'FAILED: ' + e.message.slice(0, 120))
const poll = (async () => {
  while (true) {
    await p.waitForTimeout(15000)
    const s = await p.locator('.sub, .note').first().innerText().catch(() => '')
    const btn = await p.locator('button:has-text("Converting"), button:has-text("Pass")').first().innerText().catch(() => '')
    const cur = (btn || s).replace(/\n/g, ' ').slice(0, 100)
    if (cur && cur !== last) { log('  ', cur); last = cur }
    const err = await p.locator('.note:not(.note-info)').count().catch(() => 0)
    if (err) { log('   ERROR NOTE:', (await p.locator('.note:not(.note-info)').first().innerText()).slice(0, 200)); break }
  }
})()
const r = await Promise.race([done, poll])
log(`convert -> ${r} in ${((Date.now() - t0) / 60000).toFixed(1)} min`)
if (r === 'DONE') {
  const rows = await p.locator('table.t tbody tr').allInnerTexts()
  log(`${rows.length} bundles produced. First 5:`)
  rows.slice(0, 5).forEach(x => log('   |', x.replace(/\t/g, ' · ')))
  const dl = p.waitForEvent('download', { timeout: 120000 })
  await p.locator('table.t tbody tr').first().locator('button').click()
  const d = await dl
  await d.saveAs(`${OUT}/${d.suggestedFilename()}`)
  log('downloaded', d.suggestedFilename())
}
await p.screenshot({ path: 'C:/Users/Lin/AppData/Local/Temp/labshots/atlas-done.png', fullPage: false })
await b.close()
