// Developer script — NOT part of `npm test`.
//
// Opens the 9.3 GB atlas in a real browser and reports what the app does with
// it. Needs a `vite preview` on :4174, a Playwright chromium, and the file.
import { chromium } from 'playwright-core'

const EXE = 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe'
const OUT = 'C:/Users/Lin/AppData/Local/Temp/labshots'
const FILE = process.argv[2] ?? 'E:/developing_mouse_nervous_system.h5ad'
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)

const b = await chromium.launch({ executablePath: EXE })
const page = await b.newPage({ viewport: { width: 1200, height: 1500 } })
const errs = []
page.on('pageerror', e => errs.push('PAGEERROR ' + e.message))
page.on('console', m => { if (m.type() === 'error') errs.push('CONSOLE ' + m.text()) })

await page.goto('http://localhost:4174/', { waitUntil: 'networkidle' })
log('handing over', FILE)
const t0 = Date.now()
await page.setInputFiles('input[type=file]', FILE)

// Wait for whichever comes first: the questions, or a refusal.
const ok = page.waitForSelector('text=Which column holds the cell type annotation?', { timeout: 600000 })
const bad = page.waitForSelector('.note:not(.note-info) b', { timeout: 600000 })
const first = await Promise.race([ok.then(() => 'ok'), bad.then(() => 'error')]).catch(e => 'timeout: ' + e.message)
const secs = ((Date.now() - t0) / 1000).toFixed(1)

if (first === 'ok') {
  log(`SCANNED IN THE BROWSER in ${secs} s`)
  log((await page.locator('.sub').first().innerText()).replace(/\n/g, ' '))
  const split = await page.locator('text=Split it into several').count()
  log('split card shown:', split > 0)
  for (const r of (await page.locator('button:has(.badge-file)').allInnerTexts()).slice(0, 6)) {
    log('  |', r.replace(/\n/g, ' · ').slice(0, 120))
  }
  await page.screenshot({ path: `${OUT}/atlas-review.png`, fullPage: false })
} else if (first === 'error') {
  log(`refused after ${secs} s:`, (await page.locator('.note:not(.note-info)').first().innerText()).replace(/\n/g, ' '))
} else {
  log(first, `after ${secs} s`)
}
await page.close()
await b.close()
log(errs.length ? 'ERRORS: ' + errs.slice(0, 4).join(' | ') : 'no console errors')
