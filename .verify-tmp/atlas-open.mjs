import { chromium } from 'playwright-core'
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })
const p = await b.newPage({ viewport: { width: 1500, height: 1200 } })
const errs = []
p.on('pageerror', e => { errs.push(e.message); log('[pageerror]', e.message.slice(0, 160)) })
await p.goto('http://localhost:4245/', { waitUntil: 'networkidle' })
log('dropping the 5.83 GB atlas collection')
const t0 = Date.now()
await p.setInputFiles('input[type=file]', 'C:/Users/Lin/AppData/Local/Temp/verify/atlas_collection.zip')
try {
  await p.waitForSelector('button:has-text("Overview")', { timeout: 240000 })
  log('tab bar appeared in', ((Date.now() - t0) / 1000).toFixed(1), 's')
} catch {
  log('NO TAB BAR after', ((Date.now() - t0) / 1000).toFixed(1), 's')
}
log('--- PAGE TEXT ---')
console.log((await p.locator('body').innerText()).slice(0, 1500))
log('--- END ---')
log('pageerrors:', JSON.stringify(errs.slice(0, 3)))
await p.screenshot({ path: 'C:/Users/Lin/AppData/Local/Temp/verify/shots/ATLAS-open.png', fullPage: false })
await b.close()
