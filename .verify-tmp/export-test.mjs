// Click every export control on every tab and record what actually downloads.
import { chromium } from 'playwright-core'
import { mkdirSync, statSync } from 'node:fs'
const [file, label] = process.argv.slice(2)
const dir = `C:/Users/Lin/AppData/Local/Temp/verify/exports/${label}`
mkdirSync(dir, { recursive: true })
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), `[${label}]`, ...a)
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })
const p = await b.newPage({ viewport: { width: 1500, height: 1400 }, acceptDownloads: true })
p.on('pageerror', e => log('[pageerror]', e.message.slice(0, 160)))
await p.goto('http://localhost:4245/', { waitUntil: 'networkidle' })
await p.setInputFiles('input[type=file]', file)
await p.waitForSelector('button:has-text("Overview")', { timeout: 600000 })

for (const tab of ['Markers', 'DEG table', 'Volcano', 'Gene sets']) {
  await p.locator(`button:text-is("${tab}")`).first().click()
  await p.waitForTimeout(1500)
  for (let i = 0; i < 900; i++) {
    if (!(await p.locator('[role=status]').count())) break
    await p.waitForTimeout(1000)
  }
  await p.waitForTimeout(800)
  const btns = p.locator('button:has-text("⭳")')
  const n = await btns.count()
  for (let i = 0; i < n; i++) {
    const name = (await btns.nth(i).innerText()).trim()
    try {
      const dl = p.waitForEvent('download', { timeout: 90000 })
      await btns.nth(i).click()
      const d = await dl
      const dest = `${dir}/${tab.replace(/ /g, '_')}-${i}-${d.suggestedFilename()}`
      await d.saveAs(dest)
      log(`${tab} [${name}] -> ${d.suggestedFilename()} ${statSync(dest).size} bytes`)
    } catch (e) {
      log(`${tab} [${name}] -> NO DOWNLOAD (${String(e.message).slice(0, 80)})`)
    }
  }
}
await b.close()
