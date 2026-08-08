import { chromium } from 'playwright-core'
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })
const p = await b.newPage()
p.on('pageerror', e => console.log('[pageerror]', e.message.slice(0, 200)))
await p.goto('http://localhost:4222/', { waitUntil: 'networkidle' })
const t0 = Date.now()
await p.setInputFiles('input[type=file]', process.argv[2])
const r = await Promise.race([
  // The heading the app actually renders.
  p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 300000 }).then(() => 'SCANNED'),
  p.waitForSelector('.note:not(.note-info) b', { timeout: 300000 }).then(() => 'REFUSED'),
]).catch(() => 'HANG')
console.log(`RESULT ${r} in ${((Date.now() - t0) / 1000).toFixed(1)} s`)
if (r === 'SCANNED') {
  console.log('  ' + (await p.locator('.sub').first().innerText()).replace(/\n/g, ' '))
  const split = await p.locator('text=Split it into several').count()
  console.log('  split card:', split > 0)
  for (const x of (await p.locator('button:has(.badge-file)').allInnerTexts()).slice(0, 5)) {
    console.log('   |', x.replace(/\n/g, ' · ').slice(0, 110))
  }
}
if (r === 'REFUSED') console.log('  ' + (await p.locator('.note:not(.note-info)').first().innerText()).replace(/\n/g, ' '))
await b.close()
