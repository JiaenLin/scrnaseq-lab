import { chromium } from 'playwright-core'
const b = await chromium.launch({ executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe' })
const p = await b.newPage({ viewport: { width: 1200, height: 1600 } })
await p.goto('http://localhost:4222/', { waitUntil: 'networkidle' })
await p.setInputFiles('input[type=file]', 'E:/developing_mouse_nervous_system.h5ad')
await p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 300000 })
// The chosen card is the accent-bordered one in Step 1.
const rows = await p.locator('button:has(.badge-file)').allInnerTexts()
console.log('cluster candidates, best first:')
rows.slice(0, 5).forEach(r => console.log('  |', r.replace(/\n/g, ' · ').slice(0, 120)))
const cap = await p.locator('table.t').last().innerText()
console.log('capability table:'); cap.split('\n').slice(0, 12).forEach(l => console.log('   ', l.replace(/\t/g, ' · ')))
await p.screenshot({ path: 'C:/Users/Lin/AppData/Local/Temp/labshots/atlas-review.png', fullPage: false })
await b.close()
