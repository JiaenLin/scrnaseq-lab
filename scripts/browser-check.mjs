// Developer script — NOT part of `npm test`.
//
// Drives the built app in a real browser: drops each object on the file input,
// waits for the questions, converts, downloads, and parses the download with the
// studio's reader. Needs a `vite preview` on :4174, a Playwright chromium, the
// sibling studio repo and the test objects.
import { chromium } from 'playwright-core'
import { existsSync, mkdirSync } from 'node:fs'
const EXE = 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe'
const OUT = 'C:/Users/Lin/AppData/Local/Temp/labshots'
if (!existsSync(OUT)) mkdirSync(OUT, { recursive: true })
const T = 'C:/Users/Lin/scrnaseq-studio/testdata'

const b = await chromium.launch({ executablePath: EXE })
const errors = []
const run = async (file, tag) => {
  const page = await b.newPage({ viewport: { width: 1280, height: 1600 } })
  page.on('console', m => { if (m.type() === 'error') errors.push(`[${tag}] ${m.text()}`) })
  page.on('pageerror', e => errors.push(`[${tag}] PAGEERROR ${e.message}`))
  await page.goto('http://localhost:4174/', { waitUntil: 'networkidle' })
  const t0 = Date.now()
  await page.setInputFiles('input[type=file]', `${T}/${file}`)
  await page.waitForSelector('text=Which column holds the cell type annotation?', { timeout: 300000 })
  console.log(`${tag}: scanned in ${Date.now() - t0} ms`)
  await page.screenshot({ path: `${OUT}/${tag}-review.png`, fullPage: true })
  const summary = await page.locator('.sub').first().innerText()
  console.log(`  ${summary.replace(/\n/g,' ')}`)
  // the capability table
  const rows = await page.locator('table.t').last().locator('tr').allInnerTexts()
  rows.slice(1).forEach(r => console.log('  |', r.replace(/\t/g, ' · ')))
  const t1 = Date.now()
  await page.click('button:has-text("Convert to bundle.zip")')
  await page.waitForSelector('text=Every decision the converter made', { timeout: 600000 })
  console.log(`  converted in ${Date.now() - t1} ms`)
  const dl = page.waitForEvent('download', { timeout: 60000 })
  await page.click('button:has-text("Download")')
  const d = await dl
  const path = `${OUT}/${tag}.zip`
  await d.saveAs(path)
  console.log(`  downloaded ${d.suggestedFilename()}`)
  await page.screenshot({ path: `${OUT}/${tag}-done.png`, fullPage: true })
  await page.close()
  return path
}

const a = await run('pbmc3k_final.rds', 'rds')
const c = await run('pbmc3k_processed.h5ad', 'h5ad')
await b.close()

// The downloaded zips must open in the studio's reader.
const { readFileSync } = await import('node:fs')
const { parseBundle } = await import('../../scrnaseq-studio/src/lib/bundle.ts')
for (const [p, tag] of [[a, 'rds'], [c, 'h5ad']]) {
  const buf = readFileSync(p)
  const bun = parseBundle(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength))
  console.log(`STUDIO READS ${tag}: ${bun.meta.nCells} cells · ${bun.genes.length} genes · ${bun.meta.clusters.length} clusters · pb ${bun.pseudobulk?.columns.length ?? 0}`)
}
console.log(errors.length ? `\nCONSOLE ERRORS:\n${errors.join('\n')}` : '\nNo console errors.')
