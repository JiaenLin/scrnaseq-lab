// Developer script — NOT part of `npm test`.
//
// The two format additions as a person sees them: open an object with two
// embeddings and accession-named genes, and check the page says so before the
// conversion and that the conversion's own notes say so afterwards.
//
//   npm run build && npx vite preview --port 4271 --strictPort &
//   node scripts/browser-names.mjs
//
// Every assertion is a string grepped out of src/, and a miss prints the page
// text — a silent false negative here looks exactly like a hang.

import { chromium } from 'playwright-core'

const PORT = process.argv[3] ?? 4271
const FIXTURE = process.argv[2] ?? 'C:/Users/Lin/AppData/Local/Temp/loop/accession_fixture.h5ad'
const CHROME = 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe'
const SHOT = 'C:/Users/Lin/AppData/Local/Temp/loop'

let failed = 0
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), ...a)

const b = await chromium.launch({ executablePath: CHROME })
const p = await b.newPage({ viewport: { width: 1280, height: 1600 } })
p.on('pageerror', e => { failed++; log('[pageerror]', e.message.slice(0, 300)) })
p.on('console', m => { if (m.type() === 'error') log('[console]', m.text().slice(0, 200)) })

async function says(what, fragment) {
  const text = await p.locator('body').innerText()
  const ok = text.includes(fragment)
  if (!ok) {
    failed++
    console.log(`  FAIL ${what}\n        looked for: ${fragment}\n        PAGE TEXT:\n${text.replace(/^/gm, '        | ')}`)
  } else {
    console.log(`  ok   ${what}`)
  }
}

await p.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' })
log('opening', FIXTURE)
await p.setInputFiles('input[type=file]', FIXTURE)
await p.waitForSelector('text=Which column is the cell type annotation?', { timeout: 120000 })

console.log('\nTHE REVIEW PAGE SAYS WHAT WILL TRAVEL')
await says('it counts both embeddings', '2 embeddings')
await says('the picker no longer reads as a choice between them',
  'Which one the studio opens on. All of them travel in the file.')
await says('and says so again under it', 'All 2 travel in the file; the studio can switch between them.')
await says('both are offered', 'X_tSNE')
await says('the gene naming is stated before the conversion, not after',
  'Genes are named by accessions (ENSMUSG00000000000); the symbols in var/Gene travel with them')
await p.screenshot({ path: `${SHOT}/lab-names-review.png`, fullPage: true })

console.log('\nAND THE CONVERSION SAYS WHAT IT DID')
await p.click('button:has-text("Convert")')
await p.waitForSelector('text=Every decision the converter made', { timeout: 300000 })
await says('every embedding is carried', '2 embeddings carried (X_UMAP, X_tSNE); X_UMAP opens by default')
await says('the naming is recorded', 'gene names are accessions; var/Gene carries the symbols and travels with them')
await says('genes with no symbol are counted, not dropped',
  '1 of 40 genes have no symbol in the object; those keep the name the matrix is indexed by')
await says('and genes sharing a symbol are counted, not merged',
  '2 genes share a symbol with another gene; the rows are kept separate rather than summed')
await p.screenshot({ path: `${SHOT}/lab-names-done.png`, fullPage: true })

await b.close()
console.log(failed ? `\n${failed} check(s) FAILED\n` : '\nThe browser agrees\n')
process.exit(failed ? 1 : 0)
