// Open a file in the studio and walk every tab, recording every control.
// Usage: node studio-walk.mjs <zip> <out.json> <label>
import { chromium } from 'playwright-core'
import { writeFileSync, mkdirSync } from 'node:fs'

const [file, out, label] = process.argv.slice(2)
mkdirSync('C:/Users/Lin/AppData/Local/Temp/verify/shots', { recursive: true })
const log = (...a) => console.log(new Date().toTimeString().slice(0, 8), `[${label}]`, ...a)

const b = await chromium.launch({
  executablePath: 'C:/Users/Lin/AppData/Local/ms-playwright/chromium-1228/chrome-win64/chrome.exe',
})
const p = await b.newPage({ viewport: { width: 1500, height: 1900 }, acceptDownloads: true })
const errors = []
p.on('pageerror', e => { errors.push(e.message.slice(0, 300)); log('[pageerror]', e.message.slice(0, 200)) })
p.on('console', m => { if (m.type() === 'error') { errors.push('console: ' + m.text().slice(0, 200)) } })

await p.goto('http://localhost:' + (process.env.SPORT || 4245) + '/', { waitUntil: 'networkidle' })
const t0 = Date.now()
await p.setInputFiles('input[type=file]', file)

// Wait for the tab bar to appear.
await p.waitForSelector('button:has-text("Overview")', { timeout: 900000 })
log('opened in', ((Date.now() - t0) / 1000).toFixed(1), 's')

const TABS = ['Overview', 'Cells', 'Composition', 'Markers', 'DEG table', 'Volcano',
  'Enrichment', 'Gene expression', 'Gene sets', 'Methods']

const snap = async () => {
  // Every interactive control: what it is, its label, and whether it is usable.
  return await p.evaluate(() => {
    const vis = el => {
      const r = el.getBoundingClientRect()
      return r.width > 0 && r.height > 0
    }
    const ctl = []
    for (const el of document.querySelectorAll('button, select, input, textarea, a[download], a.btn')) {
      if (!vis(el)) continue
      const t = el.tagName.toLowerCase()
      let name = (el.innerText || el.value || el.placeholder || el.getAttribute('aria-label') || el.name || '').toString().replace(/\s+/g, ' ').trim().slice(0, 70)
      const o = { t, name, disabled: !!el.disabled }
      if (t === 'select') o.options = [...el.options].map(x => x.text.slice(0, 40))
      if (t === 'input') o.type = el.type
      ctl.push(o)
    }
    return {
      controls: ctl,
      tables: document.querySelectorAll('table').length,
      rows: document.querySelectorAll('table tbody tr').length,
      svgs: document.querySelectorAll('svg').length,
      canvas: document.querySelectorAll('canvas').length,
      text: document.body.innerText,
    }
  })
}

const result = { label, file, tabs: {}, errors: [] }

for (const tab of TABS) {
  const btn = p.locator(`button:text-is("${tab}")`)
  if (await btn.count() === 0) {
    log(`TAB MISSING: ${tab}`)
    result.tabs[tab] = { present: false }
    continue
  }
  const disabled = await btn.first().isDisabled()
  await btn.first().click()
  // Let any streamed computation settle.
  await p.waitForTimeout(800)
  // Progress.tsx renders role="status" for the whole of a streaming pass.
  const tStart = Date.now()
  for (let i = 0; i < 1200; i++) {
    const busy = await p.locator('[role=status]').count().catch(() => 0)
    if (!busy) break
    if (i % 20 === 0 && i) {
      const txt = await p.locator('[role=status]').first().innerText().catch(() => '')
      log(`   ...${tab} streaming: ${txt.replace(/\n/g, ' | ').slice(0, 90)}`)
    }
    await p.waitForTimeout(1000)
  }
  const waited = (Date.now() - tStart) / 1000
  await p.waitForTimeout(500)
  const s = await snap()
  result.tabs[tab] = { present: true, tabDisabled: disabled, seconds: waited, ...s }
  log(`${tab}: ${s.controls.length} controls (${s.controls.filter(c => c.disabled).length} disabled), ${s.tables} tables/${s.rows} rows, ${s.svgs} svg, ${waited.toFixed(1)}s`)
  await p.screenshot({ path: `C:/Users/Lin/AppData/Local/Temp/verify/shots/${label}-${tab.replace(/ /g, '_')}.png`, fullPage: false })
}

result.errors = errors
writeFileSync(out, JSON.stringify(result, null, 1))
log('wrote', out, '| pageerrors:', errors.length)
await b.close()
