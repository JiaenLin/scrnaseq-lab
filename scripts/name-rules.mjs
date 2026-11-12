// The source read of the facts walk: every rule in the lab's source that could
// decide a fact from a name. A method's, a tool's or a normalisation's word in
// any of the forms a test over a name takes:
//   - a regex literal, wherever it stands (after `{`, `(`, `=`, `return`, ...)
//     and however it is applied (`.test(`, `.exec(`, `.match(`, `.search(`,
//     `.replace(`, `.split(`): its body is read, not what follows it;
//   - a string given to `.includes(`, `.startsWith(`, `.endsWith(`,
//     `.indexOf(`, `.lastIndexOf(`, `.match(`, `.search(` or `RegExp(`;
//   - a string compared with `===`, `!==`, `==` or `!=`, or a `case` label;
//   - a list of names: an array of string literals one of which holds the word
//     (what `.includes(x)`, `.indexOf(x)` or `.some(m => x === m)` tests against).
// Run alone, `node scripts/name-rules.mjs <src dir>` prints each place, one a line.

import { readFileSync, readdirSync, lstatSync } from 'node:fs'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'

export const WORD = /harmony|scanvi|scvi|mnn|bbknn|combat|scanorama|liger|integrat|doublet|scrublet|soupx|cellbender|decontx|ambient|pca|log1p|cp10k|normali/i

const STR = String.raw`'((?:\\.|[^'\\\n])*)'|"((?:\\.|[^"\\\n])*)"|` + '`((?:\\\\.|[^`\\\\\\n$])*)`'
const FORMS = [
  // a regex literal: a slash where an expression may start, a body, a slash, flags
  new RegExp(String.raw`(?:^|[=(,:!&|?{}\[;+\-*%<>~^]|\breturn|\btypeof|\bcase)\s*\/((?:\\.|\[(?:\\.|[^\]\\\n])*\]|[^/\\\n\[])+)\/[dgimsuyv]*`, 'g'),
  new RegExp(String.raw`\.(?:includes|startsWith|endsWith|indexOf|lastIndexOf|match|matchAll|search)\(\s*(?:${STR})`, 'g'),
  new RegExp(String.raw`\bRegExp\(\s*(?:${STR})`, 'g'),
  new RegExp(String.raw`[!=]==?\s*(?:${STR})`, 'g'),
  new RegExp(String.raw`(?:${STR})\s*[!=]==?`, 'g'),
  new RegExp(String.raw`\bcase\s+(?:${STR})\s*:`, 'g'),
]
// A list of names: `[` then only string literals, commas and spaces, then `]`.
const LIST = new RegExp(String.raw`\[\s*((?:(?:${STR})\s*,?\s*)+)\]`, 'g')
const ITEM = new RegExp(STR, 'g')

/** The places in one line of code (comments already removed) that test a name for a word. */
export function rulesIn(code) {
  const out = []
  for (const re of FORMS) {
    for (const m of code.matchAll(re)) {
      const body = m.slice(1).find(x => x !== undefined) ?? ''
      if (WORD.test(body)) out.push(body)
    }
  }
  for (const m of code.matchAll(LIST)) {
    const items = [...m[1].matchAll(ITEM)].map(x => x.slice(1).find(y => y !== undefined) ?? '')
    if (items.some(x => WORD.test(x))) out.push(`[${items.join(', ')}]`)
  }
  return out
}

/** Every such place under `dir`, as `file:line: what`. Symlinks are not followed. */
export function nameRules(dir) {
  const files = []
  const walk = (d) => {
    for (const n of readdirSync(d)) {
      const p = join(d, n)
      const st = lstatSync(p)
      if (st.isSymbolicLink()) continue
      if (st.isDirectory()) walk(p)
      else if (/\.tsx?$/.test(n)) files.push(p)
    }
  }
  walk(dir)
  const hits = []
  for (const f of files.sort()) {
    let block = false
    readFileSync(f, 'utf8').split('\n').forEach((line, i) => {
      let code = line
      if (block) { const end = code.indexOf('*/'); if (end < 0) return; code = code.slice(end + 2); block = false }
      code = code.replace(/\/\*.*?\*\//g, '')
      const open = code.indexOf('/*')
      if (open >= 0 && !/['"`]/.test(code.slice(0, open))) { code = code.slice(0, open); block = true }
      code = code.replace(/(^|[^:\\'"`])\/\/.*$/, '$1').trim()
      if (!code || code.startsWith('*')) return
      for (const s of rulesIn(code)) hits.push(`${f.slice(dir.length).replace(/^\/+/, '')}:${i + 1}: ${s}`)
    })
  }
  return hits
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const hits = nameRules(process.argv[2] ?? fileURLToPath(new URL('../src', import.meta.url)))
  for (const h of hits) console.log(h)
  console.log(`${hits.length} place${hits.length === 1 ? '' : 's'}`)
}
