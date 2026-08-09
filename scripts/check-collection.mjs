// Taking a finished collection apart, without loading it.
//
// Shared by the browser end-to-end script, and usable on its own against any
// container already on disk:
//
//   node scripts/check-collection.mjs <collection.zip>
//
// Everything reads through a file handle at the offsets it needs, which is the
// claim the format makes — verifying it any other way would prove nothing.
//
// Note what it does NOT use: node:fs openAsBlob reports `size` modulo 2^32, so
// a 5.83 GB file comes back claiming 1.53 GB and every offset past 4 GB lands
// in the wrong place. A perfectly good container then reads as "not a
// collection", which cost a full atlas run to work out.

import { statSync } from 'node:fs'
import { open } from 'node:fs/promises'
import { unzipSync } from 'fflate'
import { parseBundle } from '../../scrnaseq-studio/src/lib/bundle.ts'
import { readCollectionIndex, readEntry } from '../src/lib/collection.ts'
import { makeChunkCache, readGenes } from '../src/lib/chunked.ts'

/**
 * Just enough Blob for readCollectionIndex: a size, and ranges read from the
 * file at the offset asked for and nowhere else.
 */
export async function fileBlob(p) {
  const fh = await open(p, 'r')
  const make = (from, to) => ({
    size: to - from,
    slice: (a = 0, b = to - from) => make(from + a, from + Math.min(b, to - from)),
    async arrayBuffer() {
      const buf = Buffer.alloc(to - from)
      let got = 0
      while (got < buf.length) {
        const r = await fh.read(buf, got, buf.length - got, from + got)
        if (r.bytesRead === 0) break
        got += r.bytesRead
      }
      return buf.buffer.slice(buf.byteOffset, buf.byteOffset + got)
    },
  })
  return make(0, statSync(p).size)
}

/** Where a part's bytes actually begin inside the container. */
async function partPayloadStart(file, start) {
  const h = new Uint8Array(await file.slice(start, start + 30).arrayBuffer())
  const d = new DataView(h.buffer, h.byteOffset, h.byteLength)
  return start + 30 + d.getUint16(26, true) + d.getUint16(28, true)
}

/** Byte offset of a STORED entry's payload within a zip held in memory, or -1. */
function storedOffset(zip, name) {
  const d = new DataView(zip.buffer, zip.byteOffset, zip.byteLength)
  let eocd = -1
  for (let i = zip.length - 22; i >= 0; i--) {
    if (d.getUint32(i, true) === 0x06054b50) { eocd = i; break }
  }
  if (eocd < 0) return -1
  const count = d.getUint16(eocd + 10, true)
  let p = d.getUint32(eocd + 16, true)
  for (let i = 0; i < count; i++) {
    const nameLen = d.getUint16(p + 28, true)
    const got = new TextDecoder().decode(zip.subarray(p + 46, p + 46 + nameLen))
    if (got === name && d.getUint16(p + 10, true) === 0) {
      const local = d.getUint32(p + 42, true)
      return local + 30 + d.getUint16(local + 26, true) + d.getUint16(local + 28, true)
    }
    p += 46 + nameLen + d.getUint16(p + 30, true) + d.getUint16(p + 32, true)
  }
  return -1
}

/** Check one container. Returns the number of failed checks. */
export async function checkCollection(path) {
  let failed = 0
  const check = (name, got, want) => {
    const ok = JSON.stringify(got) === JSON.stringify(want)
    if (!ok) failed++
    console.log(`  ${ok ? 'ok  ' : 'FAIL'} ${name}${ok ? '' : `\n         got  ${JSON.stringify(got)}\n         want ${JSON.stringify(want)}`}`)
  }

  console.log('\nTHE ONE FILE IS A COLLECTION')
  console.log(`  ${path} — ${(statSync(path).size / 1e9).toFixed(3)} GB`)
  const blob = await fileBlob(path)
  const idx = await readCollectionIndex(blob)
  if (!idx) {
    console.log('  FAIL it is not a collection at all')
    return failed + 1
  }
  const { meta, entries } = idx
  console.log(`  ${meta.parts.length} part(s) · splitBy=${meta.splitBy} · ${meta.nCells.toLocaleString()} cells · ${meta.nGenes.toLocaleString()} genes`)
  if (meta.reason) console.log(`  reason: ${meta.reason}`)
  check('every part named in the index is in the container',
    meta.parts.every(p => entries.has(p.file)), true)
  check('the parts add up to the collection',
    meta.parts.reduce((n, p) => n + p.nCells, 0), meta.nCells)
  check('nothing but the parts and the index is in there', entries.size, meta.parts.length + 1)

  // The level orders the parts cannot be trusted to imply. Optional, and older
  // containers have none — but where one is recorded it has to be the whole
  // object's list, or the studio reconstructs a sequence the experiment never
  // had and every menu reads wrong while every count stays right.
  const orders = [
    ['clusterOrder', meta.clusterOrder],
    ['condOrder', meta.condOrder],
    ...Object.entries(meta.extraOrder ?? {}).map(([k, v]) => [`extraOrder.${k}`, v]),
  ].filter(([, v]) => v)
  for (const [name, levels] of orders) {
    check(`${name} lists every level once`, new Set(levels).size, levels.length)
    console.log(`  ${name}: ${levels.length} levels — ${levels.slice(0, 6).join(', ')}${levels.length > 6 ? ' …' : ''}`)
  }
  if (!orders.length) console.log('  no level orders recorded — the studio will reconstruct them')

  // First, middle and last. The last two sit past the 4 GB mark in an atlas,
  // which is exactly where a 32-bit offset would have wrapped.
  const pick = meta.parts.length > 2
    ? [meta.parts[0], meta.parts[meta.parts.length >> 1], meta.parts[meta.parts.length - 1]]
    : meta.parts.slice(0, 2)

  for (const p of pick) {
    console.log(`\nPART "${p.key}" — ${p.nCells.toLocaleString()} cells, ${(p.bytes / 1e6).toFixed(1)} MB`)
    const entry = entries.get(p.file)
    const raw = await readEntry(blob, entry)
    check('the entry is the size the index claims', raw.length, p.bytes)

    const bundle = parseBundle(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength))
    check('the studio parses it', bundle.meta.nCells, p.nCells)
    check('it keeps the whole gene list', bundle.genes.length, meta.nGenes)
    console.log(`  clusters: ${bundle.meta.clusters.slice(0, 6).join(', ')}${bundle.meta.clusters.length > 6 ? ` ... +${bundle.meta.clusters.length - 6}` : ''}`)

    // Now the lazy path, byte ranges and all — exactly the studio's recipe.
    const f = unzipSync(raw)

    // The two format additions. Both are optional, so absence is not a failure
    // — but when they are there they must line up with the matrix, because a
    // misaligned name list mislabels every figure without ever throwing.
    const embs = bundle.meta.embeddings
    if (embs) {
      check('the default embedding is first and lives in embed.f32',
        [embs[0].key, embs[0].file], [bundle.meta.embedding, 'embed.f32'])
      for (const e of embs) {
        check(`${e.key} is in the part as ${e.file}`, !!f[e.file], true)
        check(`${e.key} has two coordinates for every cell`,
          new Float32Array(f[e.file].slice().buffer).length, 2 * p.nCells)
      }
      console.log(`  embeddings: ${embs.map(e => e.key).join(', ')}`)
    }
    // The columns beyond the three roles. Optional, so absence is not a
    // failure — but a part that kept a level it holds no cells for would give
    // that level a different meaning here from every other part, and the union
    // in the studio would then relabel cells rather than refuse.
    const extras = bundle.meta.extras ?? []
    for (const e of extras) {
      check(`${e.key} is in the part as ${e.file}`, !!f[e.file], true)
      const codes = new Uint16Array(f[e.file].slice().buffer)
      check(`${e.key} has one level index per cell`, codes.length, p.nCells)
      let outside = 0
      for (const c of codes) if (c >= e.levels.length) outside++
      check(`every ${e.key} code names a level this part defines`, outside, 0)
      check(`${e.key} carries no level this part has no cells for`,
        new Set(codes).size, e.levels.length)
    }
    console.log(extras.length
      ? `  extra columns: ${extras.map(e => `${e.key} (${e.levels.length}: ${e.levels.slice(0, 4).join(', ')}${e.levels.length > 4 ? ' …' : ''})`).join(' · ')}`
      : '  no extra columns in this part')

    const alias = bundle.meta.geneAlias
    if (alias) {
      const names = new TextDecoder().decode(f[alias.file]).split('\n')
      check('the alias has one name per matrix row', names.length, bundle.genes.length)
      check('and none of them is blank', names.every(n => n.length > 0), true)
      console.log(`  genes are ${bundle.meta.geneIdKind}s; ${alias.kind}s from var/${alias.column}: ${bundle.genes[0]} → ${names[0]} (${alias.missing} missing, ${alias.duplicated} shared)`)
    } else {
      console.log(bundle.meta.geneIdKind
        ? `  genes are ${bundle.meta.geneIdKind}s, and the object held no second naming`
        : '  written before gene naming was recorded — genes.txt is whatever the object had')
    }
    const chunkptr = new Int32Array(f['expr.chunkptr.i32'].slice().buffer)
    const inner = storedOffset(raw, 'expr.chunk.bin')
    check('expr.chunk.bin is stored inside the part', inner >= 0, true)
    const at = await partPayloadStart(blob, entry.start) + inner
    const cache = makeChunkCache(64 << 20)
    let fetched = 0
    const getBytes = async (from, to) => {
      fetched += to - from
      return new Uint8Array(await blob.slice(at + from, at + to).arrayBuffer())
    }

    // Genes spread across the list, so several chunks are touched.
    const n = bundle.genes.length
    const want = [0, 1, (n / 7) | 0, (n / 3) | 0, (n / 2) | 0, n - 2, n - 1]
      .filter((g, i, a) => g >= 0 && g < n && a.indexOf(g) === i)
    const tg = Date.now()
    const vecs = await readGenes(getBytes, chunkptr, bundle.indptr, bundle.meta.chunkGenes, want, cache)
    const same = want.every((g, i) => {
      const lo = bundle.indptr[g], hi = bundle.indptr[g + 1]
      const v = vecs[i]
      if (v.cells.length !== hi - lo) return false
      for (let k = 0; k < v.cells.length; k++) {
        if (v.cells[k] !== bundle.indices[lo + k] || v.values[k] !== bundle.data[lo + k]) return false
      }
      return true
    })
    check(`${want.length} genes read by byte range match the flat matrix`, same, true)
    console.log(`  read them in ${Date.now() - tg} ms, touching ${(fetched / 1024).toFixed(0)} KB `
      + `of a ${(p.bytes / 1e6).toFixed(1)} MB part (${cache.inflations} chunks inflated)`)
    const g = want[2]
    const v = vecs[want.indexOf(g)]
    console.log(`  ${bundle.genes[g]}: ${v.cells.length} cells expressing it, first value ${v.values[0]}`)
  }
  return failed
}

// Run directly on a file named on the command line.
if (process.argv[1]?.endsWith('check-collection.mjs')) {
  const bad = await checkCollection(process.argv[2])
  console.log(bad ? `\n${bad} check(s) FAILED\n` : '\nEverything checked out\n')
  process.exit(bad ? 1 : 0)
}
