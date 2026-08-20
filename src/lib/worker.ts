/// <reference lib="webworker" />
//
// The reader runs here, not on the page.
//
// Two reasons, and the first is the one that matters: WORKERFS — Emscripten's
// read-only filesystem over a File — serves reads with FileReaderSync, which
// exists only inside a worker. That is what lets a 9.3 GB atlas be opened
// without copying it anywhere. The second is that a conversion pass is a minute
// of solid CPU, and on the page that is a minute of frozen scrolling.
//
// The open HDF5 handle and its closures stay here for the session. Only
// descriptions cross back.

import { scanH5adFile } from './h5ad.ts'
import { rdsCompression, RdsReader } from './rds.ts'
import { Bytes, CHUNK } from './bytes.ts'
import { GzWindow } from './gzwindow.ts'
import { scanSeurat } from './seurat.ts'
import { buildBundle, streamableMatrix } from './build.ts'
import { shardScan, type Shard } from './shard.ts'
import type { Choices, Scan, ScanInfo } from './scan.ts'
import type { CellMajor } from './matrix.ts'

/* eslint-disable @typescript-eslint/no-explicit-any */
type H5Mod = any

/**
 * The HDF5 reader, loaded once per session.
 *
 * A long detour ended here: an earlier version vendored h5wasm to public/ and
 * imported it by URL, on the belief that Vite's worker pipeline broke it. It
 * does not. `ready` settles in about 15 ms inside the bundled module worker,
 * and the apparent hang was a stale selector in the browser probe looking for
 * a heading that had been renamed. The plain import is correct.
 */
let h5mod: H5Mod = null
async function h5(): Promise<H5Mod> {
  if (!h5mod) {
    const mod = await import('h5wasm')
    h5mod = (mod as { default?: H5Mod }).default ?? mod
    if (!h5mod?.File) throw new Error('the HDF5 reader did not load')
  }
  return h5mod
}

let scan: Scan | null = null

/** Strip the closures: a Scan cannot be structured-cloned, its shape can. */
function describe(s: Scan, nnzPerCell: Int32Array | null): ScanInfo {
  return {
    format: s.format, source: s.source, nCells: s.nCells, columns: s.columns,
    matrices: s.matrices.map(m => ({
      key: m.key, geneSet: m.geneSet, nGenes: m.nGenes, nCells: m.nCells, kind: m.kind,
    })),
    embeddings: s.embeddings.map(e => ({ key: e.key, nDims: e.nDims })),
    geneSets: s.geneSets, geneAliases: s.geneAliases,
    provenance: s.provenance, notes: s.notes, nnzPerCell,
  }
}

type Msg =
  | { id: number; cmd: 'open'; file: File }
  | { id: number; cmd: 'convert'; choices: Choices; shards: Shard[] | null; passes: number[][] | null }

const post = (m: unknown, transfer: Transferable[] = []) =>
  (self as unknown as Worker).postMessage(m, transfer)

self.onmessage = async (e: MessageEvent<Msg>) => {
  const { id } = e.data
  try {
    if (e.data.cmd === 'open') {
      const file = e.data.file
      const name = file.name.toLowerCase()
      if (name.endsWith('.rds')) {
        /**
         * Streamed, in pieces.
         *
         * This was `new Uint8Array(await file.arrayBuffer())` followed by
         * `gunzipSync` — two allocations, the compressed file and the whole
         * decompressed object, each as ONE buffer. V8 caps a single buffer at
         * 2.15 GB (measured by bisection in Chrome 151), so a 2.95 GB Seurat
         * .rds threw `RangeError: Array buffer allocation failed` on the first
         * line, before a byte was parsed — a file that would not open, with
         * nothing on the page to say why.
         *
         * The cap is per buffer and not a total: 6.14 GB held fine in 512 MB
         * pieces in the same tab. So the fix is to never ask for one big one.
         * DecompressionStream gunzips incrementally and the pieces go straight
         * into a Bytes, which is what RdsReader now walks. The compressed bytes
         * are never held at all — they are pulled from the File and dropped as
         * they are consumed.
         */
        post({ id, event: 'stage', stage: 'Reading the .rds' })
        const r = new RdsReader(openRds(file, (cz, out) => post({
          id,
          event: 'stage',
          stage: `Reading the .rds — ${(cz / 1e9).toFixed(1)} of `
            + `${(file.size / 1e9).toFixed(1)} GB in, ${(out / 1e9).toFixed(1)} GB unpacked`,
        })))
        scan = scanSeurat(r, r.read(), file.name)
      } else if (name.endsWith('.h5ad') || name.endsWith('.h5')) {
        post({ id, event: 'stage', stage: 'Loading the HDF5 reader' })
        const h5wasm = await h5()
        await h5wasm.ready
        // No copy: HDF5 reads only the bytes it asks for, straight from the File.
        post({ id, event: 'stage', stage: 'Scanning obs, var and obsm' })
        scan = (await scanH5adFile(h5wasm, file)).scan
      } else {
        throw new Error(`this is a ${name.split('.').pop() ?? 'file'} — the lab reads .h5ad (Scanpy/AnnData) and .rds (Seurat).`)
      }

      post({ id, event: 'stage', stage: 'Measuring the matrix' })
      const m = scan.matrices.find(x => x.nnzPerCell)
      const nnzPerCell = m?.nnzPerCell ? m.nnzPerCell() : null
      post({ id, event: 'done', info: describe(scan, nnzPerCell) })
      return
    }

    if (e.data.cmd === 'convert') {
      if (!scan) throw new Error('no file is open')
      const { choices, shards, passes } = e.data
      const parent = scan

      // Whole object, one part. Still a part: the container is written the same
      // way either way, so the studio has one shape to understand.
      if (!shards || !passes) {
        post({ id, event: 'stage', stage: 'Reading the expression matrix' })
        const res = buildBundle(parent, choices, s => post({ id, event: 'stage', stage: s }))
        post({ id, event: 'bundle', key: choices.label || parent.source,
          zip: res.zip, meta: res.meta, pseudobulkColumns: res.pseudobulkColumns },
          [res.zip.buffer])
        post({ id, event: 'done' })
        return
      }

      // Notes from the choice go onto the parent scan, so every part carries
      // them: buildBundle copies scan.notes, and the split path never calls
      // chooseMatrices on anything that can still see the rejected matrix.
      const src = streamableMatrix(parent, s => parent.notes.push(s))
      if (!src) {
        throw new Error(
          'this object is too large to convert in one piece, and none of its matrices can be read in parts — every one is either scaled or stored in a layout that cannot be streamed.')
      }
      const total = passes.reduce((n, p) => n + p.length, 0)
      let made = 0
      for (let pi = 0; pi < passes.length; pi++) {
        const group = passes[pi].map(i => shards[i])
        // Worded as progress, not as structure. The user never chose to split
        // and will never see the pieces, so the stage line says how far along
        // the read is, not how the object was carved up.
        post({ id, event: 'stage',
          stage: `Reading the expression matrix — pass ${pi + 1} of ${passes.length}` })

        // One forward traversal per pass; every shard in it is filled together.
        const mats: CellMajor[][] = []
        if (src?.loadGroups) {
          const got = src.loadGroups(group.map(s => Int32Array.from(s.cells)),
            frac => post({ id, event: 'progress', pass: pi, frac }))
          got.forEach(g => mats.push([g]))
        } else {
          throw new Error('this matrix cannot be read in subsets, so it cannot be split')
        }

        for (let gi = 0; gi < group.length; gi++) {
          const s = group[gi]
          post({ id, event: 'stage', stage: `Writing ${made + 1} of ${total}` })
          const sub = shardScan(parent, { ...s, cells: Int32Array.from(s.cells) },
            new Map([[src.key, mats[gi][0]]]), `${parent.source} · ${choices.label}`)
          const res = buildBundle(sub, { ...choices, label: `${choices.label} ${s.key}` })
          made++
          post({ id, event: 'bundle', key: s.key,
            zip: res.zip, meta: res.meta, pseudobulkColumns: res.pseudobulkColumns },
            [res.zip.buffer])
          // The shard's matrix is dead the moment its bundle is written, and a
          // pass can hold 120 M nonzeros — a gigabyte — across its shards. Drop
          // each one as it is used rather than keeping the whole pass alive
          // while the last shard compresses.
          mats[gi].length = 0
        }
      }
      post({ id, event: 'done' })
      return
    }
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err)
    post({ id, event: 'error', message: msg && msg !== 'undefined' ? msg : 'the reader failed without saying why' })
  }
}


/**
 * A .rds as a forward-only stream, decompressed by the parser as it walks.
 *
 * Nothing is held. The old shape here decompressed the whole object into memory
 * first — which is why a 16.26 GB Seurat file could not be opened on a machine
 * with 38.7 GB of RAM: a Chrome tab tops out around 15.6 GB whatever the
 * hardware. RdsReader never rewinds, so the bytes only have to arrive in order,
 * and GzWindow drops them behind the head. Peak memory becomes the slots that
 * are kept rather than the file.
 *
 * FileReaderSync because the parser is synchronous and cannot await a slice.
 * It exists only in a worker, which is where this runs.
 */
function openRds(file: File, onRead: (compressed: number, produced: number) => void) {
  const head = new Uint8Array(new FileReaderSync().readAsArrayBuffer(file.slice(0, 2)))
  // Two bytes, not sixteen gigabytes, before refusing a format we cannot read.
  const how = rdsCompression(head[0], head[1])
  const fr = new FileReaderSync()
  const slice = (a: number, b: number) =>
    new Uint8Array(fr.readAsArrayBuffer(file.slice(a, b)))

  if (how === 'gzip') return new GzWindow(file.size, slice, onRead)

  // saveRDS(compress = FALSE). Rare, and small by the time it is: an
  // uncompressed .rds this large would not have fitted on the disk it came
  // from. Read as chunked bytes, which has no ceiling worth worrying about.
  const parts: Uint8Array[] = []
  for (let at = 0; at < file.size; at += CHUNK) {
    parts.push(slice(at, Math.min(file.size, at + CHUNK)))
    onRead(Math.min(file.size, at + CHUNK), Math.min(file.size, at + CHUNK))
  }
  // Every piece but the last must be exactly CHUNK; the loop guarantees it.
  return new Bytes(parts)
}
