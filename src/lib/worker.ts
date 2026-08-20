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
import { rdsCompression, RdsError, RdsReader } from './rds.ts'
import { Bytes, CHUNK } from './bytes.ts'
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
        post({ id, event: 'stage', stage: 'Decompressing the .rds' })
        const raw = await readRds(file, done => post({
          id, event: 'stage', stage: `Decompressing the .rds — ${(done / 1e9).toFixed(1)} GB`,
        }))
        post({ id, event: 'stage', stage: 'Walking the R object' })
        const r = new RdsReader(raw)
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
 * A .rds off disk as Bytes, decompressed on the way through.
 *
 * `DecompressionStream('gzip')` rather than fflate's gunzipSync: it is native,
 * it is incremental, and it never needs the whole compressed input or the whole
 * output in one place. The pieces it yields are repacked into CHUNK-sized
 * buffers because Bytes divides rather than searches to find one — see bytes.ts.
 *
 * An uncompressed .rds — saveRDS(compress = FALSE) — goes through the same
 * packer, so both produce the same shape and there is one path to be wrong in.
 */
async function readRds(file: File, onProgress: (bytes: number) => void): Promise<Bytes> {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer())
  // Two bytes, not three gigabytes, before refusing a format we cannot read.
  const how = rdsCompression(head[0], head[1])

  const parts: Uint8Array[] = []
  let cur = new Uint8Array(CHUNK)
  let at = 0
  let total = 0
  let told = 0
  const take = (b: Uint8Array) => {
    let off = 0
    while (off < b.length) {
      const n = Math.min(CHUNK - at, b.length - off)
      cur.set(b.subarray(off, off + n), at)
      at += n
      off += n
      if (at === CHUNK) {
        parts.push(cur)
        try {
          cur = new Uint8Array(CHUNK)
        } catch {
          /**
           * Out of room, and this is the ONE place that knows how much of the
           * object had arrived — so it is the only place that can say something
           * true about why.
           *
           * The old message was the raw "Array buffer allocation failed", which
           * is what a reader saw for a 2.95 GB .rds and says nothing about the
           * file, the size, or what to do instead. The number below is measured
           * rather than guessed: it is exactly how many bytes were decompressed
           * before the tab ran out.
           *
           * Seurat objects reach this size mostly through `scale.data`, which
           * is DENSE — genes x cells with no zeros omitted — and which no view
           * in the studio can use. Dropping it in R is usually the whole fix
           * and costs nothing anybody wanted.
           */
          const gb = (total / 1e9).toFixed(1)
          throw new RdsError(
            `this .rds is larger than the browser can hold. ${gb} GB of it had been `
            + 'decompressed when the tab ran out of memory, and an .rds has to be '
            + 'unpacked before it can be read — unlike an .h5ad, which is read in place '
            + 'a slice at a time however large it is.\n\n'
            + 'Most of a Seurat object this size is usually scale.data, which is dense '
            + 'and which the studio cannot use. In R:\n\n'
            + '    obj <- DietSeurat(obj, scale.data = FALSE, dimreducs = c("pca", "umap"))\n'
            + '    saveRDS(obj, "smaller.rds")\n\n'
            + 'Or convert to .h5ad instead, which has no size limit here:\n\n'
            + '    SeuratDisk::SaveH5Seurat(obj, "obj.h5Seurat")\n'
            + '    SeuratDisk::Convert("obj.h5Seurat", dest = "h5ad")')
        }
        at = 0
      }
    }
    total += b.length
    // Every 50 MB, not every chunk: a 512 MB piece is tens of seconds on a slow
    // disk and the stage line would sit still through all of it.
    if (total - told > 50e6) { told = total; onProgress(total) }
  }

  const src = how === 'gzip'
    ? file.stream().pipeThrough(new DecompressionStream('gzip'))
    : file.stream()
  const reader = src.getReader()
  for (;;) {
    const { done, value } = await reader.read()
    if (done) break
    if (value) take(value)
  }
  // Bytes requires every piece but the last to be exactly CHUNK.
  if (at) parts.push(cur.subarray(0, at))
  return new Bytes(parts)
}
