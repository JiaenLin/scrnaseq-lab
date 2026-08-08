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
import { decompressRds, RdsReader } from './rds.ts'
import { scanSeurat } from './seurat.ts'
import { buildBundle } from './build.ts'
import { shardScan, type Shard } from './shard.ts'
import type { Choices, Scan, ScanInfo } from './scan.ts'
import type { CellMajor } from './matrix.ts'
import { gunzipSync } from 'fflate'

let scan: Scan | null = null

/** Strip the closures: a Scan cannot be structured-cloned, its shape can. */
function describe(s: Scan, nnzPerCell: Int32Array | null): ScanInfo {
  return {
    format: s.format, source: s.source, nCells: s.nCells, columns: s.columns,
    matrices: s.matrices.map(m => ({
      key: m.key, geneSet: m.geneSet, nGenes: m.nGenes, nCells: m.nCells, kind: m.kind,
    })),
    embeddings: s.embeddings.map(e => ({ key: e.key, nDims: e.nDims })),
    geneSets: s.geneSets, provenance: s.provenance, notes: s.notes, nnzPerCell,
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
        post({ id, event: 'stage', stage: 'Decompressing the .rds' })
        const bytes = new Uint8Array(await file.arrayBuffer())
        const raw = decompressRds(bytes, b => gunzipSync(b))
        post({ id, event: 'stage', stage: 'Walking the R object' })
        const r = new RdsReader(raw)
        scan = scanSeurat(r, r.read(), file.name)
      } else if (name.endsWith('.h5ad') || name.endsWith('.h5')) {
        post({ id, event: 'stage', stage: 'Mounting the file' })
        // No copy: HDF5 reads only the bytes it asks for, straight from disk.
        const h5wasm = (await import('h5wasm')).default
        post({ id, event: 'stage', stage: 'Scanning obs, var and obsm' })
        scan = (await scanH5adFile(h5wasm as never, file)).scan
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

      // Whole object, unsplit.
      if (!shards || !passes) {
        post({ id, event: 'stage', stage: 'Reading the expression matrix' })
        const res = buildBundle(parent, choices, s => post({ id, event: 'stage', stage: s }))
        post({ id, event: 'bundle', key: null, name: `${safe(choices.label)}.zip`,
          zip: res.zip, meta: res.meta, pseudobulkColumns: res.pseudobulkColumns },
          [res.zip.buffer])
        post({ id, event: 'done' })
        return
      }

      const src = parent.matrices.find(m => m.loadGroups)
      const total = passes.reduce((n, p) => n + p.length, 0)
      let made = 0
      for (let pi = 0; pi < passes.length; pi++) {
        const group = passes[pi].map(i => shards[i])
        post({ id, event: 'stage',
          stage: `Pass ${pi + 1} of ${passes.length} — reading ${group.length} shard${group.length > 1 ? 's' : ''}` })

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
          post({ id, event: 'stage', stage: `Writing ${s.key} (${made + 1} of ${total})` })
          const sub = shardScan(parent, { ...s, cells: Int32Array.from(s.cells) },
            new Map([[src.key, mats[gi][0]]]), `${parent.source} · ${choices.label}`)
          const res = buildBundle(sub, { ...choices, label: `${choices.label} ${s.key}` })
          made++
          post({ id, event: 'bundle', key: s.key,
            name: `${safe(choices.label)}__${safe(s.key)}.zip`,
            zip: res.zip, meta: res.meta, pseudobulkColumns: res.pseudobulkColumns },
            [res.zip.buffer])
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

const safe = (s: string) => (s || 'bundle').replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 80)
