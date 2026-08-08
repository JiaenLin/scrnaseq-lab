// Talking to the reader worker.
//
// One worker for the session, because it holds the open file: re-opening a
// 9.3 GB atlas for every question would undo the point of mounting it.
//
// This is also where the parts become one file. The worker writes bundles and
// hands each one over; nothing downstream of here knows there was more than
// one, because `convert` resolves with a single collection blob.

import type { Choices, ScanInfo } from './scan.ts'
import type { Shard } from './shard.ts'
import type { BundleMeta } from './build.ts'
import { COLLECTION_SCHEMA, INDEX_NAME, type CollectionMeta, type PartInfo } from './collection.ts'

/** One part, as it comes back from the worker. */
export interface BuiltPart {
  /** The split level this came from, or the label when the object was not split. */
  key: string
  /** Entry name inside the container. */
  file: string
  bytes: Uint8Array
  meta: BundleMeta
  pseudobulkColumns: number
}

/**
 * What the user downloads: one file, whatever happened inside.
 *
 * The parts are kept as they came, not composed into a Blob. A Blob is how this
 * used to work and it cannot carry an atlas — see collectionPieces — so the
 * container only ever exists as it is being written to disk. `bytes` is what
 * that file will weigh, counted rather than measured.
 */
export interface BuiltCollection {
  name: string
  meta: CollectionMeta
  parts: { file: string; bytes: Uint8Array }[]
  bytes: number
}

interface Handlers {
  onStage?: (s: string) => void
  onProgress?: (pass: number, frac: number) => void
  /** Fired per part, for progress only — the user is never shown the parts. */
  onPart?: (p: BuiltPart, total: number) => void
}

/** How the object is to be cut up. Decided by the lab, never asked. */
export interface SplitOrder {
  column: string
  reason: string
  shards: Shard[]
  passes: number[][]
}

let worker: Worker | null = null
let nextId = 1

function ensure(): Worker {
  worker ??= new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' })
  return worker
}

/** Drop the worker, and with it the open file handle and its memory. */
export function release(): void {
  worker?.terminate()
  worker = null
}

function call<T>(
  msg: Record<string, unknown>, h: Handlers, onPart: (p: Omit<BuiltPart, 'file'>) => void,
  finish: (m: Record<string, unknown>) => T,
): Promise<T> {
  const w = ensure()
  const id = nextId++
  return new Promise<T>((resolve, reject) => {
    const onMessage = (e: MessageEvent) => {
      const d = e.data
      // Reported by the worker's unhandledrejection hook; it belongs to
      // whatever call is in flight, so it is not filtered by id.
      if (d?.event === 'fatal') {
        w.removeEventListener('message', onMessage)
        reject(new Error(`${d.message}${d.stack ? ` — ${d.stack}` : ''}`))
        return
      }
      if (d?.id !== id) return
      switch (d.event) {
        case 'stage': h.onStage?.(d.stage); break
        case 'progress': h.onProgress?.(d.pass, d.frac); break
        case 'bundle':
          onPart({
            key: d.key, bytes: d.zip as Uint8Array,
            meta: d.meta, pseudobulkColumns: d.pseudobulkColumns,
          })
          break
        case 'done':
          w.removeEventListener('message', onMessage)
          resolve(finish(d))
          break
        case 'error':
          w.removeEventListener('message', onMessage)
          reject(new Error(d.message))
          break
      }
    }
    // A worker that dies while loading never answers, and without this the call
    // sits pending behind a spinner that means nothing. Any failure inside the
    // worker is a failure of this call.
    const onError = (ev: ErrorEvent) => {
      w.removeEventListener('message', onMessage)
      w.removeEventListener('error', onError)
      release()
      reject(new Error(`the reader failed to start: ${ev.message || 'no reason given'}`))
    }
    w.addEventListener('message', onMessage)
    w.addEventListener('error', onError)
    w.postMessage({ ...msg, id })
  })
}

/** Open and scan a file. The worker keeps it open afterwards. */
export const openFile = (file: File, h: Handlers = {}): Promise<ScanInfo> =>
  call<ScanInfo>({ cmd: 'open', file }, h, () => {}, d => d.info as ScanInfo)

const safe = (s: string) => (s || 'part').replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 80)

/**
 * Convert, and return the one file the user takes away.
 *
 * `split` is the lab's own decision (see shard.ts `chooseSplit`), not a user
 * setting. Whether it is null or names forty-three shards, the result has the
 * same shape: one collection.zip, so the studio has exactly one thing to open.
 */
export async function convert(
  choices: Choices, split: SplitOrder | null, h: Handlers = {},
): Promise<BuiltCollection> {
  const parts: BuiltPart[] = []
  const partInfo: PartInfo[] = []
  const used = new Set<string>()
  const expected = split ? split.shards.length : 1

  await call<void>({
    cmd: 'convert', choices,
    shards: split ? split.shards : null,
    passes: split ? split.passes : null,
  }, h, p => {
    // Distinct split levels can flatten onto the same safe name ("10x 3'" and
    // "10x 5'" both become "10x_"), and two entries with one name would make a
    // container whose index points twice at the same bytes.
    let file = `parts/${safe(p.key)}.zip`
    for (let i = 2; used.has(file); i++) file = `parts/${safe(p.key)}-${i}.zip`
    used.add(file)
    const part: BuiltPart = { ...p, file }
    parts.push(part)
    partInfo.push({
      key: p.key, file, nCells: p.meta.nCells, nnz: p.meta.nnz, bytes: p.bytes.length,
    })
    h.onPart?.(part, expected)
  }, () => undefined)

  if (parts.length === 0) throw new Error('the conversion produced nothing')

  const first = parts[0].meta
  const meta: CollectionMeta = {
    schema: COLLECTION_SCHEMA,
    label: choices.label || first.label,
    source: first.source,
    splitBy: split ? split.column : null,
    reason: split ? split.reason : null,
    nCells: parts.reduce((n, p) => n + p.meta.nCells, 0),
    nGenes: first.nGenes,
    parts: partInfo,
    // The conversion notes are decisions about the object as a whole — which
    // matrix, which normalization, what was missing — so they are identical in
    // every part and belong to the collection.
    notes: first.notes,
  }

  // What the container will weigh: every part, its local header, its central
  // directory record, and the end record. Counted here so the finish screen can
  // say the size without anything having been assembled.
  const bytes = partInfo.reduce(
    (n, p) => n + p.bytes + 2 * (30 + p.file.length) + 16, 0)
    + 2 * (30 + INDEX_NAME.length) + JSON.stringify(meta, null, 1).length + 22

  return {
    name: `${safe(choices.label || first.label)}.zip`,
    meta,
    parts: parts.map(p => ({ file: p.file, bytes: p.bytes })),
    bytes,
  }
}
