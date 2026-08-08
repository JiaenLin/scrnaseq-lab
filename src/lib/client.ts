// Talking to the reader worker.
//
// One worker for the session, because it holds the open file: re-opening a
// 9.3 GB atlas for every question would undo the point of mounting it.

import type { Choices, ScanInfo } from './scan.ts'
import type { Shard } from './shard.ts'
import type { BundleMeta } from './build.ts'

export interface BuiltBundle {
  /** The split level this came from, or null when the object was not split. */
  key: string | null
  name: string
  blob: Blob
  meta: BundleMeta
  pseudobulkColumns: number
}

interface Handlers {
  onStage?: (s: string) => void
  onProgress?: (pass: number, frac: number) => void
  onBundle?: (b: BuiltBundle) => void
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
  msg: Record<string, unknown>, h: Handlers, finish: (m: Record<string, unknown>) => T,
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
          h.onBundle?.({
            key: d.key, name: d.name,
            blob: new Blob([d.zip], { type: 'application/zip' }),
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

/**
 * Where the worker should load the HDF5 reader from.
 *
 * Resolved here rather than in the worker: the worker's own URL is
 * /src/lib/worker.ts in dev and /assets/worker-<hash>.js in a build, so a path
 * relative to it lands somewhere different each time. `document.baseURI`
 * follows Vite's base, so this is also correct under a GitHub Pages subpath.
 */
const READER_URL = new URL('vendor/h5wasm/hdf5_hl.js', document.baseURI).href

/** Open and scan a file. The worker keeps it open afterwards. */
export const openFile = (file: File, h: Handlers = {}): Promise<ScanInfo> =>
  call<ScanInfo>({ cmd: 'open', file, readerUrl: READER_URL }, h, d => d.info as ScanInfo)

/**
 * Convert, either whole or split into shards.
 *
 * `passes` groups shard indices into forward traversals of the file, so the
 * caller decides how much memory a pass may use and the worker just follows it.
 */
export const convert = (
  choices: Choices, split: { shards: Shard[]; passes: number[][] } | null, h: Handlers = {},
): Promise<void> =>
  call<void>({
    cmd: 'convert', choices,
    shards: split ? split.shards : null,
    passes: split ? split.passes : null,
  }, h, () => undefined)
