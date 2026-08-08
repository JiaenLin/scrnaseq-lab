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
    w.addEventListener('message', onMessage)
    w.postMessage({ ...msg, id })
  })
}

/** Open and scan a file. The worker keeps it open afterwards. */
export const openFile = (file: File, h: Handlers = {}): Promise<ScanInfo> =>
  call<ScanInfo>({ cmd: 'open', file }, h, d => d.info as ScanInfo)

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
