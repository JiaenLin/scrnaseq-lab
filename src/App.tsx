import { useState, type ReactNode } from 'react'
import Landing from './components/Landing.tsx'
import Review from './components/Review.tsx'
import { Card, Flow } from './components/Ui.tsx'
import {
  convert as runConvert, openFile, release, type BuiltCollection,
} from './lib/client.ts'
import {
  collectionPieces, writeCollection, type CollectionMeta,
} from './lib/collection.ts'
import { chooseSplit, type SplitDecision } from './lib/shard.ts'
import {
  candidates, guessCluster, guessEmbedding, guessRole,
  type Choices, type ScanInfo,
} from './lib/scan.ts'

const STUDIO = 'https://jiaenlin.github.io/scrnaseq-studio/'

/** One numbered step of the handoff out of this app. */
function Handoff({ n, title, done, children }: {
  n: number; title: string; done?: boolean; children: ReactNode
}) {
  return (
    <div className="flex gap-3 rounded-xl px-3.5 py-3" style={{ background: 'var(--sunk)' }}>
      <div
        className="grid h-[22px] w-[22px] flex-none place-items-center rounded-[7px] text-[11px] font-bold text-white"
        style={{ background: done ? 'var(--good)' : 'var(--accent)' }}
      >{done ? '✓' : n}</div>
      <div className="min-w-0 flex-1">
        <div className="mb-2 text-[13px] font-semibold">{title}</div>
        {children}
      </div>
    </div>
  )
}

/**
 * The other groupings this object carries — offered, never assumed.
 *
 * The three questions below this are roles: every view in the studio is per
 * cell type, across groups, within samples. This is not a fourth role. It is
 * the rest of what the object already knows about each cell — a dissection
 * beside an age, a Class above a Subclass — carried under its own name so the
 * studio can put any of it on either axis of a composition figure. Nothing is
 * chosen by default: a bundle with none of these is the bundle the lab has
 * always written, byte for byte.
 */
function Extras({ scan, choices, setChoices }: {
  scan: ScanInfo; choices: Choices; setChoices: (c: Choices) => void
}) {
  const usable = candidates(scan, 'extra',
    [choices.cluster, choices.sample, choices.condition])
  const picked = choices.extras ?? []
  // A column chosen as one of the three roles drops off the list above, and the
  // builder will not carry it twice — but if it were also ticked here it would
  // vanish while still ticked, with no way to untick it. So a stale pick stays
  // on screen until the user takes it off.
  const shown = [...usable, ...scan.columns.filter(c =>
    picked.includes(c.name) && !usable.some(u => u.name === c.name))]
  if (!shown.length) return null
  const toggle = (name: string) => setChoices({
    ...choices,
    extras: picked.includes(name) ? picked.filter(x => x !== name) : [...picked, name],
  })

  return (
    <div className="wrap pt-7">
      <Card
        eyebrow="Extra columns · optional"
        title="Anything else worth breaking a figure down by"
        sub={`Carried per cell under the object's own name. The studio then offers every
          combination of these with the cell type, the group and the sample — brain region
          by age, cell type by region. Pick none and nothing changes.`}
      >
        <div className="mt-3.5 grid gap-2"
          style={{ gridTemplateColumns: 'repeat(auto-fill,minmax(230px,1fr))' }}>
          {shown.map(c => {
            const on = picked.includes(c.name)
            return (
              <button key={c.name} onClick={() => toggle(c.name)}
                className="w-full rounded-xl px-3 py-2.5 text-left"
                style={{
                  border: `1.5px solid ${on ? 'var(--accent)' : 'var(--line)'}`,
                  background: on ? 'var(--accent-soft)' : 'var(--surface)',
                }}>
                <div className="flex items-baseline gap-2">
                  <span className="mono text-[12.5px] font-semibold">{c.name}</span>
                  <span className="badge badge-file">{c.levels.length} levels</span>
                </div>
                <div className="mt-1.5 text-[11.5px]" style={{ color: 'var(--ink-2)' }}>
                  {c.levels.slice(0, 3).join(' · ')}
                  {c.levels.length > 3 && ` · +${c.levels.length - 3} more`}
                </div>
              </button>
            )
          })}
        </div>
        {picked.length > 1 && (
          <p className="sub mt-2.5">
            {picked.length} extra columns. Every level of each one travels with every cell;
            the cost is 2 bytes per cell per column, against a matrix measured in gigabytes.
          </p>
        )}
      </Card>
    </div>
  )
}

/**
 * The level orders a split would otherwise have to guess.
 *
 * A part is written with the levels it uses and no others, so no part carries
 * the whole object's order — and the studio reconstructs what it can and falls
 * back to collation for the rest. Collation compares the digit run after the
 * dot as a number, which puts e16.5 before e16.25: no count changes, but every
 * menu, contrast picker and axis then offers a sequence the experiment never
 * had. The order is still here, in the scan, so it is written down.
 *
 * It lives here rather than in client.ts for the same reason `clusterOrder`
 * takes a parameter: the client talks to the worker and has never seen the
 * scan. These two want the same parameter list.
 */
function recordOrders(
  built: BuiltCollection, scan: ScanInfo, choices: Choices,
): BuiltCollection {
  const levelsOf = (name: string | null | undefined): string[] | undefined =>
    (name ? scan.columns.find(c => c.name === name)?.levels : undefined)

  const extraOrder: Record<string, string[]> = {}
  for (const name of choices.extras ?? []) {
    const levels = levelsOf(name)
    if (levels) extraOrder[name] = levels
  }
  const meta: CollectionMeta & {
    condOrder?: string[]
    extraOrder?: Record<string, string[]>
  } = {
    ...built.meta,
    condOrder: levelsOf(choices.condition),
    extraOrder: Object.keys(extraOrder).length ? extraOrder : undefined,
  }
  // The index is written from this object, so the byte total the save button
  // counts down has to follow it.
  const grew = JSON.stringify(meta, null, 1).length
    - JSON.stringify(built.meta, null, 1).length
  return { ...built, meta, bytes: built.bytes + grew }
}

/**
 * A sentence for the user, whatever was thrown.
 *
 * Emscripten's ErrnoError carries no `message`, so reading `.message` blindly
 * printed "That file did not open. undefined".
 */
const explain = (e: unknown): string => {
  const msg = e instanceof Error ? e.message : String(e)
  return msg && msg !== 'undefined'
    ? msg
    : 'the reader failed without saying why. If the object is very large, try splitting it.'
}

/**
 * The part of the File System Access API used here, which TypeScript's DOM
 * library does not declare.
 *
 * Written out rather than left to inference. `let handle` with no annotation is
 * an evolving `any`, so every member reached through it — `abort`, `name`,
 * `move` — typechecked whether or not it existed, on the one code path in this
 * app that cannot be exercised by a test. `move` is optional because it is
 * Chromium-only and the code has to cope with its absence.
 */
interface FileWriter {
  // Uint8Array, not BufferSource: a Uint8Array may be backed by a
  // SharedArrayBuffer, which BufferSource excludes, and these views come
  // straight off the generator that builds the container.
  write: (d: Uint8Array | Blob) => Promise<void>
  close: () => Promise<void>
  abort: () => Promise<void>
}
interface SaveHandle {
  readonly name: string
  createWritable: () => Promise<FileWriter>
  move?: (to: string) => Promise<void>
}
interface SavePicker {
  showSaveFilePicker?: (o: {
    suggestedName?: string
    types?: { description: string; accept: Record<string, string[]> }[]
  }) => Promise<SaveHandle>
}

/**
 * Put the one file on the user's disk.
 *
 * Two things rule out the obvious route of composing a Blob and pointing an
 * `<a download>` at a `blob:` URL, both measured in Chromium rather than
 * assumed. A Blob past a few hundred megabytes spills to disk and then throws
 * NotReadableError on every read, so the container cannot be built that way at
 * all. And a `blob:` download is cancelled somewhere between 256 MB and
 * 512 MB, while a 1 GB download over plain HTTP in the same browser saves fine.
 * An atlas collection is 5.8 GB.
 *
 * So the zip is generated piece by piece and written straight into the file the
 * user names, and nothing larger than one part is ever held twice. The Blob
 * route stays as the fallback for browsers without showSaveFilePicker, where it
 * is the only option and where the small objects it can carry are the ones
 * those users convert anyway.
 */
/**
 * Why a save failed, in terms of what to do about it.
 *
 * Chromium hands through the operating system's own sentence — "a file or
 * directory which could not be modified due to the state of the underlying
 * filesystem" — which names no cause and suggests no action. The causes worth
 * naming are few and all recoverable, and after a conversion that took a
 * quarter of an hour the one thing that must not happen is a dead end.
 */
/**
 * A save that failed, carrying how far it got.
 *
 * The distinction is the whole point. Chromium writes through a swap file and
 * renames it at the end, so a failure DURING the writes leaves a partial file
 * worth discarding, while a failure AT THE RENAME leaves a complete one worth
 * keeping. Told apart, the second is not a failure at all — it is a finished
 * file under the wrong name.
 */
/**
 * The extension the file is WRITTEN under, which is not the one it keeps.
 *
 * Chromium runs a Safe Browsing check over a finished file before renaming its
 * swap into place, and only for extensions its download protection treats as
 * risky — .zip among them, because an archive can carry an executable. On a
 * 5.80 GB container that check does not come back: close() rejects with
 * "Failed to perform Safe Browsing check" AFTER every byte is already on disk,
 * and the whole conversion is stranded in a .crswap file the page cannot reach.
 * Measured on this object, twice.
 *
 * A name its download protection does not claim is never checked, so the write
 * finishes. The file is still a zip byte for byte — the studio reads the
 * container, never the name — and `promote` below asks for the .zip name back
 * once the bytes are safely down, where failing costs nothing.
 */
const DATA_EXT = '.scstudio'

/**
 * Ask for the familiar name, and do not mind being refused.
 *
 * Renaming is where the check we just avoided would run again, so this is
 * attempted only after close() has succeeded and the object is safe on disk. A
 * refusal leaves a file that works and reads oddly; treating it as an error
 * would throw away a finished conversion over a file extension.
 */
async function promote(handle: SaveHandle): Promise<string> {
  const want = handle.name.replace(new RegExp(`${DATA_EXT}$`), '.zip')
  if (want === handle.name || typeof handle.move !== 'function') return handle.name
  try {
    await handle.move(want)
    return want
  } catch {
    return handle.name
  }
}

class SaveFailed extends Error {
  readonly written: number
  readonly total: number
  readonly swap: string
  constructor(e: unknown, written: number, total: number, swap: string) {
    super(e instanceof Error ? e.message : String(e))
    this.name = e instanceof Error ? e.name : 'Error'
    this.written = written
    this.total = total
    this.swap = swap
  }
}

function whyNotWritten(e: unknown, bytes: number): string {
  const raw = e instanceof Error ? e.message : String(e)
  const gb = (bytes / 1e9).toFixed(1)
  const name = e instanceof Error ? e.name : ''
  const locked = /NoModificationAllowed|state of the underlying filesystem|in use|locked/i.test(raw)
  const space = /space|quota|disk full|ENOSPC/i.test(raw) || name === 'QuotaExceededError'

  // Every byte is already on disk and only the last step failed. Chrome runs a
  // Safe Browsing check over an archive before putting it in place, and on a
  // file this size that check does not come back — measured on a 5.80 GB
  // conversion whose swap file was complete to the byte. Converting again would
  // take another quarter of an hour and produce the same bytes that are already
  // sitting there, so the thing to say is where they are.
  if (e instanceof SaveFailed && e.written >= e.total) {
    const scan = /safe browsing/i.test(raw)
    return `Every byte was written — all ${gb} GB of it — and only the last step failed`
      + (scan ? ', the check Chrome runs over an archive before putting it in place.' : `: ${raw}.`)
      + ` Your file is finished and sitting beside where you chose to save it, named`
      + ` ${e.swap}. Rename it to ${e.swap.replace(/\.crswap$/, '')} and it is ready for the`
      + ` studio — there is no need to convert again. Pressing the button also works if you`
      + ` would rather have the browser do it.`
  }

  const why = space
    ? `There is not enough room where you chose to save it. This file is ${gb} GB, and while it is being written the browser needs about that much again beside it for a temporary copy.`
    : locked
      ? `That location would not accept the file. The usual reasons are a folder that syncs or is not a plain local disk — OneDrive, a network share, a phone or a memory stick — or a file of the same name already open in another tab or program. The browser also needs about ${gb} GB free there for a temporary copy while it writes.`
      : raw
  // The browser's own words, kept even when they have been translated into the
  // sentence above. Those two branches each cover several causes, so without
  // this a second occurrence is indistinguishable from the first and there is
  // nothing to tell them apart by — which is exactly the position this left us
  // in when it happened twice.
  const said = [name, raw].filter(Boolean).join(': ')
  const detail = said && said !== why ? ` (the browser said: ${said})` : ''
  return `The file was not written. ${why} Nothing was lost — the converted data is still here, so choose a plain folder on your own disk and press the button again.${detail}`
}

/** The name the file ended up with, or null if the picker was dismissed. */
async function saveCollection(
  result: BuiltCollection, onProgress: (done: number) => void,
): Promise<string | null> {
  const pick = (window as unknown as SavePicker).showSaveFilePicker
  if (pick) {
    let handle: SaveHandle
    try {
      handle = await pick({
        suggestedName: result.name.replace(/\.zip$/, DATA_EXT),
        types: [{ description: 'scRNA-seq Studio dataset', accept: { 'application/octet-stream': [DATA_EXT] } }],
      })
    } catch {
      return null   // the user closed the picker; that is not an error
    }
    const w = await handle.createWritable()
    // A part can be hundreds of megabytes; hand it over in 64 MB writes so the
    // stream keeps flowing and there is a number to put on the button.
    const STEP = 64 << 20
    const swap = `${handle.name}.crswap`
    let done = 0
    try {
      for (const piece of collectionPieces(result.meta, result.parts)) {
        for (let at = 0; at < piece.length; at += STEP) {
          await w.write(piece.subarray(at, Math.min(piece.length, at + STEP)))
          done += Math.min(piece.length, at + STEP) - at
          onProgress(done)
        }
      }
    } catch (e) {
      // Failed part way through the writes. The swap file holds a truncated
      // container that is no use to anyone, and an abandoned stream keeps both
      // it and the destination locked — so the next attempt would fail for a
      // different reason than this one. Abort discards it, which is right here
      // and would be very wrong below.
      try { await w.abort() } catch { /* the stream was already gone */ }
      throw new SaveFailed(e, done, result.bytes, swap)
    }

    // Every byte is on disk; all that is left is the rename. Keep the FIRST
    // error if a retry also fails — a second close on an already-errored stream
    // reports the stream's state rather than what went wrong, and what went
    // wrong is the only useful thing here.
    let closed: unknown = null
    try {
      await w.close()
    } catch (first) {
      closed = first
      try {
        await new Promise(r => setTimeout(r, 1500))
        await w.close()
        closed = null
      } catch { /* keep `first` */ }
    }
    // Deliberately NOT aborting: abort() discards the swap file, and the swap
    // file is now the finished object. Leaving it costs the user a rename;
    // discarding it costs them the conversion.
    if (closed) throw new SaveFailed(closed, result.bytes, result.bytes, swap)
    return promote(handle)
  }

  const url = URL.createObjectURL(writeCollection(result.meta, result.parts))
  const a = document.createElement('a')
  a.href = url
  a.download = result.name
  a.click()
  // Long enough for the browser to have started writing; revoking on a short
  // timer is what cancelled this download when it was the primary path.
  setTimeout(() => URL.revokeObjectURL(url), 600_000)
  return result.name
}

export default function App() {
  const [scan, setScan] = useState<ScanInfo | null>(null)
  const [choices, setChoices] = useState<Choices | null>(null)
  // How the lab decided to store the object. Never a question, never a control:
  // the studio opens the result as one dataset with all its cells either way.
  const [split, setSplit] = useState<SplitDecision | null>(null)
  const [result, setResult] = useState<BuiltCollection | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [stage, setStage] = useState('')
  const [frac, setFrac] = useState(0)
  const [busy, setBusy] = useState(false)
  // The name the file actually ended up with, not merely that it was saved —
  // the two differ whenever Chromium refuses to rename a large object into an
  // archive extension, and the user has to be told which one to look for.
  const [saved, setSaved] = useState<string | null>(null)
  // Writing 5.8 GB to disk takes a while and the button must not be pressable
  // twice while it happens.
  const [saving, setSaving] = useState(false)
  const [wrote, setWrote] = useState(0)

  async function open(file: File) {
    setError(null)
    setResult(null)
    setBusy(true)
    try {
      setStage(`Reading ${file.name}`)
      const s = await openFile(file, { onStage: setStage })
      setScan(s)
      setChoices({
        cluster: guessCluster(s) ?? '',
        sample: guessRole(s, 'sample'),
        condition: guessRole(s, 'condition'),
        embedding: guessEmbedding(s) ?? '',
        label: file.name.replace(/\.(h5ad|rds|h5)$/i, ''),
        extras: [],
      })
      setSplit(chooseSplit(s.columns, s.nnzPerCell ?? null, s.nCells))
    } catch (e) {
      setError(explain(e))
    } finally {
      setBusy(false)
      setStage('')
    }
  }

  async function convert() {
    if (!scan || !choices) return
    setBusy(true)
    setError(null)
    setResult(null)
    setSaved(null)
    try {
      const order = split
        ? {
            column: split.column.name,
            reason: split.reason,
            shards: split.plan.shards,
            passes: split.plan.passes.map(p => p.map(s => split.plan.shards.indexOf(s))),
          }
        : null
      // The whole object's cluster order, before any part drops the levels it
      // has no cells for. It decides colour in the studio.
      const clusterOrder = scan.columns.find(c => c.name === choices.cluster)?.levels
      const built = await runConvert(choices, order, {
        onStage: setStage,
        onProgress: (_pass, f) => setFrac(f),
      }, clusterOrder)
      setResult(recordOrders(built, scan, choices))
    } catch (e) {
      setError(explain(e))
    } finally {
      setBusy(false)
      setStage('')
      setFrac(0)
    }
  }

  function restart() {
    release()
    setScan(null)
    setChoices(null)
    setSplit(null)
    setResult(null)
    setError(null)
    setSaved(null)
  }

  if (result) {
    const { meta, name, bytes } = result
    return (
      <div className="wrap py-10">
        <Card
          eyebrow="Converted — the lab's job is done"
          title={`${name} is ready for scRNA-seq Studio`}
          sub={`${meta.nCells.toLocaleString()} cells · ${meta.nGenes.toLocaleString()} genes · ${(bytes / 1e6).toFixed(1)} MB`}
        >
          <Flow at="done" />

          <div className="mt-5 grid gap-2.5">
            <Handoff n={1} title="Save the file" done={!!saved}>
              <button
                className={`btn ${saved ? '' : 'btn-primary'}`}
                disabled={saving}
                onClick={async () => {
                  // The last failure's message goes with the attempt it belongs
                  // to. Left standing, it sits under a button that reads
                  // "Writing…" and tells the user their file was not written
                  // while it is being written — and the message itself is what
                  // told them to press again.
                  setError(null)
                  setSaving(true)
                  setWrote(0)
                  try {
                    const wrote = await saveCollection(result, setWrote)
                    if (wrote) setSaved(wrote)
                  } catch (e) {
                    setError(whyNotWritten(e, result.bytes))
                  } finally {
                    setSaving(false)
                  }
                }}
              >{saving
                ? `Writing… ${(wrote / 1e9).toFixed(1)} of ${(bytes / 1e9).toFixed(1)} GB`
                : saved ? 'Saved ✓ — download again' : `Download ${name}`}</button>
              {error && <p className="sub mt-2" style={{ color: 'var(--bad)' }}>{error}</p>}
              {/* Only when the name is not the one that was offered. Chromium
                  will not rename a file this size into an archive extension —
                  see DATA_EXT — so on a large object the file keeps the name it
                  was written under, and the studio opens it either way. Saying
                  nothing would leave the user hunting for a .zip that is not
                  there. */}
              {saved && !saved.endsWith('.zip') && (
                <p className="sub mt-2">Saved as <b>{saved}</b>. It is a zip inside, and the studio
                  opens it as it is — Chrome will not put a file this large under a{' '}
                  <code>.zip</code> name, because it cannot finish the safety check it runs on an
                  archive. Rename it if you want to look inside with an archive tool.</p>
              )}
              {/* Stated, not offered. The user made no decision here and has
                  nothing to do about it; it is on the page because a 2 GB file
                  that took forty minutes deserves an explanation. */}
              {meta.reason && (
                <p className="sub mt-2">{meta.reason}. It still opens as one dataset with all
                  {' '}{meta.nCells.toLocaleString()} cells.</p>
              )}
            </Handoff>

            <Handoff n={2} title="Open it in scRNA-seq Studio">
              <a className={`btn ${saved ? 'btn-primary' : ''}`} href={STUDIO}
                target="_blank" rel="noreferrer">Open scRNA-seq Studio →</a>
              <p className="sub mt-2">
                Drop it onto the panel — it stays on your machine.
              </p>
            </Handoff>
          </div>

          <div className="mt-4">
            <button className="btn btn-ghost" onClick={restart}>Convert another object</button>
          </div>

          <div className="mt-5">
            <div className="eyebrow mb-2">Every decision the converter made</div>
            <ul className="grid gap-1">
              {meta.notes.map((n, i) => (
                <li key={i} className="text-[12.5px]" style={{ color: 'var(--ink-2)' }}>· {n}</li>
              ))}
            </ul>
            <p className="sub mt-2">
              These travel inside the file and appear again on the studio&rsquo;s Overview tab.
            </p>
          </div>
        </Card>
      </div>
    )
  }

  if (scan && choices) {
    return (
      <>
        {error && (
          <div className="wrap pt-6">
            <div className="note"><b>That did not convert.</b> {error}</div>
          </div>
        )}
        {/* Above the questions rather than after them: the Convert button is the
            last thing on the Review card, and an optional step below it is a
            step nobody sees. */}
        <Extras scan={scan} choices={choices} setChoices={setChoices} />
        <Review
          scan={scan} choices={choices} setChoices={setChoices}
          onBuild={convert} busy={busy}
          stage={stage + (frac > 0 ? ` · ${(frac * 100).toFixed(0)}%` : '')}
        />
        <div className="wrap pb-10">
          <button className="btn btn-ghost" onClick={restart}>Open a different file</button>
        </div>
      </>
    )
  }

  return <Landing onFile={open} error={error} busy={busy} stage={stage} />
}
