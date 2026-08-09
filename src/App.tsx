import { useState, type ReactNode } from 'react'
import Landing from './components/Landing.tsx'
import Review from './components/Review.tsx'
import { Card, Flow } from './components/Ui.tsx'
import {
  convert as runConvert, openFile, release, type BuiltCollection,
} from './lib/client.ts'
import { collectionPieces, writeCollection } from './lib/collection.ts'
import { chooseSplit, type SplitDecision } from './lib/shard.ts'
import { guessCluster, guessEmbedding, guessRole, type Choices, type ScanInfo } from './lib/scan.ts'

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
 * Minimal shape of the File System Access API, which TypeScript's DOM library
 * does not declare. Only the three calls used below.
 */
interface SavePicker {
  showSaveFilePicker?: (o: {
    suggestedName?: string
    types?: { description: string; accept: Record<string, string[]> }[]
  }) => Promise<{ createWritable: () => Promise<{
    write: (d: BufferSource | Blob) => Promise<void>; close: () => Promise<void>
  }> }>
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
async function saveCollection(
  result: BuiltCollection, onProgress: (done: number) => void,
): Promise<boolean> {
  const pick = (window as unknown as SavePicker).showSaveFilePicker
  if (pick) {
    let handle
    try {
      handle = await pick({
        suggestedName: result.name,
        types: [{ description: 'scRNA-seq Studio dataset', accept: { 'application/zip': ['.zip'] } }],
      })
    } catch {
      return false   // the user closed the picker; that is not an error
    }
    const w = await handle.createWritable()
    // A part can be hundreds of megabytes; hand it over in 64 MB writes so the
    // stream keeps flowing and there is a number to put on the button.
    const STEP = 64 << 20
    let done = 0
    for (const piece of collectionPieces(result.meta, result.parts)) {
      for (let at = 0; at < piece.length; at += STEP) {
        await w.write(piece.subarray(at, Math.min(piece.length, at + STEP)))
        done += Math.min(piece.length, at + STEP) - at
        onProgress(done)
      }
    }
    await w.close()
    return true
  }

  const url = URL.createObjectURL(writeCollection(result.meta, result.parts))
  const a = document.createElement('a')
  a.href = url
  a.download = result.name
  a.click()
  // Long enough for the browser to have started writing; revoking on a short
  // timer is what cancelled this download when it was the primary path.
  setTimeout(() => URL.revokeObjectURL(url), 600_000)
  return true
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
  const [saved, setSaved] = useState(false)
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
    setSaved(false)
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
      setResult(await runConvert(choices, order, {
        onStage: setStage,
        onProgress: (_pass, f) => setFrac(f),
      }, clusterOrder))
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
    setSaved(false)
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
            <Handoff n={1} title="Save the file" done={saved}>
              <button
                className={`btn ${saved ? '' : 'btn-primary'}`}
                disabled={saving}
                onClick={async () => {
                  setSaving(true)
                  setWrote(0)
                  try {
                    if (await saveCollection(result, setWrote)) setSaved(true)
                  } catch (e) {
                    setError(`the file could not be written. ${explain(e)}`)
                  } finally {
                    setSaving(false)
                  }
                }}
              >{saving
                ? `Writing… ${(wrote / 1e9).toFixed(1)} of ${(bytes / 1e9).toFixed(1)} GB`
                : saved ? 'Saved ✓ — download again' : `Download ${name}`}</button>
              {error && <p className="sub mt-2" style={{ color: 'var(--bad)' }}>{error}</p>}
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
