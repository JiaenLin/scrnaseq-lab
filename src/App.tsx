import { useState, type ReactNode } from 'react'
import Landing from './components/Landing.tsx'
import Review from './components/Review.tsx'
import Split from './components/Split.tsx'
import { Card, Flow } from './components/Ui.tsx'
import { convert as runConvert, openFile, release, type BuiltBundle } from './lib/client.ts'
import { COMFORTABLE_NNZ, planSplit, splitCandidates } from './lib/shard.ts'
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

function saveBlob(b: Blob, name: string) {
  const url = URL.createObjectURL(b)
  const a = document.createElement('a')
  a.href = url
  a.download = name
  a.click()
  setTimeout(() => URL.revokeObjectURL(url), 10_000)
}

export default function App() {
  const [scan, setScan] = useState<ScanInfo | null>(null)
  const [choices, setChoices] = useState<Choices | null>(null)
  const [splitBy, setSplitBy] = useState<string | null>(null)
  const [bundles, setBundles] = useState<BuiltBundle[]>([])
  const [finished, setFinished] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [stage, setStage] = useState('')
  const [frac, setFrac] = useState(0)
  const [busy, setBusy] = useState(false)
  const [saved, setSaved] = useState<Set<string>>(new Set())

  async function open(file: File) {
    setError(null)
    setBundles([])
    setFinished(false)
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
      // Pre-select the best split when the object is too large to be one
      // bundle. Split renders nothing at all when it already fits, so a normal
      // object never sees any of this.
      let total = 0
      for (const v of s.nnzPerCell ?? []) total += v
      setSplitBy(total > COMFORTABLE_NNZ
        ? splitCandidates(s.columns, s.nnzPerCell ?? null, s.nCells)[0]?.column.name ?? null
        : null)
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
    setBundles([])
    setFinished(false)
    try {
      const col = splitBy ? scan.columns.find(c => c.name === splitBy) : null
      const plan = col ? planSplit(col, scan.nnzPerCell ?? null, scan.nCells) : null
      const split = plan
        ? {
            shards: plan.shards,
            passes: plan.passes.map(p => p.map(s => plan.shards.indexOf(s))),
          }
        : null
      await runConvert(choices, split, {
        onStage: setStage,
        onProgress: (_pass, f) => setFrac(f),
        onBundle: b => setBundles(prev => [...prev, b]),
      })
      setFinished(true)
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
    setBundles([])
    setFinished(false)
    setError(null)
    setSaved(new Set())
  }

  function save(b: BuiltBundle) {
    saveBlob(b.blob, b.name)
    setSaved(prev => new Set(prev).add(b.name))
  }

  if (finished && bundles.length > 0) {
    const totalMb = bundles.reduce((n, b) => n + b.blob.size, 0) / 1e6
    const cells = bundles.reduce((n, b) => n + b.meta.nCells, 0)
    const allSaved = bundles.every(b => saved.has(b.name))
    const many = bundles.length > 1
    return (
      <div className="wrap py-10">
        <Card
          eyebrow="Converted — the lab's job is done"
          title={many
            ? `${bundles.length} bundles are ready for scRNA-seq Studio`
            : `${bundles[0].name} is ready for scRNA-seq Studio`}
          sub={`${cells.toLocaleString()} cells · ${bundles[0].meta.nGenes.toLocaleString()} genes · ${totalMb.toFixed(1)} MB in total`}
        >
          <Flow at="done" />

          <div className="mt-5 grid gap-2.5">
            <Handoff n={1} title={many ? 'Save the bundles' : 'Save the bundle'} done={allSaved}>
              {many && (
                <button
                  className={`btn mb-2.5 ${allSaved ? '' : 'btn-primary'}`}
                  onClick={() => bundles.forEach((b, i) => setTimeout(() => save(b), i * 350))}
                >
                  Download all {bundles.length}
                </button>
              )}
              <div className="scrollx" style={{ maxHeight: 320 }}>
                <table className="t">
                  <thead>
                    <tr>
                      <th>Bundle</th><th className="num">Cells</th>
                      <th className="num">Values</th><th className="num">Size</th><th />
                    </tr>
                  </thead>
                  <tbody>
                    {bundles.map(b => (
                      <tr key={b.name}>
                        <td className="font-semibold">{b.key ?? b.meta.label}</td>
                        <td className="num">{b.meta.nCells.toLocaleString()}</td>
                        <td className="num" style={{ color: 'var(--ink-2)' }}>
                          {(b.meta.nnz / 1e6).toFixed(1)} M
                        </td>
                        <td className="num" style={{ color: 'var(--ink-2)' }}>
                          {(b.blob.size / 1e6).toFixed(1)} MB
                        </td>
                        <td className="num">
                          <button className="chip" onClick={() => save(b)}>
                            {saved.has(b.name) ? 'saved ✓' : 'download'}
                          </button>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </Handoff>

            <Handoff n={2} title={many ? 'Open them in scRNA-seq Studio' : 'Open it in scRNA-seq Studio'}>
              <a className={`btn ${allSaved ? 'btn-primary' : ''}`} href={STUDIO}
                target="_blank" rel="noreferrer">Open scRNA-seq Studio →</a>
              <p className="sub mt-2">
                {many
                  ? 'Drop one bundle at a time onto its panel. Each opens with every tab working on its own cells, and the gene list and metadata columns are identical across all of them.'
                  : 'Drop it onto the panel — it stays on your machine.'}
              </p>
            </Handoff>
          </div>

          <div className="mt-4">
            <button className="btn btn-ghost" onClick={restart}>Convert another object</button>
          </div>

          <div className="mt-5">
            <div className="eyebrow mb-2">Every decision the converter made</div>
            <ul className="grid gap-1">
              {bundles[0].meta.notes.map((n, i) => (
                <li key={i} className="text-[12.5px]" style={{ color: 'var(--ink-2)' }}>· {n}</li>
              ))}
            </ul>
            <p className="sub mt-2">
              These travel inside every bundle and appear again on the studio&rsquo;s Overview tab.
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
        <div className="wrap pt-7">
          <Split scan={scan} splitBy={splitBy} onChange={setSplitBy} />
        </div>
        <Review
          scan={scan} choices={choices} setChoices={setChoices}
          onBuild={convert} busy={busy}
          stage={stage + (frac > 0 ? ` · ${(frac * 100).toFixed(0)}%` : '')}
        />
        {busy && bundles.length > 0 && (
          <div className="wrap pb-4">
            <div className="note note-info">
              <b>{bundles.length} written so far</b> — {bundles.map(b => b.key).join(', ')}
            </div>
          </div>
        )}
        <div className="wrap pb-10">
          <button className="btn btn-ghost" onClick={restart}>Open a different file</button>
        </div>
      </>
    )
  }

  return <Landing onFile={open} error={error} busy={busy} stage={stage} />
}
