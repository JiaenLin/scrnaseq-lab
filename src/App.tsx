import { useState, type ReactNode } from 'react'
import { gunzipSync } from 'fflate'
import Landing from './components/Landing.tsx'
import Review from './components/Review.tsx'
import { Card, Flow } from './components/Ui.tsx'
import { decompressRds, RdsReader } from './lib/rds.ts'
import { scanSeurat } from './lib/seurat.ts'
import { scanH5ad } from './lib/h5ad.ts'
import { buildBundle, type BundleMeta } from './lib/build.ts'
import { guessCluster, guessEmbedding, guessRole, type Choices, type Scan } from './lib/scan.ts'

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

/** Let the browser paint the stage label before the next blocking step. */
const yieldToPaint = () =>
  new Promise<void>(r => requestAnimationFrame(() => setTimeout(r, 0)))

interface Done {
  blob: Blob
  name: string
  meta: BundleMeta
  pseudobulkColumns: number
}

export default function App() {
  const [scan, setScan] = useState<Scan | null>(null)
  const [choices, setChoices] = useState<Choices | null>(null)
  const [done, setDone] = useState<Done | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [stage, setStage] = useState('')
  const [busy, setBusy] = useState(false)
  const [saved, setSaved] = useState(false)

  async function open(file: File) {
    setError(null)
    setDone(null)
    setBusy(true)
    try {
      const name = file.name.toLowerCase()
      setStage(`Reading ${file.name}`)
      await yieldToPaint()
      const bytes = new Uint8Array(await file.arrayBuffer())

      let s: Scan
      if (name.endsWith('.rds')) {
        setStage('Decompressing the .rds')
        await yieldToPaint()
        const raw = decompressRds(bytes, b => gunzipSync(b))
        setStage('Walking the R object')
        await yieldToPaint()
        const r = new RdsReader(raw)
        s = scanSeurat(r, r.read(), file.name)
      } else if (name.endsWith('.h5ad') || name.endsWith('.h5')) {
        setStage('Opening the HDF5 file')
        await yieldToPaint()
        // Loaded on demand: h5wasm is a WASM build of HDF5 and has no business
        // being in the bundle of someone who only ever opens .rds files.
        const h5wasm = (await import('h5wasm')).default
        setStage('Scanning obs, var and obsm')
        await yieldToPaint()
        s = await scanH5ad(h5wasm, bytes, file.name)
      } else {
        throw new Error(`this is a ${name.split('.').pop() ?? 'file'} — the lab reads .h5ad (Scanpy/AnnData) and .rds (Seurat).`)
      }

      setScan(s)
      setChoices({
        cluster: guessCluster(s) ?? '',
        sample: guessRole(s, 'sample'),
        condition: guessRole(s, 'condition'),
        embedding: guessEmbedding(s) ?? '',
        label: file.name.replace(/\.(h5ad|rds|h5)$/i, ''),
      })
    } catch (e) {
      setError((e as Error).message)
    } finally {
      setBusy(false)
      setStage('')
    }
  }

  async function convert() {
    if (!scan || !choices) return
    setBusy(true)
    setError(null)
    try {
      let out: ReturnType<typeof buildBundle> | null = null
      // buildBundle reports each stage as it starts; yielding between them is
      // what makes those labels visible rather than one long freeze.
      const stages: string[] = []
      await new Promise<void>((resolve, reject) => {
        const run = async () => {
          try {
            out = buildBundle(scan, choices, st => {
              stages.push(st)
              setStage(st)
            })
            resolve()
          } catch (e) { reject(e) }
        }
        setStage('Starting')
        requestAnimationFrame(() => void run())
      })
      const res = out! as ReturnType<typeof buildBundle>
      const blobPart = new Uint8Array(res.zip)
      setDone({
        blob: new Blob([blobPart as unknown as BlobPart], { type: 'application/zip' }),
        name: `${(choices.label || 'bundle').replace(/[^\w.-]+/g, '_')}.zip`,
        meta: res.meta,
        pseudobulkColumns: res.pseudobulkColumns,
      })
    } catch (e) {
      setError((e as Error).message)
    } finally {
      setBusy(false)
      setStage('')
    }
  }

  const download = () => {
    if (!done) return
    const url = URL.createObjectURL(done.blob)
    const a = document.createElement('a')
    a.href = url
    a.download = done.name
    a.click()
    setSaved(true)
    setTimeout(() => URL.revokeObjectURL(url), 5000)
  }

  const restart = () => {
    setScan(null)
    setChoices(null)
    setDone(null)
    setError(null)
    setSaved(false)
  }

  if (done) {
    const mb = done.blob.size / 1e6
    return (
      <div className="wrap py-10">
        <Card
          eyebrow="Converted — the lab's job is done"
          title={`${done.name} is ready for scRNA-seq Studio`}
          sub={`${done.meta.nCells.toLocaleString()} cells · ${done.meta.nGenes.toLocaleString()} genes · ${done.meta.nnz.toLocaleString()} stored values · ${mb.toFixed(1)} MB`}
        >
          <Flow at="done" />

          {/* Two steps, and only one of them is live at a time: the studio cannot
              open a file that has not been saved yet, so it stays secondary
              until the download has actually happened. */}
          <div className="mt-5 grid gap-2.5">
            <Handoff n={1} title="Save the bundle" done={saved}>
              <button className={`btn ${saved ? '' : 'btn-primary'}`} onClick={download}>
                {saved ? `Download ${done.name} again` : `Download ${done.name}`}
              </button>
            </Handoff>
            <Handoff n={2} title="Open it in scRNA-seq Studio">
              <a className={`btn ${saved ? 'btn-primary' : ''}`} href={STUDIO}
                target="_blank" rel="noreferrer">
                Open scRNA-seq Studio →
              </a>
              <p className="sub mt-2">
                Drop <span className="mono">{done.name}</span> onto its panel — it stays on your
                machine.
              </p>
            </Handoff>
          </div>

          <div className="mt-4">
            <button className="btn btn-ghost" onClick={restart}>Convert another object</button>
          </div>

          <div className="mt-5">
            <div className="eyebrow mb-2">What went in</div>
            <div className="scrollx">
              <table className="t">
                <tbody>
                  <tr><td className="font-semibold">Cell types</td>
                    <td>{done.meta.clusters.length} — {done.meta.clusters.join(', ')}</td></tr>
                  <tr><td className="font-semibold">Conditions</td>
                    <td>{done.meta.conditions.join(', ')}</td></tr>
                  <tr><td className="font-semibold">Samples</td>
                    <td>{done.meta.samples.length} — {done.meta.samples.slice(0, 8).map(s => s.id).join(', ')}
                      {done.meta.samples.length > 8 ? ' …' : ''}</td></tr>
                  <tr><td className="font-semibold">Embedding</td><td className="mono">{done.meta.embedding}</td></tr>
                  <tr><td className="font-semibold">Expression</td><td className="mono">{done.meta.expression}</td></tr>
                  <tr><td className="font-semibold">Pseudobulk</td>
                    <td>{done.pseudobulkColumns
                      ? `${done.pseudobulkColumns} sample × cluster columns of summed raw counts`
                      : 'not written — no raw counts in the object'}</td></tr>
                </tbody>
              </table>
            </div>
          </div>

          <div className="mt-5">
            <div className="eyebrow mb-2">Every decision the converter made</div>
            <ul className="grid gap-1">
              {done.meta.notes.map((n, i) => (
                <li key={i} className="text-[12.5px]" style={{ color: 'var(--ink-2)' }}>· {n}</li>
              ))}
            </ul>
            <p className="sub mt-2">
              These travel inside the bundle and show again on the studio&rsquo;s Overview tab.
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
          onBuild={convert} busy={busy} stage={stage}
        />
        <div className="wrap pb-10">
          <button className="btn btn-ghost" onClick={restart}>Open a different file</button>
        </div>
      </>
    )
  }

  return <Landing onFile={open} error={error} busy={busy} stage={stage} />
}
