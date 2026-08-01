import { useRef, useState, type ReactNode } from 'react'

const STUDIO = 'https://jiaenlin.github.io/scrnaseq-studio/'

/** What the object has to contain, and what happens when it does not. */
const REQUIREMENTS: [string, ReactNode, ReactNode][] = [
  ['A cell annotation', 'any categorical column — cell_type, seurat_annotations, louvain, leiden',
    <>Required. You pick which one after the scan.</>],
  ['An embedding', 'UMAP, t-SNE or PCA in obsm / @reductions',
    <>Required by the Cells and Feature plot tabs.</>],
  ['Expression', 'log-normalized values, or raw counts to build them from',
    <>Required. Scaled values are z-scores; the converter refuses them rather than plotting
      violins that look normal and are wrong.</>],
  ['A sample column', 'donor, animal, run — orig.ident and the like',
    <>Optional. Without it composition cannot show between-animal spread.</>],
  ['A condition column', 'the experimental group',
    <>Optional. Without it the object opens single-condition.</>],
]

export default function Landing({ onFile, error, busy, stage }: {
  onFile: (f: File) => void
  error: string | null
  busy: boolean
  stage: string
}) {
  const input = useRef<HTMLInputElement>(null)
  const [over, setOver] = useState(false)

  const take = (list: FileList | null) => {
    const f = list?.[0]
    if (f) onFile(f)
  }

  return (
    <div className="grid min-h-screen place-items-center px-5 py-10">
      <div className="w-full max-w-[700px]">
        <div
          className="rounded-[18px] px-[30px] py-[34px] text-center"
          style={{
            border: `1.5px dashed ${over ? 'var(--accent)' : 'var(--line-2)'}`,
            background: over ? 'var(--accent-soft)' : 'var(--surface)',
            transition: 'border-color 160ms ease, background-color 160ms ease',
          }}
          onDragOver={e => { e.preventDefault(); setOver(true) }}
          onDragLeave={() => setOver(false)}
          onDrop={e => { e.preventDefault(); setOver(false); take(e.dataTransfer.files) }}
        >
          <div
            className="mx-auto mb-3.5 grid h-[30px] w-[30px] place-items-center rounded-[9px] text-xs font-bold text-white"
            style={{ background: 'linear-gradient(150deg, var(--accent), #a855f7)' }}
          >sc</div>
          <h1 className="text-[19px] tracking-[-0.02em]">scRNA-seq Lab</h1>
          <p className="mb-[18px] mt-2 text-[13.5px]" style={{ color: 'var(--ink-2)' }}>
            Turn a Scanpy <code className="mono">.h5ad</code> or a Seurat{' '}
            <code className="mono">.rds</code> into a studio bundle. It reads the object here,
            in this tab — nothing is uploaded, and there is no Python or R to install.
          </p>

          <input
            ref={input} type="file" accept=".h5ad,.rds,.h5,.RDS" className="hidden"
            onChange={e => take(e.target.files)}
          />
          <button className="btn btn-primary" disabled={busy} onClick={() => input.current?.click()}>
            {busy ? (stage || 'Reading…') : 'Choose a .h5ad or .rds'}
          </button>
          <p className="mt-2 text-[12px]" style={{ color: 'var(--ink-3)' }}>
            or drop it anywhere on this panel
          </p>

          {busy && stage && (
            <div className="note note-info mt-4 text-left"><b>{stage}</b></div>
          )}

          {error && (
            <div className="note mt-4 text-left">
              <b>That file did not open.</b> {error}
            </div>
          )}

          <div className="mt-5 rounded-xl px-4 py-3 text-left" style={{ background: 'var(--sunk)' }}>
            <div className="eyebrow">How this works</div>
            <p className="mt-1.5 text-[12.5px]" style={{ color: 'var(--ink-2)' }}>
              The lab scans the object&rsquo;s metadata, shows you every column it found, and asks
              which one is the cell type annotation and which — if any — is the condition to group
              by. Then it writes <code className="mono">bundle.zip</code>, which you open in{' '}
              <a className="underline" href={STUDIO} target="_blank" rel="noreferrer">
                scRNA-seq Studio</a>.
            </p>
            <p className="mt-2 text-[12.5px]" style={{ color: 'var(--ink-2)' }}>
              Large objects are fine: a 288&nbsp;MB Seurat object scans in about a second, because
              the reader walks the file without expanding the matrices it is not going to use.
            </p>
          </div>

          <div className="eyebrow mt-6 text-left">What the object needs</div>
          <div className="scrollx mt-2">
            <table className="t">
              <thead><tr><th>Needs</th><th>Which is</th><th>If missing</th></tr></thead>
              <tbody>
                {REQUIREMENTS.map(([name, what, missing]) => (
                  <tr key={name}>
                    <td className="whitespace-nowrap font-semibold">{name}</td>
                    <td style={{ color: 'var(--ink-2)' }}>{what}</td>
                    <td style={{ color: 'var(--ink-2)' }}>{missing}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <p className="mt-4 text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
            Raw FASTQ instead? Start at <b>rnaseq-service</b>. A bulk count matrix goes to{' '}
            <b>rnaseq-lab</b>, then <b>rnaseq-studio</b>.
          </p>
        </div>
      </div>
    </div>
  )
}
