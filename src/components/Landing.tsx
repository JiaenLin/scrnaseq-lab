import { useRef, useState } from 'react'
import { Flow } from './Ui.tsx'

const STUDIO = 'https://jiaenlin.github.io/scrnaseq-studio/'

/**
 * What the object has to contain.
 *
 * Folded away by default. It is reference material — the app enforces every row
 * at conversion time and says which one failed — so putting it between someone
 * and the file picker made them read a five-row table to do nothing at all.
 */
const NEEDS: [string, boolean, string][] = [
  ['Cell annotation', true, 'any categorical column — you pick which one after the scan'],
  ['Embedding', true, 'UMAP, t-SNE or PCA'],
  ['Expression', true, 'log-normalized, or raw counts to build it from — scaled z-scores are refused'],
  ['Sample', false, 'donor, animal or run — without it, composition cannot show spread'],
  ['Condition', false, 'the experimental group — without it, the comparison tabs stay empty'],
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
      <div className="w-full max-w-[540px]">
        <header className="text-center">
          <div
            className="mx-auto mb-3 grid h-[30px] w-[30px] place-items-center rounded-[9px] text-xs font-bold text-white"
            style={{ background: 'linear-gradient(150deg, var(--accent), #a855f7)' }}
          >sc</div>
          <h1 className="text-[20px] tracking-[-0.02em]">scRNA-seq Lab</h1>
          <p className="mx-auto mt-1.5 max-w-[430px] text-[13px]" style={{ color: 'var(--ink-2)' }}>
            Turns an annotated single-cell object into the file{' '}
            {/* nowrap: the hyphen in "scRNA-seq" is a break opportunity, so the
                product name was splitting across two lines mid-word. */}
            <a className="underline" href={STUDIO} target="_blank" rel="noreferrer"
              style={{ whiteSpace: 'nowrap' }}>scRNA-seq Studio</a>{' '}opens.
          </p>
          <div className="mt-3.5"><Flow at="convert" /></div>
        </header>

        <div
          className="mt-6 rounded-2xl px-6 py-9 text-center"
          style={{
            border: `1.5px dashed ${over ? 'var(--accent)' : 'var(--line-2)'}`,
            background: over ? 'var(--accent-soft)' : 'var(--surface)',
            transition: 'border-color 160ms ease, background-color 160ms ease',
          }}
          onDragOver={e => { e.preventDefault(); setOver(true) }}
          onDragLeave={() => setOver(false)}
          onDrop={e => { e.preventDefault(); setOver(false); take(e.dataTransfer.files) }}
        >
          <input
            ref={input} type="file" accept=".h5ad,.rds,.h5,.RDS" className="hidden"
            onChange={e => take(e.target.files)}
          />
          <button className="btn btn-primary" disabled={busy} onClick={() => input.current?.click()}>
            {busy ? (stage || 'Reading…') : 'Choose a file'}
          </button>
          <p className="mt-2.5 text-[12px]" style={{ color: 'var(--ink-3)' }}>
            or drop it here — <code className="mono">.h5ad</code> or{' '}
            <code className="mono">.rds</code>
          </p>
        </div>

        {busy && stage && <div className="note note-info mt-3"><b>{stage}</b></div>}
        {error && (
          <div className="note mt-3"><b>That file did not open.</b> {error}</div>
        )}

        <p className="mt-3 text-center text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
          Read in this tab — nothing is uploaded, and there is no Python or R to install.
        </p>

        <details className="mt-5 rounded-xl px-3.5 py-2.5" style={{ background: 'var(--sunk)' }}>
          <summary className="cursor-pointer text-[12.5px] font-semibold">
            What the object needs
          </summary>
          <div className="mt-2.5 grid gap-1.5">
            {NEEDS.map(([name, required, what]) => (
              <div key={name} className="flex items-baseline gap-2">
                <span className="w-[104px] flex-none text-[12px] font-semibold">{name}</span>
                <span className={`badge flex-none ${required ? 'badge-here' : 'badge-none'}`}>
                  {required ? 'required' : 'optional'}
                </span>
                <span className="text-[11.5px]" style={{ color: 'var(--ink-2)' }}>{what}</span>
              </div>
            ))}
          </div>
          <p className="mt-2.5 text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
            The lab does not cluster or annotate — bring an object you have already processed.
            After the scan it asks which column is the cell annotation and which is the condition.
          </p>
        </details>

        <p className="mt-4 text-center text-[11px]" style={{ color: 'var(--ink-3)' }}>
          Bulk RNA-seq? Counts go to <b>rnaseq-lab</b>; raw FASTQ starts at <b>rnaseq-service</b>.
        </p>
      </div>
    </div>
  )
}
