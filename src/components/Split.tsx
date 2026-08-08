import { useMemo } from 'react'
import type { ScanInfo } from '../lib/scan.ts'
import { COMFORTABLE_NNZ, planSplit, splitCandidates, type SplitPlan } from '../lib/shard.ts'
import { Card } from './Ui.tsx'

const M = (n: number) => `${(n / 1e6).toFixed(n < 1e7 ? 1 : 0)} M`

/** Roughly what a bundle of this many nonzeros weighs on disk. */
const sizeOf = (nnz: number) => `${Math.max(1, Math.round(nnz * 8 / 1.9e6))} MB`

export default function Split({ scan, splitBy, onChange }: {
  scan: ScanInfo
  splitBy: string | null
  onChange: (col: string | null) => void
}) {
  const nnzPerCell = scan.nnzPerCell ?? null
  const totalNnz = useMemo(() => {
    if (!nnzPerCell) return 0
    let t = 0
    for (const v of nnzPerCell) t += v
    return t
  }, [nnzPerCell])

  const options = useMemo(
    () => splitCandidates(scan.columns, nnzPerCell, scan.nCells).slice(0, 6),
    [scan, nnzPerCell])

  // Nothing to decide: the object already fits, and splitting a small object
  // into pieces only makes it harder to look at.
  if (!totalNnz || totalNnz <= COMFORTABLE_NNZ) return null

  const chosen: SplitPlan | null = splitBy
    ? planSplit(scan.columns.find(c => c.name === splitBy)!, nnzPerCell, scan.nCells)
    : null

  return (
    <Card
      eyebrow="This object is too large for one bundle"
      title="Split it into several"
      sub={`${M(totalNnz)} stored values — the studio holds the matrix in memory, so it opens about ${M(COMFORTABLE_NNZ)} comfortably. Split along a column and each piece opens on its own, with every tab working.`}
    >
      <div className="mt-3.5 grid gap-2">
        {options.map(({ column: c, plan }) => {
          const on = splitBy === c.name
          const ok = plan.oversized.length === 0
          return (
            <button
              key={c.name} onClick={() => onChange(c.name)}
              className="w-full rounded-xl px-3 py-2.5 text-left"
              style={{
                border: `1.5px solid ${on ? 'var(--accent)' : 'var(--line)'}`,
                background: on ? 'var(--accent-soft)' : 'var(--surface)',
              }}
            >
              <div className="flex flex-wrap items-baseline gap-2">
                <span className="mono text-[12.5px] font-semibold">{c.name}</span>
                <span className="badge badge-file">{plan.shards.length} bundles</span>
                <span className={`badge ${ok ? 'badge-here' : 'badge-none'}`}>
                  {ok ? 'all open' : `${plan.oversized.length} still too big`}
                </span>
                <span className="text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
                  largest {M(plan.shards[0].nnz)} · {plan.passes.length} pass
                  {plan.passes.length > 1 ? 'es' : ''} over the file
                </span>
              </div>
            </button>
          )
        })}
        <button
          onClick={() => onChange(null)}
          className="w-full rounded-xl px-3 py-2.5 text-left"
          style={{
            border: `1.5px solid ${splitBy === null ? 'var(--accent)' : 'var(--line)'}`,
            background: splitBy === null ? 'var(--accent-soft)' : 'var(--surface)',
          }}
        >
          <div className="text-[12.5px] font-semibold">Do not split</div>
          <div className="mt-1 text-[11.5px]" style={{ color: 'var(--ink-2)' }}>
            One bundle of {M(totalNnz)} values, about {sizeOf(totalNnz)}. It will convert, and the
            studio will most likely fail to open it.
          </div>
        </button>
      </div>

      {chosen && (
        <div className="mt-4">
          <div className="eyebrow mb-2">
            {chosen.shards.length} bundles from {chosen.column} · {chosen.passes.length} pass
            {chosen.passes.length > 1 ? 'es' : ''}, about a minute each
          </div>
          <div className="scrollx" style={{ maxHeight: 240 }}>
            <table className="t">
              <thead><tr><th>Bundle</th><th className="num">Cells</th><th className="num">Values</th><th className="num">≈ size</th></tr></thead>
              <tbody>
                {chosen.shards.map(s => (
                  <tr key={s.key}>
                    <td className="font-semibold">{s.key}</td>
                    <td className="num">{s.cells.length.toLocaleString()}</td>
                    <td className="num" style={{ color: s.nnz > COMFORTABLE_NNZ ? 'var(--bad)' : 'var(--ink-2)' }}>
                      {M(s.nnz)}
                    </td>
                    <td className="num" style={{ color: 'var(--ink-2)' }}>{sizeOf(s.nnz)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <p className="sub mt-2">
            Sizes are exact, counted from the row offsets before any values were read. Every
            bundle keeps the full gene list and all {scan.columns.length} metadata columns — only
            the cells differ.
          </p>
        </div>
      )}
    </Card>
  )
}
