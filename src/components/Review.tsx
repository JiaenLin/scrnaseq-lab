import { useMemo } from 'react'
import {
  candidates, column, crossCheck, duplicateOf, replicatesPerCondition,
  samePartition, unusable,
  type Choices, type Column, type Role, type Scan,
} from '../lib/scan.ts'
import { Card, Field, LevelBar, Mono } from './Ui.tsx'

const MIN_REPS_PB = 4

const counts = (c: Column): number[] => {
  const out = new Array(c.levels.length).fill(0)
  for (const v of c.codes) if (v >= 0) out[v]++
  return out
}

/** One selectable column, with the shape of its grouping shown inline. */
function Option({ col, chosen, sameAs, onPick }: {
  col: Column; chosen: boolean; sameAs?: string; onPick: () => void
}) {
  const n = useMemo(() => counts(col), [col])
  const total = col.codes.length
  const preview = col.levels
    .map((l, i) => ({ l, n: n[i] }))
    .sort((a, b) => b.n - a.n)
    .slice(0, 4)

  return (
    <button
      onClick={onPick}
      className="w-full rounded-xl px-3 py-2.5 text-left"
      style={{
        border: `1.5px solid ${chosen ? 'var(--accent)' : 'var(--line)'}`,
        background: chosen ? 'var(--accent-soft)' : 'var(--surface)',
      }}
    >
      <div className="flex items-baseline gap-2">
        <span className="mono text-[12.5px] font-semibold">{col.name}</span>
        <span className="badge badge-file">{col.levels.length} levels</span>
        {col.missing > 0 && (
          <span className="badge badge-none">{col.missing} NA</span>
        )}
      </div>
      {sameAs && (
        <div className="mt-1 text-[11px]" style={{ color: 'var(--ink-3)' }}>
          same grouping as <span className="mono">{sameAs}</span>, relabelled
        </div>
      )}
      <LevelBar counts={n} total={total} />
      <div className="mt-1.5 text-[11.5px]" style={{ color: 'var(--ink-2)' }}>
        {preview.map(p => `${p.l} (${p.n})`).join(' · ')}
        {col.levels.length > 4 && ` · +${col.levels.length - 4} more`}
      </div>
    </button>
  )
}

/** A column already spoken for by another question, and what it was picked as. */
interface Taken { name: string | null; as: string }

function Picker({ scan, role, value, onChange, allowNone, noneLabel, emptyNote, exclude = [] }: {
  scan: Scan; role: Role; value: string | null
  onChange: (v: string | null) => void
  allowNone?: boolean; noneLabel?: string; emptyNote?: string; exclude?: Taken[]
}) {
  const usable = candidates(scan, role, exclude.map(t => t.name))
  const dup = useMemo(() => duplicateOf(scan.columns), [scan])

  /** Why this column is not on offer here — the role clash, or its own shape. */
  const why = (c: Column): string | null => {
    for (const t of exclude) {
      const other = column(scan, t.name)
      if (!other) continue
      if (other.name === c.name) return `already chosen as the ${t.as}`
      if (samePartition(other, c)) {
        return `the same grouping as ${other.name}, already chosen as the ${t.as}`
      }
    }
    return unusable(c, role, scan.nCells)
  }
  const rejected = scan.columns
    .map(c => ({ c, why: why(c) }))
    .filter(x => x.why)

  const rejectedList = rejected.length > 0 && (
    <details className="mt-2">
      <summary className="cursor-pointer text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
        {rejected.length} column{rejected.length > 1 ? 's' : ''} cannot be used here
      </summary>
      <div className="mt-1.5 grid gap-1">
        {rejected.map(({ c, why }) => (
          <div key={c.name} className="text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
            <span className="mono">{c.name}</span> — {why}
          </div>
        ))}
      </div>
    </details>
  )

  // Nothing to choose between: one quiet line, not a full-width card offering
  // the only option there is plus a warning saying the same thing again.
  if (!usable.length) {
    return (
      <>
        <p className="text-[12.5px]" style={{ color: 'var(--ink-2)' }}>
          No column here can act as{' '}
          {role === 'cluster' ? 'a cell annotation' : `a ${role}`}
          {emptyNote ? ` — ${emptyNote}` : '.'}
        </p>
        {rejectedList}
      </>
    )
  }

  return (
    <>
      {/* auto-fill, not auto-fit: with a single candidate, auto-fit collapses the
          empty tracks and stretches that one card across the whole row. */}
      <div className="grid gap-2" style={{ gridTemplateColumns: 'repeat(auto-fill,minmax(240px,1fr))' }}>
        {allowNone && (
          <button
            onClick={() => onChange(null)}
            className="w-full rounded-xl px-3 py-2.5 text-left"
            style={{
              border: `1.5px solid ${value === null ? 'var(--accent)' : 'var(--line)'}`,
              background: value === null ? 'var(--accent-soft)' : 'var(--surface)',
            }}
          >
            <div className="text-[12.5px] font-semibold">{noneLabel ?? 'Not in this object'}</div>
            <div className="mt-1.5 text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
              Better than picking something that is not it.
            </div>
          </button>
        )}
        {usable.map(c => (
          <Option key={c.name} col={c} chosen={value === c.name} sameAs={dup.get(c.name)}
            onPick={() => onChange(c.name)} />
        ))}
      </div>
      {rejected.length > 0 && (
        <details className="mt-2">
          <summary className="cursor-pointer text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
            {rejected.length} column{rejected.length > 1 ? 's' : ''} cannot be used here
          </summary>
          <div className="mt-1.5 grid gap-1">
            {rejected.map(({ c, why }) => (
              <div key={c.name} className="text-[11.5px]" style={{ color: 'var(--ink-3)' }}>
                <span className="mono">{c.name}</span> — {why}
              </div>
            ))}
          </div>
        </details>
      )}
    </>
  )
}

export default function Review({ scan, choices, setChoices, onBuild, busy, stage }: {
  scan: Scan
  choices: Choices
  setChoices: (c: Choices) => void
  onBuild: () => void
  busy: boolean
  stage: string
}) {
  // Choosing a column for one role releases it from the others, rather than
  // leaving a stale answer that the pickers no longer offer.
  const set = (patch: Partial<Choices>) => {
    const next = { ...choices, ...patch }
    const collides = (a: string | null, b: string | null) => {
      if (!a || !b) return false
      if (a === b) return true
      const ca = column(scan, a), cb = column(scan, b)
      return !!ca && !!cb && samePartition(ca, cb)
    }
    if (collides(next.cluster, next.condition)) next.condition = null
    if (collides(next.cluster, next.sample) || collides(next.condition, next.sample)) {
      next.sample = null
    }
    setChoices(next)
  }

  const clusterCol = column(scan, choices.cluster)
  const sampleCol = column(scan, choices.sample)
  const condCol = column(scan, choices.condition)
  const clash = crossCheck(sampleCol, condCol)
  const reps = replicatesPerCondition(sampleCol, condCol)
  const repValues = Object.values(reps)
  const pbOK = condCol != null && repValues.length > 1 &&
    repValues.every(v => v >= MIN_REPS_PB)

  const hasCounts = scan.matrices.some(m => m.kind === 'counts' || m.kind === 'log-counts')
  const usableMatrix = scan.matrices.some(m => m.kind !== 'scaled')

  return (
    <div className="wrap py-7">
      <div className="flex items-baseline gap-2.5">
        <h1 className="text-[19px] tracking-[-0.02em]">{scan.source}</h1>
        <span className="badge badge-file">{scan.format === 'h5ad' ? 'AnnData' : 'Seurat'}</span>
      </div>
      <p className="sub mt-1.5">
        {scan.nCells.toLocaleString()} cells · {scan.columns.length} metadata columns ·{' '}
        {scan.matrices.length} matri{scan.matrices.length === 1 ? 'x' : 'ces'} ·{' '}
        {scan.embeddings.length} embedding{scan.embeddings.length === 1 ? '' : 's'}
      </p>

      {!usableMatrix && (
        <div className="note mt-4">
          <b>Every matrix here is scaled data.</b> The studio plots expression, not z-scores —
          re-save including the raw counts or the log-normalized matrix.
        </div>
      )}

      <Card
        eyebrow="Step 1 · required"
        title="Which column is the cell type annotation?"
        sub="Ranked by name, but the name is only a hint — check the levels."
      >
        <div className="mt-3.5">
          <Picker scan={scan} role="cluster"
            value={choices.cluster} onChange={v => set({ cluster: v ?? '' })} />
        </div>
      </Card>

      <Card
        eyebrow="Step 2 · optional"
        title="Group by — the experimental condition"
        sub="Treatment, genotype, timepoint. Skip it and the object opens as one group."
      >
        <div className="mt-3.5">
          <Picker scan={scan} role="condition" allowNone
            noneLabel="No condition — one group"
            emptyNote="this object opens as one group"
            exclude={[{ name: choices.cluster, as: 'cell annotation' }]}
            value={choices.condition} onChange={v => set({ condition: v })} />
        </div>
      </Card>

      <Card
        eyebrow="Step 3 · optional"
        title="Sample — animal, donor or run"
        sub="Makes replication visible. Skip it and every cell counts as one sample."
      >
        <div className="mt-3.5">
          <Picker scan={scan} role="sample" allowNone
            noneLabel="No sample — treat as one"
            emptyNote="every cell counts as one sample"
            exclude={[{ name: choices.cluster, as: 'cell annotation' },
                      { name: choices.condition, as: 'condition' }]}
            value={choices.sample} onChange={v => set({ sample: v })} />
        </div>
      </Card>

      {clash && <div className="note mt-4"><b>These two do not nest.</b> {clash}</div>}

      <Card eyebrow="Step 4" title="Embedding and label">
        <Field label="Embedding" hint="What the Cells and Feature plot tabs draw on.">
          <div className="seg">
            {scan.embeddings.map(e => (
              <button key={e.key} className={e.key === choices.embedding ? 'sel' : ''}
                onClick={() => set({ embedding: e.key })}>
                {e.key} <span style={{ opacity: 0.6 }}>{e.nDims}D</span>
              </button>
            ))}
          </div>
          {!scan.embeddings.length && (
            <div className="note mt-2">
              No 2D embedding in this object. Run a UMAP before converting.
            </div>
          )}
          {/pca/i.test(choices.embedding) && (
            <p className="sub mt-2">Principal components only — a coarser picture than a UMAP.</p>
          )}
        </Field>

        <Field label="Label" hint="Shown at the top of the studio.">
          <input className="inp w-full max-w-[420px]" value={choices.label}
            onChange={e => set({ label: e.target.value })} />
        </Field>
      </Card>

      <Card
        eyebrow="Before you convert"
        title="What this bundle will support"
        sub="Follows from the answers above — nothing turns up empty later."
      >
        <div className="scrollx mt-3">
          <table className="t">
            <thead><tr><th>In the studio</th><th>Status</th><th>Because</th></tr></thead>
            <tbody>
              <tr>
                <td className="font-semibold">Clusters, markers, gene search</td>
                <td><span className="badge badge-here">yes</span></td>
                <td style={{ color: 'var(--ink-2)' }}>
                  {clusterCol
                    ? `${clusterCol.levels.length} cell types from ${clusterCol.name}`
                    : 'pick a cell annotation above'}
                </td>
              </tr>
              <tr>
                <td className="font-semibold">Condition comparison</td>
                <td>
                  <span className={`badge ${condCol ? 'badge-here' : 'badge-none'}`}>
                    {condCol ? 'yes' : 'no'}
                  </span>
                </td>
                <td style={{ color: 'var(--ink-2)' }}>
                  {condCol
                    ? `${condCol.levels.length} groups: ${condCol.levels.slice(0, 4).join(', ')}`
                    : 'no condition column chosen'}
                </td>
              </tr>
              <tr>
                <td className="font-semibold">Composition spread</td>
                <td>
                  <span className={`badge ${sampleCol ? 'badge-here' : 'badge-none'}`}>
                    {sampleCol ? 'yes' : 'no'}
                  </span>
                </td>
                <td style={{ color: 'var(--ink-2)' }}>
                  {sampleCol
                    ? `${sampleCol.levels.length} samples`
                    : 'needs a sample column to have anything to vary across'}
                </td>
              </tr>
              <tr>
                <td className="font-semibold">Pseudobulk DESeq2</td>
                <td>
                  <span className={`badge ${pbOK && hasCounts ? 'badge-here' : 'badge-none'}`}>
                    {pbOK && hasCounts ? 'yes' : 'no'}
                  </span>
                </td>
                <td style={{ color: 'var(--ink-2)' }}>
                  {!hasCounts
                    ? 'this object carries no raw counts to sum'
                    : !condCol ? 'needs a condition to compare'
                    : repValues.length < 2 ? 'needs at least two groups'
                    : pbOK ? `${repValues.join(' v ')} samples per group`
                    : `${repValues.join(' v ')} samples per group — fewer than ${MIN_REPS_PB}, so only the per-cell Wilcoxon is offered`}
                </td>
              </tr>
            </tbody>
          </table>
        </div>
        <p className="sub mt-2.5">
          Matrices: {scan.matrices.map(m => `${m.key} (${m.kind})`).join(', ')} — the converter
          takes the log-normalized one or builds it from counts, never <Mono>scaled</Mono>.
        </p>
      </Card>

      <div className="mt-5 flex items-center gap-3">
        <button className="btn btn-primary" disabled={busy || !choices.cluster || !choices.embedding}
          onClick={onBuild}>
          {busy ? (stage || 'Converting…') : 'Convert to bundle.zip'}
        </button>
        {busy && <span className="sub">{stage}</span>}
        {!busy && !choices.cluster && <span className="sub">Pick a cell annotation first.</span>}
      </div>
    </div>
  )
}
