import type { ReactNode } from 'react'

export const Mono = ({ children }: { children: ReactNode }) =>
  <code className="mono">{children}</code>

export function Card({ eyebrow, title, sub, right, children }: {
  eyebrow?: string; title?: string; sub?: ReactNode; right?: ReactNode; children: ReactNode
}) {
  return (
    <section className="card">
      {(eyebrow || title || right) && (
        <div className="flex items-start gap-3">
          <div className="min-w-0 flex-1">
            {eyebrow && <div className="eyebrow">{eyebrow}</div>}
            {title && <h2 className="mt-1 text-[15px] font-semibold tracking-[-0.01em]">{title}</h2>}
            {sub && <p className="sub mt-1">{sub}</p>}
          </div>
          {right && <div className="flex-none">{right}</div>}
        </div>
      )}
      {children}
    </section>
  )
}

/** A horizontal bar of level proportions — the shape of a grouping at a glance. */
export function LevelBar({ counts, total }: { counts: number[]; total: number }) {
  const shown = counts.slice(0, 24)
  return (
    <div className="mt-1.5 flex h-[6px] w-full overflow-hidden rounded-full"
      style={{ background: 'var(--sunk)' }}>
      {shown.map((c, i) => (
        <div key={i} style={{
          width: `${(c / total) * 100}%`,
          background: `hsl(${(i * 47) % 360} 58% ${i % 2 ? 62 : 52}%)`,
        }} />
      ))}
    </div>
  )
}

/**
 * Where this page sits between the object you have and the app you want.
 *
 * The lab does exactly one thing, and the fastest way to say so is to show the
 * whole chain with only one link lit up — anyone who lands here from a search
 * result needs to know in one glance that this is a converter, not a viewer.
 */
export function Flow({ at }: { at: 'convert' | 'done' }) {
  const steps: [string, string][] = [
    ['.h5ad / .rds', 'your annotated object'],
    ['scRNA-seq Lab', at === 'done' ? 'converted it' : 'converts it — you are here'],
    ['bundle.zip', at === 'done' ? 'ready to download' : "the studio's input format"],
    ['scRNA-seq Studio', at === 'done' ? 'open it there next' : 'where you explore it'],
  ]
  const lit = at === 'done' ? 3 : 1
  return (
    <ol className="mt-4 flex flex-wrap items-stretch justify-center gap-1.5">
      {steps.map(([name, what], i) => (
        <li key={name} className="flex items-stretch gap-1.5">
          <div
            className="rounded-lg px-2.5 py-1.5 text-left"
            style={{
              background: i === lit ? 'var(--accent-soft)' : 'var(--sunk)',
              border: `1px solid ${i === lit ? 'var(--accent)' : 'transparent'}`,
            }}
          >
            <div className="text-[11.5px] font-semibold leading-tight"
              style={{ color: i === lit ? 'var(--accent-ink)' : 'var(--ink)' }}>{name}</div>
            <div className="text-[10.5px] leading-tight" style={{ color: 'var(--ink-3)' }}>{what}</div>
          </div>
          {i < steps.length - 1 && (
            <span className="self-center text-[11px]" style={{ color: 'var(--ink-3)' }}>→</span>
          )}
        </li>
      ))}
    </ol>
  )
}

export function Field({ label, hint, children }: {
  label: string; hint?: ReactNode; children: ReactNode
}) {
  return (
    <div className="mt-4">
      <div className="eyebrow">{label}</div>
      {hint && <p className="sub mt-1 mb-2">{hint}</p>}
      {children}
    </div>
  )
}
