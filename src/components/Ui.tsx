import { Fragment, type ReactNode } from 'react'

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
 * The chain, in one line: where this page sits between the object you have and
 * the app you want.
 *
 * The four names carry the whole story on their own — an object goes in, a
 * bundle comes out, the studio opens it — so this used to also caption each
 * step, which said the same thing a second time and cost four lines of height
 * above the button people came to press.
 */
export function Flow({ at }: { at: 'convert' | 'done' }) {
  const steps = ['.h5ad / .rds', 'scRNA-seq Lab', 'bundle.zip', 'scRNA-seq Studio']
  const lit = at === 'done' ? 3 : 1
  return (
    <div className="flex flex-wrap items-center justify-center gap-x-1 gap-y-1">
      {steps.map((name, i) => (
        <Fragment key={name}>
          {i > 0 && <span className="text-[11px]" style={{ color: 'var(--ink-3)' }}>→</span>}
          <span
            className="mono rounded-md px-2 py-[3px] text-[11px]"
            style={i === lit
              ? { background: 'var(--accent-soft)', color: 'var(--accent-ink)', fontWeight: 600 }
              : { background: 'var(--sunk)', color: 'var(--ink-3)' }}
          >{name}</span>
        </Fragment>
      ))}
    </div>
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
