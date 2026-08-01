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
