// Reading an .rds without ever holding it.
//
// The .rds path used to decompress the whole object into memory. That works
// until it does not: a reported 2.95 GB Seurat file decompresses to 16.26 GB —
// measured — and a Chrome tab tops out around 15.6 GB, so it failed four
// percent short. Chunking the buffers (lib/bytes.ts) lifted V8's 2.15 GB
// per-buffer cap; this file removes the ceiling from the question entirely.
//
// It rests on one property of the parser, checked before any of this was
// written: RdsReader NEVER rewinds. `this.p` is only ever advanced. So the
// bytes do not have to be addressable, only to arrive in order, and the ones
// already passed can be thrown away.
//
// Which is the other half. Most of a large Seurat object is `scale.data`: a
// DENSE genes x cells matrix with no zeros omitted, which the studio cannot use
// — every view refuses scaled values and `chooseMatrices` will not pick it.
// Skipping its bytes as they stream past is what DietSeurat does in R, done
// here instead, so nobody has to open R to convert their own file.
//
// What is kept is materialised as it streams: a payload the reader will want
// becomes a typed array immediately and the bytes are dropped like the rest.
// Peak memory is the pending queue plus the slots kept, not the whole object.

import { Gunzip } from 'fflate'

/** Compressed input pushed per pump. Output per push is whatever gzip yields. */
const FEED = 2 * 1024 * 1024

/**
 * A forward-only view of a gzip stream, with everything behind the head dropped.
 *
 * Output is held as the queue of chunks gzip produced, never as one buffer.
 * That is deliberate: `Gunzip.push` emits however much it likes and cannot be
 * back-pressured, so a fixed window is a buffer overrun waiting for a file with
 * an unusual compression ratio. A queue absorbs any push and drops whole chunks
 * as the head moves past them.
 */
export class GzWindow {
  /** Absolute offset of the first byte still queued. */
  private start = 0
  /** Pending output, in order, starting at `start`. */
  private q: Uint8Array[] = []
  /** Bytes across `q`. `start + queued` is the end of what can be served. */
  private queued = 0
  private gz: Gunzip
  /** Compressed input consumed, and whether the stream has ended. */
  private at = 0
  private ended = false
  /** Bytes of a payload still to be discarded unread as they arrive. */
  private dropping = 0
  /** Scratch for a value that spans two chunks. */
  private edge = new DataView(new ArrayBuffer(8))

  /** Compressed length, and a synchronous reader for a slice of it. */
  private readonly size: number
  private readonly slice: (a: number, b: number) => Uint8Array
  /** Called as compressed input is consumed, for the stage line. */
  private readonly onRead: (compressed: number, produced: number) => void
  private told = 0
  private readonly feed: number

  constructor(
    size: number,
    slice: (a: number, b: number) => Uint8Array,
    onRead: (compressed: number, produced: number) => void = () => {},
    /**
     * Compressed bytes per push. Injectable so a test can set it to 64 and put
     * a chunk boundary through every read this class makes — at the default a
     * small fixture arrives in one piece and none of the seam handling runs.
     */
    feed: number = FEED,
  ) {
    this.size = size
    this.slice = slice
    this.onRead = onRead
    this.feed = feed
    this.gz = new Gunzip(chunk => this.absorb(chunk))
  }

  /** Decompressed bytes produced so far. Final once the stream has ended. */
  get seen(): number { return this.start + this.queued }

  private absorb(chunk: Uint8Array) {
    if (!chunk.length) return
    if (this.dropping > 0) {
      const n = Math.min(this.dropping, chunk.length)
      this.dropping -= n
      this.start += n
      if (n === chunk.length) return
      chunk = chunk.subarray(n)
    }
    // fflate reuses its output buffer between pushes, so this must be a copy.
    this.q.push(chunk.slice())
    this.queued += chunk.length
  }

  /** Push compressed input until `end` has been produced, or the file ends. */
  private pump(end: number) {
    while (this.start + this.queued < end && this.pumpOnce()) { /* until produced */ }
  }

  /**
   * Move the head to `p`, discarding what it passes.
   *
   * It must PUMP to get there, not simply add to `start`. The first version did
   * the latter, and it was wrong in a way that only shows up when the head
   * jumps forward over bytes that have not been decompressed yet: `start` moved
   * to 1000, the next chunk to arrive still began at absolute 0, and every
   * subsequent read was off by however far the jump had been. Reads that walk
   * forward one field at a time — which is most of the parser — never saw it,
   * because the bytes were always already queued. `subarray(1000, 1016)`
   * returned the first sixteen bytes of the object.
   */
  private advance(p: number) {
    while (this.start < p) {
      if (!this.q.length) {
        if (this.ended) { this.start = p; return }
        this.pumpOnce()
        continue
      }
      const head = this.q[0]
      const need = p - this.start
      if (head.length <= need) {
        this.start += head.length
        this.queued -= head.length
        this.q.shift()
      } else {
        this.q[0] = head.subarray(need)
        this.start += need
        this.queued -= need
      }
    }
  }

  /** One push of compressed input. Returns false at end of file. */
  private pumpOnce(): boolean {
    if (this.ended) return false
    const to = Math.min(this.size, this.at + this.feed)
    const last = to >= this.size
    this.gz.push(this.slice(this.at, to), last)
    this.at = to
    if (last) this.ended = true
    if (this.at - this.told > 100e6 || last) {
      this.told = this.at
      this.onRead(this.at, this.start + this.queued)
    }
    return true
  }

  /** Make [p, p+n) available, dropping everything before p. */
  private ensure(p: number, n: number) {
    if (p < this.start) {
      throw new Error(
        `GzWindow: read at ${p} is behind the head at ${this.start}. The RDS parser `
        + 'only moves forward; this is a bug, not a bad file.')
    }
    this.advance(p)
    this.pump(p + n)
    if (this.start + this.queued < p + n) {
      throw new Error('this .rds ends before the object does — the file is truncated.')
    }
  }

  /** Copy [p, p+n) out of the queue. `p` must already be ensured. */
  private grab(p: number, n: number, out: Uint8Array) {
    let need = n
    let at = 0
    let seek = p - this.start
    for (const c of this.q) {
      if (need === 0) break
      if (seek >= c.length) { seek -= c.length; continue }
      const take = Math.min(c.length - seek, need)
      out.set(c.subarray(seek, seek + take), at)
      at += take
      need -= take
      seek = 0
    }
  }

  byte(p: number): number {
    this.ensure(p, 1)
    const seek = p - this.start
    let s = seek
    for (const c of this.q) {
      if (s < c.length) return c[s]
      s -= c.length
    }
    return 0
  }

  getInt32(p: number): number {
    this.ensure(p, 4)
    this.grab(p, 4, new Uint8Array(this.edge.buffer, 0, 4))
    return this.edge.getInt32(0, false)
  }

  getFloat64(p: number): number {
    this.ensure(p, 8)
    this.grab(p, 8, new Uint8Array(this.edge.buffer, 0, 8))
    return this.edge.getFloat64(0, false)
  }

  subarray(a: number, b: number): Uint8Array {
    if (b <= a) return new Uint8Array(0)
    this.ensure(a, b - a)
    const out = new Uint8Array(b - a)
    this.grab(a, b - a, out)
    return out
  }

  /**
   * `n` big-endian values from `off` — the bulk path, in bites.
   *
   * Never assembles the payload: a 200 million-element vector is converted a
   * few megabytes at a time, so the queue stays small however large the slot.
   */
  read(out: Float64Array | Int32Array | Float32Array, off: number, n: number, wide: boolean) {
    const size = wide ? 8 : 4
    const per = Math.floor((4 * 1024 * 1024) / size)
    const scratch = new Uint8Array(per * size)
    const dv = new DataView(scratch.buffer)
    let i = 0
    while (i < n) {
      const take = Math.min(per, n - i)
      const at = off + i * size
      this.ensure(at, take * size)
      this.grab(at, take * size, scratch)
      if (wide) for (let k = 0; k < take; k++) out[i + k] = dv.getFloat64(k * 8, false)
      else for (let k = 0; k < take; k++) out[i + k] = dv.getInt32(k * 4, false)
      i += take
    }
  }

  /**
   * A numeric payload, materialised now or thrown away now.
   *
   * Called at the moment the parser reaches it, which is the only moment the
   * bytes exist. `keep` false is the diet: the bytes are discarded inside
   * `absorb` as they decompress, so the queue never grows to hold what is being
   * thrown away, and the slot is marked so that asking for it later says why
   * rather than handing back zeros.
   */
  take(p: number, n: number, wide: boolean, keep: boolean): Float64Array | Int32Array | null {
    const bytes = n * (wide ? 8 : 4)
    if (keep) {
      const out = wide ? new Float64Array(n) : new Int32Array(n)
      this.read(out, p, n, wide)
      return out
    }
    this.advance(p)
    // Whatever is already queued can go straight away; the rest is dropped as
    // it arrives, without ever being stored.
    const have = Math.min(this.queued, bytes)
    this.advance(p + have)
    this.dropping = bytes - have
    this.pump(p + bytes)
    this.dropping = 0
    return null
  }
}
