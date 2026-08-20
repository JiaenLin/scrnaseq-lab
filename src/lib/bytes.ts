// A byte array too big to be one byte array.
//
// V8 caps a single ArrayBuffer. Measured in Chrome 151 on this machine by
// bisection: the largest that allocates is 2.15 GB, and 2.95 GB throws
// `RangeError: Array buffer allocation failed`. The cap is PER BUFFER and not a
// total — 6.14 GB held fine in 512 MB pieces in the same tab — which is the
// whole reason this file can exist.
//
// It was reported as a 2.95 GB Seurat .rds that would not open. The first line
// of the .rds path was `new Uint8Array(await file.arrayBuffer())`, so it threw
// before a byte was parsed, and the decompressed object is larger still. The
// .h5ad path never had this problem because HDF5 is random-access and the lab
// reads it in place, a slice at a time; RDS is a serial format and has to be
// walked from the front, so the bytes have to be somewhere. They are here, in
// pieces small enough to allocate.
//
// One piece is the overwhelmingly common case and it costs nothing: `parts`
// holds a single view, `crosses` is never true, and every read is one DataView
// call on the same object the old code used.

/**
 * 512 MB. Comfortably inside the measured 2.15 GB cap with room for the cap to
 * be lower on another machine, and large enough that a 7 GB object is fourteen
 * pieces rather than a list long enough to matter.
 */
export const CHUNK = 512 * 1024 * 1024

export class Bytes {
  readonly length: number
  private readonly parts: Uint8Array[]
  private readonly views: DataView[]
  /** Staging for a value that straddles two pieces. Reused; never re-entered. */
  private readonly edge = new DataView(new ArrayBuffer(8))

  /** The piece size these parts were packed at. Injectable so a test can use a
   *  chunk of 64 bytes and put a seam through every read this class has. */
  private readonly chunk: number

  constructor(parts: Uint8Array[], chunk: number = CHUNK) {
    this.chunk = chunk
    this.parts = parts.length ? parts : [new Uint8Array(0)]
    this.views = this.parts.map(p => new DataView(p.buffer, p.byteOffset, p.byteLength))
    this.length = this.parts.reduce((n, p) => n + p.byteLength, 0)
    // Every piece but the last must be exactly CHUNK, because `where` divides
    // rather than searching. Building them is this file's job, so this is an
    // invariant check and not an error a user can reach.
    for (let i = 0; i < this.parts.length - 1; i++) {
      if (this.parts[i].byteLength !== chunk) throw new Error('Bytes: ragged chunk')
    }
  }

  /**
   * Which piece a byte is in, and where.
   *
   * Math.floor rather than a shift: `p` runs past 2^31 on the objects this
   * exists for, and JavaScript's bitwise operators are 32-bit signed, so
   * `p >>> 29` silently returns the wrong piece for anything over 2 GB — which
   * is every file that needs this class and no file that tests it by accident.
   */
  private at(p: number): number { return Math.floor(p / this.chunk) }

  /** One byte. There is no numeric indexer on a class, so this is the `b[i]`. */
  byte(p: number): number {
    const c = this.at(p)
    return this.parts[c][p - c * this.chunk]
  }

  getInt32(p: number): number {
    const c = this.at(p)
    const off = p - c * this.chunk
    if (off + 4 <= this.chunk) return this.views[c].getInt32(off, false)
    return this.straddle(p, 4).getInt32(0, false)
  }

  getFloat64(p: number): number {
    const c = this.at(p)
    const off = p - c * this.chunk
    if (off + 8 <= this.chunk) return this.views[c].getFloat64(off, false)
    return this.straddle(p, 8).getFloat64(0, false)
  }

  /** The n bytes at p, copied into `edge`, for a value spanning a seam. */
  private straddle(p: number, n: number): DataView {
    for (let i = 0; i < n; i++) {
      const q = p + i
      const c = this.at(q)
      this.edge.setUint8(i, this.parts[c][q - c * this.chunk])
    }
    return this.edge
  }

  /** The bytes from `a` to `b`, copied only when they cross a seam. */
  subarray(a: number, b: number): Uint8Array {
    // An empty run first. R writes zero-length CHARSXPs — `""` is a perfectly
    // ordinary level name — and when one falls at the very end of the object
    // `at(a)` is one piece past the last, so `parts[c]` is undefined. The flat
    // path never noticed because a Uint8Array subarray past the end is empty
    // rather than an error.
    if (b <= a) return new Uint8Array(0)
    const c = this.at(a)
    const off = a - c * this.chunk
    if (off + (b - a) <= this.chunk) return this.parts[c].subarray(off, off + (b - a))
    const out = new Uint8Array(b - a)
    let at = 0
    let p = a
    while (p < b) {
      const k = this.at(p)
      const from = p - k * this.chunk
      const take = Math.min(this.chunk - from, b - p)
      out.set(this.parts[k].subarray(from, from + take), at)
      at += take
      p += take
      }
    return out
  }

  /**
   * `n` big-endian values from `off` into `out`.
   *
   * The bulk path, and the reason it is not just a loop over getFloat64: the
   * piece is computed once per RUN rather than once per element, so a 100
   * million-element vector pays fourteen divisions instead of a hundred
   * million. On a single-piece Bytes this is exactly the loop it replaces.
   */
  read(out: Float64Array | Int32Array | Float32Array, off: number, n: number, wide: boolean) {
    const size = wide ? 8 : 4
    let i = 0
    let p = off
    while (i < n) {
      const c = this.at(p)
      const from = p - c * this.chunk
      // Whole values left in this piece. A value straddling the seam is read by
      // the slow path and then the run resumes in the next piece.
      const room = Math.floor((this.chunk - from) / size)
      if (room <= 0) {
        out[i++] = wide ? this.getFloat64(p) : this.getInt32(p)
        p += size
        continue
      }
      const take = Math.min(room, n - i)
      const dv = this.views[c]
      if (wide) for (let k = 0; k < take; k++) out[i + k] = dv.getFloat64(from + k * 8, false)
      else for (let k = 0; k < take; k++) out[i + k] = dv.getInt32(from + k * 4, false)
      i += take
      p += take * size
    }
  }
}

/** One buffer, as a Bytes. The path every small file takes. */
export const oneChunk = (b: Uint8Array): Bytes => new Bytes([b])
