// Reading R's serialization format, directly.
//
// A Seurat object is an S4 object inside an .rds, and .rds is R's own binary
// serialization — documented, stable since 2015, and perfectly parseable without
// R. The alternative was shipping webR, which is a ~30 MB download and a whole
// R runtime to read one file. This is ~400 lines and reads the file in place.
//
// The one thing that makes it practical on a 288 MB object: numeric vectors are
// NOT materialized while walking. A Seurat object's `scale.data` slot alone can
// be a gigabyte of doubles that nobody asked for. Walking has to visit every
// byte — the format is a stream with no index — but visiting is not allocating,
// so numeric payloads are recorded as (offset, length) and skipped. Only the
// slots the caller actually pulls are ever turned into arrays.

import { Bytes, oneChunk } from './bytes.ts'
import type { GzWindow } from './gzwindow.ts'

/**
 * The offset of a payload that was read past and thrown away.
 *
 * Not -1: that is a legitimate `extra` index (-1 - 0). Any value below every
 * possible one works, and a named constant beats a magic number in a file where
 * a wrong offset is a silently wrong answer.
 */
const DROPPED = -1e9

/**
 * Elements above which a payload is passed over rather than held.
 *
 * 200 million. A Float64Array of that many is 1.6 GB, comfortably inside V8's
 * ~2.15 GB per-array cap with room for `floats()` or `numbers()` to widen a
 * narrower payload on the way out.
 */
const DEFER_ABOVE = 200_000_000

/** SEXP type codes, from R's Rinternals.h. Only the ones that appear here. */
const NILSXP = 0, SYMSXP = 1, LISTSXP = 2, CLOSXP = 3, ENVSXP = 4, PROMSXP = 5,
  LANGSXP = 6, CHARSXP = 9, LGLSXP = 10, INTSXP = 13, REALSXP = 14, CPLXSXP = 15,
  STRSXP = 16, DOTSXP = 17, VECSXP = 19, EXPRSXP = 20, BCODESXP = 21,
  EXTPTRSXP = 22, WEAKREFSXP = 23, RAWSXP = 24, S4SXP = 25,
  BCREPREF = 243, BCREPDEF = 244,
  ALTREP_SXP = 238, ATTRLISTSXP = 239, ATTRLANGSXP = 240, EMPTYENV_SXP = 242,
  BASEENV_SXP = 241, GLOBALENV_SXP = 253, UNBOUNDVALUE_SXP = 252,
  MISSINGARG_SXP = 251, BASENAMESPACE_SXP = 250, NAMESPACESXP = 249,
  PACKAGESXP = 248, PERSISTSXP = 247, CLASSREFSXP = 246, GENERICREFSXP = 245,
  NILVALUE_SXP = 254, REFSXP = 255

/** A numeric vector left in the buffer: read it with `numbers()`. */
export interface RNum {
  /** True when the payload was skipped by the diet; see RdsReader#payload. */
  dropped?: boolean
  /**
   * True when the payload was too large to hold and was passed over instead.
   *
   * `off` still points at it in the decompressed stream, and `n` is still its
   * length, so a caller that can re-open the stream can read it in pieces —
   * which is what a dgCMatrix's `loadGroups` does. What cannot be done is
   * `numbers()`, because the answer would not fit in one array.
   */
  deferred?: boolean
  /** A strided sample of a deferred payload, for classify(). */
  sample?: Float64Array
  t: 'num'
  /** 'double' | 'int' | 'logical' — logical is stored as int32 with NA = INT_MIN. */
  kind: 'double' | 'int' | 'logical'
  n: number
  /** Byte offset of the first element in the decompressed buffer. */
  off: number
  attr?: RAttr
}

export type RValue =
  | { t: 'null' }
  | RNum
  | { t: 'str'; v: (string | null)[]; attr?: RAttr }
  | { t: 'list'; v: RValue[]; attr?: RAttr }
  | { t: 'sym'; v: string }
  | { t: 's4'; attr?: RAttr }
  | { t: 'raw'; n: number; off: number; attr?: RAttr }
  | { t: 'lang'; car: RValue; cdr: RValue; tag: string | null; attr?: RAttr }
  | { t: 'env' }
  | { t: 'other'; code: number; attr?: RAttr }

export type RAttr = Record<string, RValue>

export class RdsError extends Error {}

/** Strip the gzip/bzip2/xz wrapper, or say plainly which one we cannot strip. */
/**
 * What a .rds is compressed with, from its first two bytes.
 *
 * Split out of decompressRds so the streaming path can refuse bzip2 and xz with
 * the same sentence, having read two bytes rather than three gigabytes.
 */
export function rdsCompression(b0: number, b1: number): 'gzip' | 'none' | never {
  if (b0 === 0x1f && b1 === 0x8b) return 'gzip'
  if (b0 === 0x42 && b1 === 0x5a) {
    throw new RdsError(
      'this .rds is bzip2-compressed. Re-save it with saveRDS(obj, "out.rds", compress = "gzip") — gzip is R\'s default, so this file was saved with compress = "bzip2" on purpose.')
  }
  if (b0 === 0xfd && b1 === 0x37) {
    throw new RdsError(
      'this .rds is xz-compressed. Re-save it with saveRDS(obj, "out.rds", compress = "gzip").')
  }
  return 'none'
}

export function decompressRds(
  bytes: Uint8Array, gunzip: (b: Uint8Array) => Uint8Array,
): Uint8Array {
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) return gunzip(bytes)
  if (bytes[0] === 0x42 && bytes[1] === 0x5a) {
    throw new RdsError(
      'this .rds is bzip2-compressed. Re-save it with saveRDS(obj, "out.rds", compress = "gzip") — gzip is R\'s default, so this file was saved with compress = "bzip2" on purpose.')
  }
  if (bytes[0] === 0xfd && bytes[1] === 0x37) {
    throw new RdsError(
      'this .rds is xz-compressed. Re-save it with saveRDS(obj, "out.rds", compress = "gzip").')
  }
  return bytes
}

export class RdsReader {
  readonly buf: Bytes | GzWindow
  private readonly dv: Bytes | GzWindow
  /** Set when the source is a stream, which is when the diet applies. */
  private readonly stream: GzWindow | null
  private p = 0
  /** Symbols and environments, referenced later by index. Order matters. */
  private refs: RValue[] = []
  readonly version: number
  readonly writer: string

  /**
   * A fresh stream over the same file.
   *
   * A deferred payload is read after the walk has gone past it, and the walk's
   * own stream cannot go back. So the caller hands over the means to start
   * another one — see `dgcGroups`, which uses it to read a matrix straight into
   * the per-shard arrays rather than into one that would not fit.
   */
  reopen: (() => GzWindow) | null = null

  /**
   * Elements above which a payload is passed over. DEFER_ABOVE in production;
   * a test sets it to single digits so the streamed path runs on a fixture that
   * fits in a file, instead of only on objects nobody can put in a test.
   */
  private readonly deferAbove: number

  constructor(
    raw: Uint8Array | Bytes | GzWindow,
    reopen?: () => GzWindow,
    deferAbove: number = DEFER_ABOVE,
  ) {
    this.deferAbove = deferAbove
    // A Uint8Array is still accepted and is still the common case — oneChunk
    // wraps it with no copy. Bytes exists for the files that cannot BE a
    // Uint8Array: V8 caps a single buffer at 2.15 GB on this machine, and a
    // 2.95 GB .rds threw on the allocation before it was ever parsed. See
    // lib/bytes.ts.
    this.buf = raw instanceof Uint8Array ? oneChunk(raw) : raw
    this.dv = this.buf
    // `take` is what a stream has and a byte array does not: the chance to
    // materialise or discard a payload at the moment it goes past.
    this.stream = 'take' in this.buf ? (this.buf as GzWindow) : null
    this.reopen = reopen ?? null

    // this.buf, not `raw`: a Bytes is a class and has no numeric indexer, so
    // `raw[0]` on one is undefined and every file large enough to need chunking
    // was rejected here as "does not look like an .rds file" — the one error
    // message guaranteed to send the reader looking in the wrong place.
    const fmt = String.fromCharCode(this.buf.byte(0), this.buf.byte(1))
    if (fmt === 'A\n' || fmt === 'B\n') {
      throw new RdsError(
        `this .rds uses the ${fmt[0] === 'A' ? 'ASCII' : 'native binary'} format. Re-save it with saveRDS(obj, "out.rds") — the default XDR format is what every other tool writes.`)
    }
    if (fmt !== 'X\n') {
      throw new RdsError(
        'this does not look like an .rds file. If it came from save() rather than saveRDS(), it is an .RData holding several objects — load it and re-save the Seurat object with saveRDS().')
    }
    this.p = 2
    this.version = this.int()
    const w = this.int(), _min = this.int()
    void _min
    this.writer = `${w >> 16}.${(w >> 8) & 255}.${w & 255}`
    // Version 3 carries the native encoding as a length-prefixed string.
    if (this.version >= 3) {
      const n = this.int()
      this.p += n
    }
  }

  private int(): number {
    const v = this.dv.getInt32(this.p)
    this.p += 4
    return v
  }

  /** Lengths ≥ 2^31 are written as -1 then a 64-bit pair. */
  private length(): number {
    const n = this.int()
    if (n >= 0) return n
    const hi = this.int(), lo = this.int()
    return hi * 4294967296 + (lo >>> 0)
  }

  private str(n: number): string {
    const s = new TextDecoder('utf-8', { fatal: false })
      .decode(this.buf.subarray(this.p, this.p + n))
    this.p += n
    return s
  }

  /** The entry point: read the whole object graph. */
  read(): RValue { return this.item() }

  private item(): RValue {
    const flags = this.int()
    const type = flags & 255
    const hasAttr = (flags & (1 << 9)) !== 0
    const hasTag = (flags & (1 << 10)) !== 0

    switch (type) {
      case NILVALUE_SXP:
      case NILSXP:
        return { t: 'null' }

      case REFSXP: {
        // The index is packed into the flags unless it is too big to fit.
        const idx = (flags >> 8) || this.int()
        const v = this.refs[idx - 1]
        if (!v) throw new RdsError(`this .rds references object ${idx}, which is not in the file`)
        return v
      }

      case SYMSXP: {
        const name = this.item()
        const sym: RValue = { t: 'sym', v: name.t === 'str' ? (name.v[0] ?? '') : '' }
        this.refs.push(sym)
        return sym
      }

      case PACKAGESXP:
      case NAMESPACESXP: {
        const v = this.item()
        this.refs.push(v)
        return v
      }

      case ENVSXP: {
        // Registered before its contents are read, because the contents may
        // refer back to it. Seurat objects carry environments in `commands`;
        // we walk past them rather than modelling scopes.
        const env: RValue = { t: 'env' }
        this.refs.push(env)
        this.int()          // locked
        this.item()         // enclosure
        this.item()         // frame
        this.item()         // hash table
        this.item()         // attributes
        return env
      }

      case LISTSXP: case LANGSXP: case CLOSXP: case PROMSXP: case DOTSXP:
      case ATTRLISTSXP: case ATTRLANGSXP: {
        // Order is fixed by R: attributes, then tag, then CAR, then CDR.
        const attr = hasAttr || type === ATTRLISTSXP || type === ATTRLANGSXP
          ? this.attrList() : undefined
        let tag: string | null = null
        if (hasTag) {
          const t = this.item()
          tag = t.t === 'sym' ? t.v : null
        }
        /**
         * The slot name, remembered while its value is read.
         *
         * This is the whole of the diet's machinery. R writes a tagged pairlist
         * node as attributes, tag, CAR, CDR in that fixed order, so by the time
         * the value is parsed its name is known — and `scale.data` is the one
         * name worth acting on. Saved and restored because the walk is
         * recursive and a nested list must not inherit its parent's name.
         */
        const outer = this.slotTag
        if (tag) this.slotTag = tag
        const car = this.item()
        this.slotTag = outer
        const cdr = this.item()
        return { t: 'lang', car, cdr, tag, attr }
      }

      case CHARSXP: {
        const n = this.int()
        return { t: 'str', v: [n < 0 ? null : this.str(n)] }
      }

      case LGLSXP: case INTSXP: {
        const n = this.length()
        const off = this.payload(n, false)
        this.p += n * 4
        const v: RNum = { t: 'num', kind: type === LGLSXP ? 'logical' : 'int', n, off }
        if (off === DROPPED) v.dropped = true
        if (this.deferredHere) { v.deferred = true; v.sample = this.deferredHere }
        if (hasAttr) v.attr = this.attrList()
        return v
      }

      case REALSXP: {
        const n = this.length()
        const off = this.payload(n, true)
        this.p += n * 8
        const v: RNum = { t: 'num', kind: 'double', n, off }
        if (off === DROPPED) v.dropped = true
        if (this.deferredHere) { v.deferred = true; v.sample = this.deferredHere }
        if (hasAttr) v.attr = this.attrList()
        return v
      }

      case CPLXSXP: {
        const n = this.length()
        const off = this.p
        this.p += n * 16
        const v: RValue = { t: 'other', code: CPLXSXP }
        if (hasAttr) v.attr = this.attrList()
        void off
        return v
      }

      case RAWSXP: {
        const n = this.length()
        const off = this.p
        this.p += n
        const v: RValue = { t: 'raw', n, off }
        if (hasAttr) v.attr = this.attrList()
        return v
      }

      case STRSXP: {
        const n = this.length()
        const v: (string | null)[] = new Array(n)
        for (let i = 0; i < n; i++) {
          const s = this.item()
          v[i] = s.t === 'str' ? s.v[0] : null
        }
        const out: RValue = { t: 'str', v }
        if (hasAttr) out.attr = this.attrList()
        return out
      }

      case VECSXP: case EXPRSXP: {
        const n = this.length()
        const v: RValue[] = new Array(n)
        for (let i = 0; i < n; i++) v[i] = this.item()
        const out: RValue = { t: 'list', v }
        if (hasAttr) out.attr = this.attrList()
        return out
      }

      case S4SXP: {
        // An S4 object is nothing but its attributes: every slot is one.
        const out: RValue = { t: 's4' }
        if (hasAttr) out.attr = this.attrList()
        return out
      }

      case ALTREP_SXP: {
        const info = this.item()
        const state = this.item()
        const attr = this.item()
        return this.altrep(info, state, attr)
      }

      case EXTPTRSXP: case WEAKREFSXP: {
        const v: RValue = { t: 'other', code: type }
        this.refs.push(v)
        this.item()
        this.item()
        if (hasAttr) v.attr = this.attrList()
        return v
      }

      case BCODESXP: {
        // Every Seurat object carries these: `commands` records the call that
        // produced each step, and a recorded call is a compiled function. We
        // have no use for it, but byte-code has its own sub-format, so it must
        // be walked exactly or every offset after it is wrong.
        this.int()          // how many BCREPDEF slots follow; we hold none
        return this.bc1()
      }

      case EMPTYENV_SXP: case BASEENV_SXP: case GLOBALENV_SXP:
      case BASENAMESPACE_SXP:
        return { t: 'env' }

      case UNBOUNDVALUE_SXP: case MISSINGARG_SXP:
        return { t: 'null' }

      case PERSISTSXP: case CLASSREFSXP: case GENERICREFSXP: {
        const v = this.item()
        this.refs.push(v)
        return v
      }

      default:
        throw new RdsError(
          `this .rds contains an R type this reader does not handle (code ${type}). If the object is a Seurat object, please open an issue with the R version that wrote it.`)
    }
  }

  /**
   * Byte-code. Consumed for its length only.
   *
   * The layout is R's own and does not follow the usual flags-first rule: inside
   * the constant pool a bare type integer precedes each entry, and language
   * nodes read attributes, tag, CAR, CDR in that fixed order. Get one of those
   * ints wrong and the parse does not fail here — it fails hundreds of megabytes
   * later, somewhere unrelated.
   */
  private bc1(): RValue {
    this.item()        // the code vector
    this.bcConsts()
    return { t: 'other', code: BCODESXP }
  }

  private bcConsts(): void {
    const n = this.int()
    for (let i = 0; i < n; i++) {
      const type = this.int()
      if (type === BCODESXP) this.bc1()
      else if (type === LANGSXP || type === LISTSXP || type === BCREPDEF ||
               type === BCREPREF || type === ATTRLANGSXP || type === ATTRLISTSXP) {
        this.bcLang(type)
      } else this.item()
    }
  }

  private bcLang(type: number): void {
    if (type === BCREPREF) { this.int(); return }
    if (type === BCREPDEF || type === LANGSXP || type === LISTSXP ||
        type === ATTRLANGSXP || type === ATTRLISTSXP) {
      let t = type
      if (t === BCREPDEF) { this.int(); t = this.int() }
      if (t === ATTRLANGSXP || t === ATTRLISTSXP) this.item()   // attributes
      this.item()                                               // tag
      this.bcLang(this.int())                                   // CAR
      this.bcLang(this.int())                                   // CDR
      return
    }
    this.item()
  }

  /** ALTREP: a compact representation R expands on access. We expand it too. */
  private altrep(info: RValue, state: RValue, attr: RValue): RValue {
    // info is a pairlist whose CAR is the class symbol.
    let cls = ''
    if (info.t === 'lang' && info.car.t === 'sym') cls = info.car.v
    else if (info.t === 'sym') cls = info.v

    // compact_intseq/realseq store (length, start, step) as three doubles.
    if (cls === 'compact_intseq' || cls === 'compact_realseq') {
      if (state.t === 'num' && state.n >= 3) {
        const [n, start, step] = this.numbers(state) as unknown as number[]
        const out = new Float64Array(n)
        for (let i = 0; i < n; i++) out[i] = start + i * step
        return this.inline(out, cls === 'compact_intseq' ? 'int' : 'double', attr)
      }
    }
    // wrap_* and deferred_string keep the real payload as the head of a list.
    if (state.t === 'lang') return this.withAttr(state.car, attr)
    if (state.t === 'list' && state.v.length > 0) return this.withAttr(state.v[0], attr)
    return this.withAttr(state, attr)
  }

  private refuseDeferred(v: RNum): never {
    throw new RdsError(
      `this slot holds ${v.n.toLocaleString()} values, more than one array can, so it was `
      + 'read past rather than held. It can still be read in parts — that is what the '
      + 'split path does — but not all at once.')
  }

  private refuseDropped(): never {
    throw new RdsError(
      'this slot is scale.data, which was skipped while reading so that an object '
      + 'too large to hold could be opened at all. It holds z-scores, which no view '
      + 'in the studio can use — pick the counts or data matrix instead.')
  }

  /** Which slot the walk is inside, for the diet. See LISTSXP above. */
  private slotTag: string | null = null

  /**
   * What to do with a numeric payload the walk has just reached.
   *
   * On a byte array: nothing. The bytes stay where they are and the offset is
   * recorded, exactly as before — this is the path every file that fits takes,
   * and it allocates nothing.
   *
   * On a stream the bytes exist only now, so the choice has to be made here.
   * Everything is kept except `scale.data`, which is a DENSE genes x cells
   * matrix — no zeros omitted — that the studio cannot use for anything:
   * `chooseMatrices` will not pick a scaled matrix and every view refuses one.
   * On the reported 16.26 GB object it is most of the file. Dropping it is what
   * `DietSeurat(obj, scale.data = FALSE)` does in R, done here so that nobody
   * has to open R to convert their own object.
   */
  /** Set by `payload` when it deferred; consumed by the case that called it. */
  private deferredHere: Float64Array | null = null

  private payload(n: number, wide: boolean): number {
    this.deferredHere = null
    if (!this.stream) return this.p
    const keep = this.slotTag !== 'scale.data'
    // A vector that cannot BE a vector.
    //
    // V8 caps one typed array at 2.15 GB (measured). The reported object holds
    // 618,618,142 non-zeros in each of two matrices: 4.95 GB as doubles, and
    // still 2.47 GB as float32, so neither `x` nor a narrowed copy of it fits
    // in a single array however the reading is done. Streaming solved holding
    // the FILE; it cannot solve holding one slot. Said here, where the number
    // is known, rather than as a bare allocation failure from three frames up.
    /**
     * Too big to hold: passed over, and read later in pieces.
     *
     * V8 caps one typed array at about 2.15 GB, and the object this was built
     * for holds 618,618,142 non-zeros in each of two matrices — 2.47 GB of
     * indices and 4.95 GB of values apiece. No amount of careful unpacking puts
     * that in one array, so it is not put in one: the offset and length are
     * recorded, a few thousand values are sampled on the way past so the scan
     * can still say what KIND of matrix it is, and the payload itself is read
     * later, straight into the per-shard arrays that the split path wants
     * anyway. That is how the .h5ad side has always survived an atlas.
     *
     * The threshold is well under the cap because the values may be widened:
     * `floats()` makes a Float32Array and `numbers()` a Float64Array, so an
     * Int32 payload of 1.2 GB becomes 2.4 GB the moment anybody asks for it as
     * doubles.
     */
    // `p` and `Dim` are never passed over, whatever their size. They are the
    // structure — the column pointers are what ATTRIBUTES a non-zero to a cell,
    // so deferring them would make the shard read impossible rather than merely
    // slow. In practice they are tiny (one int per cell), and the exemption
    // matters only where the threshold is small, which is a test and a
    // pathological object.
    const structural = this.slotTag === 'p' || this.slotTag === 'Dim'
    if (keep && !structural && n > this.deferAbove) {
      const sample = this.stream.sample(this.p, n, wide, 20000)
      this.deferredHere = sample
      return this.p
    }
    const got = this.stream.take(this.p, n, wide, keep)
    if (!got) return DROPPED
    const slot = this.extra.length
    this.extra.push(got)
    return -1 - slot
  }

  /** Expanded ALTREP data does not live in the buffer, so it is parked here. */
  private extra: (Float64Array | Int32Array)[] = []

  private inline(data: Float64Array | Int32Array, kind: 'int' | 'double', attr: RValue): RValue {
    const slot = this.extra.length
    this.extra.push(data)
    const v: RNum = { t: 'num', kind, n: data.length, off: -1 - slot }
    const a = this.asAttr(attr)
    if (a) v.attr = a
    return v
  }

  private withAttr(v: RValue, attr: RValue): RValue {
    const a = this.asAttr(attr)
    if (a && v.t !== 'sym' && v.t !== 'null' && v.t !== 'env') {
      return { ...v, attr: { ...(v as { attr?: RAttr }).attr, ...a } } as RValue
    }
    return v
  }

  private asAttr(v: RValue): RAttr | undefined {
    if (v.t !== 'lang') return undefined
    const out: RAttr = {}
    let node: RValue = v
    while (node.t === 'lang') {
      if (node.tag) out[node.tag] = node.car
      node = node.cdr
    }
    return Object.keys(out).length ? out : undefined
  }

  /** Attributes arrive as a tagged pairlist; a plain object is far easier to use. */
  private attrList(): RAttr | undefined {
    return this.asAttr(this.item())
  }

  /** Materialize a numeric vector. This is where the bytes are finally read. */
  numbers(v: RNum): Float64Array | Int32Array {
    if (v.off === DROPPED) this.refuseDropped()
    if (v.deferred) this.refuseDeferred(v)
    if (v.off < 0) return this.extra[-1 - v.off]
    if (v.kind === 'double') {
      const out = new Float64Array(v.n)
      this.dv.read(out, v.off, v.n, true)
      return out
    }
    const out = new Int32Array(v.n)
    this.dv.read(out, v.off, v.n, false)
    return out
  }

  /** One element, without materializing the vector — for sampling a huge slot. */
  at(v: RNum, i: number): number {
    if (v.off === DROPPED) this.refuseDropped()
    if (v.deferred) {
      // Mapped into the sample, not clamped to it. `stridedSample` walks
      // 0..n-1, so clamping returned the last sampled value for every index
      // past 20 000 — and `classify` would have been judging one number.
      const sm = v.sample!
      return sm[Math.min(sm.length - 1, Math.floor((i / Math.max(1, v.n)) * sm.length))]
    }
    if (v.off < 0) return this.extra[-1 - v.off][i]
    return v.kind === 'double'
      ? this.dv.getFloat64(v.off + i * 8)
      : this.dv.getInt32(v.off + i * 4)
  }

  /** Materialize as float32 — half the memory, and the bundle stores f32 anyway. */
  floats(v: RNum): Float32Array {
    if (v.off === DROPPED) this.refuseDropped()
    if (v.deferred) this.refuseDeferred(v)
    if (v.off < 0) return Float32Array.from(this.extra[-1 - v.off])
    const out = new Float32Array(v.n)
    this.dv.read(out, v.off, v.n, v.kind === 'double')
    return out
  }
}

// ---------------------------------------------------------------------------
// Small helpers for walking what we just parsed.

export const attrOf = (v: RValue | undefined, name: string): RValue | undefined =>
  v && 'attr' in v ? v.attr?.[name] : undefined

export const classOf = (v: RValue | undefined): string[] => {
  const c = attrOf(v, 'class')
  return c?.t === 'str' ? (c.v.filter(Boolean) as string[]) : []
}

export const namesOf = (v: RValue | undefined): string[] => {
  const n = attrOf(v, 'names')
  return n?.t === 'str' ? n.v.map(s => s ?? '') : []
}

/** A named element of an R list, by name. */
export function slot(v: RValue | undefined, name: string): RValue | undefined {
  if (!v) return undefined
  if ('attr' in v && v.attr && name in v.attr) return v.attr[name]
  if (v.t === 'list') {
    const i = namesOf(v).indexOf(name)
    if (i >= 0) return v.v[i]
  }
  return undefined
}

export const strings = (v: RValue | undefined): string[] =>
  v?.t === 'str' ? v.v.map(s => s ?? '') : []

/**
 * An R factor as plain labels. Factors are int codes (1-based, NA = INT_MIN)
 * plus a `levels` attribute — reading the codes as numbers would silently turn
 * a cell type annotation into 1..8.
 */
export function factorLabels(r: RdsReader, v: RValue | undefined): string[] | null {
  if (v?.t !== 'num' || v.kind !== 'int') return null
  const levels = strings(attrOf(v, 'levels'))
  if (!levels.length) return null
  const codes = r.numbers(v) as Int32Array
  const out: string[] = new Array(v.n)
  for (let i = 0; i < v.n; i++) {
    const c = codes[i]
    out[i] = c >= 1 && c <= levels.length ? levels[c - 1] : 'NA'
  }
  return out
}

/** Any column of a data.frame as strings, whatever R type it happens to be. */
export function columnAsStrings(r: RdsReader, v: RValue | undefined): string[] | null {
  if (!v) return null
  const f = factorLabels(r, v)
  if (f) return f
  if (v.t === 'str') return v.v.map(s => s ?? 'NA')
  if (v.t === 'num') {
    const nums = r.numbers(v)
    const out: string[] = new Array(v.n)
    for (let i = 0; i < v.n; i++) {
      const x = nums[i]
      out[i] = v.kind === 'logical'
        ? (x === 0 ? 'FALSE' : x === 1 ? 'TRUE' : 'NA')
        : String(x)
    }
    return out
  }
  return null
}
