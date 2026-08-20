// RDS reader regressions.
//
// The reader walks a byte stream with no index, so a single miscounted integer
// does not fail where it happened — it fails hundreds of megabytes later, in an
// unrelated slot, or worse, silently returns the wrong numbers. So these tests
// write R's format themselves and read it back: if the two ever disagree about
// where a field ends, it shows up here rather than in someone's data.

import { RdsReader, RdsError, decompressRds, attrOf, classOf, columnAsStrings,
         factorLabels, namesOf, slot, strings } from '../src/lib/rds.ts'
import { Bytes } from '../src/lib/bytes.ts'
import { classify, stridedSample } from '../src/lib/matrix.ts'
import { GzWindow } from '../src/lib/gzwindow.ts'
import { gzipSync } from 'zlib'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const rejects = (name, build, fragment) => {
  try {
    build()
    failed++
    console.log(`  FAIL ${name}\n        it was accepted`)
  } catch (e) {
    const ok = String(e.message).toLowerCase().includes(fragment.toLowerCase())
    if (!ok) failed++
    console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        message was: ${e.message}`}`)
  }
}

// ---------------------------------------------------------------------------
// A minimal writer for R's XDR serialization — the mirror of the reader.

const NILVALUE = 254, SYMSXP = 1, LISTSXP = 2, CHARSXP = 9, INTSXP = 13,
  REALSXP = 14, STRSXP = 16, VECSXP = 19, BCODESXP = 21, S4SXP = 25, ALTREP = 238

class W {
  constructor() { this.b = [] }
  u8(v) { this.b.push(v & 255) }
  int(v) { this.u8(v >> 24); this.u8(v >> 16); this.u8(v >> 8); this.u8(v) }
  dbl(v) {
    const dv = new DataView(new ArrayBuffer(8))
    dv.setFloat64(0, v, false)
    for (let i = 0; i < 8; i++) this.u8(dv.getUint8(i))
  }
  bytes(s) { for (const c of new TextEncoder().encode(s)) this.u8(c) }
  flags(type, { attr = false, tag = false, obj = false } = {}) {
    this.int(type | (obj ? 1 << 8 : 0) | (attr ? 1 << 9 : 0) | (tag ? 1 << 10 : 0))
  }
  charsxp(s) {
    this.flags(CHARSXP)
    const enc = new TextEncoder().encode(s)
    this.int(enc.length)
    for (const c of enc) this.u8(c)
  }
  sym(name) { this.flags(SYMSXP); this.charsxp(name) }
  /** Attributes as R writes them: a tagged pairlist ending in NILVALUE. */
  attrs(pairs) {
    for (const [tag, write] of pairs) {
      this.flags(LISTSXP, { tag: true })
      this.sym(tag)
      write(this)
      // CDR continues the chain; the loop writes the next node's flags next.
    }
    this.flags(NILVALUE)
  }
  real(vals, attrs) {
    this.flags(REALSXP, { attr: !!attrs })
    this.int(vals.length)
    for (const v of vals) this.dbl(v)
    if (attrs) this.attrs(attrs)
  }
  intv(vals, attrs) {
    this.flags(INTSXP, { attr: !!attrs })
    this.int(vals.length)
    for (const v of vals) this.int(v)
    if (attrs) this.attrs(attrs)
  }
  str(vals, attrs) {
    this.flags(STRSXP, { attr: !!attrs })
    this.int(vals.length)
    for (const v of vals) this.charsxp(v)
    if (attrs) this.attrs(attrs)
  }
  list(writers, attrs) {
    this.flags(VECSXP, { attr: !!attrs })
    this.int(writers.length)
    for (const w of writers) w(this)
    if (attrs) this.attrs(attrs)
  }
  s4(attrs) { this.flags(S4SXP, { attr: true, obj: true }); this.attrs(attrs) }
  nil() { this.flags(NILVALUE) }
  done(version = 3) {
    const head = new W()
    head.bytes('X\n')
    head.int(version)
    head.int((4 << 16) | (6 << 8) | 0)   // written by R 4.6.0
    head.int((3 << 16) | (5 << 8) | 0)
    if (version >= 3) { head.int(5); head.bytes('UTF-8') }
    return new Uint8Array([...head.b, ...this.b])
  }
}

const read = (w, version) => {
  const r = new RdsReader(w.done(version))
  return [r, r.read()]
}

console.log('\nHEADER AND VERSIONS')
{
  const mk = v => { const w = new W(); w.real([1, 2, 3]); return read(w, v) }
  for (const v of [2, 3]) {
    const [r, o] = mk(v)
    check(`version ${v} parses`, [r.version, Array.from(r.numbers(o))], [v, [1, 2, 3]])
  }
  check('the writer version is reported', mk(3)[0].writer, '4.6.0')
}
rejects('an ASCII .rds is refused with advice',
  () => new RdsReader(new TextEncoder().encode('A\n0000')), 'ASCII')
rejects('a native-binary .rds is refused',
  () => new RdsReader(new TextEncoder().encode('B\n0000')), 'native binary')
rejects('a non-rds file is refused',
  () => new RdsReader(new TextEncoder().encode('PK\x03\x04....')), 'does not look like an .rds')

console.log('\nCOMPRESSION IS DETECTED BY MAGIC, NOT BY EXTENSION')
check('gzip is unwrapped',
  decompressRds(new Uint8Array([0x1f, 0x8b, 9, 9]), () => new Uint8Array([1, 2])), new Uint8Array([1, 2]))
check('an uncompressed file passes through',
  decompressRds(new Uint8Array([0x58, 0x0a, 1]), () => { throw new Error('should not gunzip') }),
  new Uint8Array([0x58, 0x0a, 1]))
rejects('bzip2 says how to re-save', () => decompressRds(new Uint8Array([0x42, 0x5a, 0]), x => x), 'gzip')
rejects('xz says how to re-save', () => decompressRds(new Uint8Array([0xfd, 0x37, 0]), x => x), 'gzip')

console.log('\nNUMERIC VECTORS ARE READ LAZILY BUT EXACTLY')
{
  const w = new W()
  w.real([0, -1.5, 3.25, 1e300, 0.1])
  const [r, o] = read(w)
  check('values survive the round trip', Array.from(r.numbers(o)), [0, -1.5, 3.25, 1e300, 0.1])
  check('at() agrees with numbers()', [r.at(o, 1), r.at(o, 4)], [-1.5, 0.1])
  check('nothing was materialized while walking', o.t === 'num' && o.off > 0, true)
  check('float32 view rounds, as the bundle will', +r.floats(o)[2].toFixed(3), 3.25)
}
{
  const w = new W()
  w.intv([-2147483648, -1, 0, 7])
  const [r, o] = read(w)
  check('int32 including the NA sentinel', Array.from(r.numbers(o)), [-2147483648, -1, 0, 7])
}

console.log('\nSTRINGS AND NA')
{
  const w = new W()
  w.flags(STRSXP)
  w.int(3)
  w.charsxp('Gene1')
  w.flags(CHARSXP); w.int(-1)          // NA_character_
  w.charsxp('Ünïcøde')
  const [, o] = read(w)
  check('an NA string stays null, not "NA"', o.v, ['Gene1', null, 'Ünïcøde'])
  check('strings() renders NA as text for display', strings(o), ['Gene1', '', 'Ünïcøde'])
}

console.log('\nFACTORS ARE LABELS, NOT CODES')
// The failure this guards against is silent: read the codes as numbers and a
// cell type annotation becomes 1..9, which still renders, still clusters, and
// is wrong everywhere.
{
  const w = new W()
  w.intv([1, 3, 2, 1, -2147483648], [
    ['levels', x => x.str(['B cell', 'NK', 'Monocyte'])],
    ['class', x => x.str(['factor'])],
  ])
  const [r, o] = read(w)
  check('levels come back', strings(attrOf(o, 'levels')), ['B cell', 'NK', 'Monocyte'])
  check('class is readable', classOf(o), ['factor'])
  check('codes become labels, 1-based',
    factorLabels(r, o), ['B cell', 'Monocyte', 'NK', 'B cell', 'NA'])
  check('a plain int vector is not a factor', factorLabels(r, { t: 'num', kind: 'int', n: 1, off: 0 }), null)
  check('any column can be rendered as strings',
    columnAsStrings(r, o), ['B cell', 'Monocyte', 'NK', 'B cell', 'NA'])
}
{
  const w = new W()
  w.real([1.5, 2.5])
  const [r, o] = read(w)
  check('a numeric column renders as strings too', columnAsStrings(r, o), ['1.5', '2.5'])
}

console.log('\nDATA FRAMES AND NAMED LISTS')
{
  const w = new W()
  w.list([
    x => x.real([10, 20]),
    x => x.str(['a', 'b']),
  ], [
    ['names', x => x.str(['nCount_RNA', 'orig.ident'])],
    ['class', x => x.str(['data.frame'])],
    ['row.names', x => x.str(['cell1', 'cell2'])],
  ])
  const [r, o] = read(w)
  check('column names', namesOf(o), ['nCount_RNA', 'orig.ident'])
  check('a column by name', Array.from(r.numbers(slot(o, 'nCount_RNA'))), [10, 20])
  check('a string column by name', strings(slot(o, 'orig.ident')), ['a', 'b'])
  check('an absent name is undefined, not an exception', slot(o, 'nope'), undefined)
  check('class survives', classOf(o), ['data.frame'])
}

console.log('\nS4 OBJECTS ARE THEIR ATTRIBUTES')
{
  const w = new W()
  w.s4([
    ['i', x => x.intv([0, 2])],
    ['p', x => x.intv([0, 1, 2])],
    ['Dim', x => x.intv([3, 2])],
    ['x', x => x.real([1.5, 2.5])],
    ['Dimnames', x => x.list([y => y.str(['G1', 'G2', 'G3']), y => y.str(['c1', 'c2'])])],
    ['class', x => x.str(['dgCMatrix'])],
  ])
  const [r, o] = read(w)
  check('every slot is reachable', classOf(o), ['dgCMatrix'])
  check('Dim', Array.from(r.numbers(slot(o, 'Dim'))), [3, 2])
  check('the column pointers', Array.from(r.numbers(slot(o, 'p'))), [0, 1, 2])
  check('gene names off Dimnames', strings(slot(o, 'Dimnames').v[0]), ['G1', 'G2', 'G3'])
}

console.log('\nBYTE-CODE IS WALKED, NOT CHOKED ON')
// Every Seurat object carries compiled calls in @commands. Byte-code has its
// own sub-format, and mis-walking it corrupts every offset that follows — so
// what matters is that the value written after it still reads correctly.
{
  const w = new W()
  w.flags(VECSXP)
  w.int(2)
  w.flags(BCODESXP)
  w.int(0)                    // no BCREPDEF slots
  w.intv([1, 2, 3])           // the code vector
  w.int(0)                    // no constants
  w.str(['after the bytecode'])
  const [, o] = read(w)
  check('the element after byte-code is intact', strings(o.v[1]), ['after the bytecode'])
}
{
  // Byte-code whose constant pool holds a language node — the branch with the
  // bare type integer before each entry.
  const w = new W()
  w.flags(VECSXP)
  w.int(2)
  w.flags(BCODESXP)
  w.int(0)
  w.intv([9])
  w.int(1)                    // one constant
  w.int(LISTSXP)              // ...of type LISTSXP, read by the bcLang path
  w.nil()                     // tag
  w.int(NILVALUE); w.nil()    // CAR
  w.int(NILVALUE); w.nil()    // CDR
  w.str(['still aligned'])
  const [, o] = read(w)
  check('a language constant does not shift the stream', strings(o.v[1]), ['still aligned'])
}

console.log('\nALTREP IS EXPANDED')
{
  // compact_intseq: R writes row.names this way for any sizeable data.frame.
  const w = new W()
  w.flags(ALTREP)
  w.flags(LISTSXP, { tag: false })
  w.sym('compact_intseq')     // CAR: the class
  w.nil()                     // CDR
  w.real([5, 1, 1])           // state: length, start, step
  w.nil()                     // attributes
  const [r, o] = read(w)
  check('a compact sequence becomes real values', Array.from(r.numbers(o)), [1, 2, 3, 4, 5])
  check('at() works on expanded data too', r.at(o, 3), 4)
}

console.log('\nREFERENCES')
{
  // Symbols are added to the reference table as they appear; a REFSXP that
  // points past the end is corruption, and must say so rather than read junk.
  const w = new W()
  w.flags(REALSXP, { attr: true })
  w.int(1); w.dbl(1)
  w.flags(LISTSXP, { tag: true })
  w.int(255 | (9 << 8))       // REFSXP to entry 9, which does not exist
  w.nil(); w.nil()
  rejects('a dangling reference is refused', () => read(w), 'not in the file')
}

console.log('\nUNKNOWN TYPES FAIL LOUDLY')
{
  const w = new W()
  w.int(200)                  // no such SEXP type
  rejects('an unhandled R type names its code', () => read(w), 'code 200')
}
check('RdsError is the error type', new RdsError('x') instanceof Error, true)


console.log('\nTHE SAME BYTES, IN PIECES')
{
  // Bytes exists because V8 caps a single ArrayBuffer at 2.15 GB (measured by
  // bisection in Chrome 151) while holding 6.14 GB fine in 512 MB pieces, so a
  // 2.95 GB .rds could not be loaded at all. In production the pieces are
  // 512 MB and no test can afford one. So the piece size is injectable and the
  // test uses tiny ones: at 3, 5 and 7 bytes every int32, every float64 and
  // every string in the fixture straddles at least one seam, which is the only
  // part of this class that can be wrong.
  const chop = (b, k) => {
    const out = []
    for (let i = 0; i < b.length; i += k) out.push(b.subarray(i, Math.min(i + k, b.length)))
    // Bytes requires every piece but the last to be exactly k. subarray of the
    // tail is short by construction, which is the shape it wants.
    return out
  }

  // A fixture with one of everything the reader materialises.
  const mk = () => {
    const w = new W()
    w.list([
      x => x.real([1.5, -2.25, 1e300, 0.125, 3]),
      x => x.intv([1, -1, 2147483647, -2147483648, 0]),
      x => x.str(['alpha', 'beta', 'a longer symbol name', '']),
    ])
    return w.done(3)
  }

  const bytes = mk()
  const flat = new RdsReader(bytes)
  const flatObj = flat.read()
  const want = {
    nums: Array.from(flat.numbers(flatObj.v[0])),
    ints: Array.from(flat.numbers(flatObj.v[1])),
    names: strings(flatObj.v[2]),
  }
  check('the flat reader still reads the fixture', want.nums.length, 5)

  let agree = true
  let firstBad = null
  for (const k of [3, 5, 7, 8, 16, 64, 4096]) {
    const r = new RdsReader(new Bytes(chop(bytes, k), k))
    const o = r.read()
    const got = {
      nums: Array.from(r.numbers(o.v[0])),
      ints: Array.from(r.numbers(o.v[1])),
      names: strings(o.v[2]),
    }
    if (JSON.stringify(got) !== JSON.stringify(want)) {
      agree = false
      if (!firstBad) firstBad = { k, got }
    }
  }
  check('every piece size gives the identical object',
    agree ? true : `chunk ${firstBad.k}: ${JSON.stringify(firstBad.got)}`, true)

  // The header itself straddles: version, writer and encoding are read before
  // anything else and sit in the first twenty bytes, so a piece of 3 puts a
  // seam inside each of them.
  const tiny = new RdsReader(new Bytes(chop(bytes, 3), 3))
  check('the header survives a seam through it', [tiny.version, tiny.writer], [3, '4.6.0'])

  // A single element read at a straddle, which is the `at` path rather than the
  // bulk `read` path — used when a slot is too big to materialise.
  const one = new RdsReader(new Bytes(chop(bytes, 5), 5))
  const oo = one.read()
  check('one element at a time agrees too',
    [0, 1, 2, 3, 4].map(i => one.at(oo.v[0], i)), want.nums)
}


console.log('\nWHAT KIND OF MATRIX IT IS SURVIVES THE PIECES')
{
  // The lab never trusts a slot name — it classifies a matrix by its VALUES,
  // sampled with a stride through r.at(). That accessor is the one the chunked
  // reader had to change, and it is the one every downstream decision hangs
  // off: whether pseudobulk is offered at all, whether the studio is handed
  // counts to log-normalize or values already logged, and whether the object is
  // refused as scaled. If chunking shifted a sample by one element the kind
  // could flip, and nothing else in the suite would notice.
  const kinds = {
    counts: [0, 3, 1, 17, 0, 42, 5, 1, 0, 9, 231, 4],
    lognorm: [0, 0.693, 1.0986, 2.3026, 0, 1.6094, 0.4055, 3.1781, 0, 0.9163],
    'log-counts': [0, Math.log1p(3), Math.log1p(17), Math.log1p(1), 0, Math.log1p(42)],
    scaled: [-1.2, 0.4, 2.9, -0.7, 1.1, -2.2, 0.3],
    linear: [0, 120.5, 3300.25, 47.75, 0, 918.5, 12000.125],
  }
  let allAgree = true
  const report = []
  for (const [want, vals] of Object.entries(kinds)) {
    const w = new W()
    w.real(vals)
    const bytes = w.done(3)

    const flat = new RdsReader(bytes)
    const fo = flat.read()
    const flatKind = classify(stridedSample(vals.length, 1000, i => flat.at(fo, i)))
    report.push(`${want}->${flatKind}`)
    if (flatKind !== want) allAgree = false

    for (const k of [3, 5, 7, 8, 64]) {
      const chop = []
      for (let i = 0; i < bytes.length; i += k) {
        chop.push(bytes.subarray(i, Math.min(i + k, bytes.length)))
      }
      const r = new RdsReader(new Bytes(chop, k))
      const o = r.read()
      const kind = classify(stridedSample(vals.length, 1000, i => r.at(o, i)))
      if (kind !== flatKind) { allAgree = false; report.push(`chunk ${k}: ${kind}`) }
    }
  }
  check('every kind is recognised, and identically in pieces',
    allAgree ? true : report.join(' '), true)
  check('and the five kinds are distinct', report.length, 5)
}


console.log('\nSTREAMED, NEVER HELD')
{
  // GzWindow reads an .rds without ever holding it: the parser drives the
  // gunzip and the bytes behind the head are dropped. It exists because a
  // 16.26 GB Seurat object cannot be materialised in a tab that tops out around
  // 15.6 GB, whatever the machine's RAM. Everything below reads the SAME
  // fixture flat and streamed and demands the same answer, at feed sizes small
  // enough that a chunk boundary falls inside every read the class makes.
  const mk = () => {
    const w = new W()
    w.list([
      x => x.real([1.5, -2.25, 1e300, 0.125, 3, 7.75, -0.5]),
      x => x.intv([1, -1, 2147483647, -2147483648, 0, 99]),
      x => x.str(['alpha', 'beta', 'a much longer symbol name than the others', '']),
      x => x.real(Array.from({ length: 5000 }, (_v, i) => i * 0.5)),
    ])
    return w.done(3)
  }
  const bytes = mk()
  const gz = new Uint8Array(gzipSync(Buffer.from(bytes)))
  const src = feed => new GzWindow(
    gz.length,
    (a, b) => gz.subarray(a, b),
    () => {},
    feed)

  const flat = new RdsReader(bytes)
  const fo = flat.read()
  const want = {
    nums: Array.from(flat.numbers(fo.v[0])),
    ints: Array.from(flat.numbers(fo.v[1])),
    names: strings(fo.v[2]),
    big: Array.from(flat.numbers(fo.v[3])).slice(0, 40),
    bigLen: flat.numbers(fo.v[3]).length,
  }

  let agree = true
  let bad = null
  for (const feed of [16, 64, 257, 4096, 1 << 20]) {
    const r = new RdsReader(src(feed))
    const o = r.read()
    const got = {
      nums: Array.from(r.numbers(o.v[0])),
      ints: Array.from(r.numbers(o.v[1])),
      names: strings(o.v[2]),
      big: Array.from(r.numbers(o.v[3])).slice(0, 40),
      bigLen: r.numbers(o.v[3]).length,
    }
    if (JSON.stringify(got) !== JSON.stringify(want)) {
      agree = false
      if (!bad) bad = { feed, got }
    }
  }
  check('streamed gives the same object as held, at every feed size',
    agree ? true : `feed ${bad.feed}`, true)
  check('including a vector long enough to span many pushes', want.bigLen, 5000)

  // The header is read before anything else and is four separate int32s in the
  // first twenty bytes; at feed 16 each arrives in its own push.
  const tiny = new RdsReader(src(16))
  check('the header survives arriving in pieces', [tiny.version, tiny.writer], [3, '4.6.0'])

  // Single-element access, which is the path the matrix classifier uses.
  const one = new RdsReader(src(64))
  const oo = one.read()
  check('one element at a time agrees', [0, 1, 2, 3, 4, 5, 6].map(i => one.at(oo.v[0], i)),
    want.nums)
  check('and the kind is classified the same either way',
    classify(stridedSample(7, 100, i => one.at(oo.v[0], i))),
    classify(stridedSample(7, 100, i => flat.at(fo.v[0], i))))
}

console.log('\nTHE DIET: scale.data GOES PAST WITHOUT BEING KEPT')
{
  // What DietSeurat does in R, done while reading so nobody has to open R.
  // scale.data is a DENSE genes x cells matrix the studio cannot use — every
  // view refuses scaled values — and on a large object it is most of the file.
  const w = new W()
  w.s4([
    ['counts', x => x.real([1, 2, 3, 4])],
    ['scale.data', x => x.real([-1.5, 0.5, 2.5, -0.25, 9, 9, 9, 9])],
    ['data', x => x.real([0.5, 0.7, 1.1, 1.3])],
  ])
  const bytes = w.done(3)
  const gz = new Uint8Array(gzipSync(Buffer.from(bytes)))

  const held = new RdsReader(bytes)
  const ho = held.read()
  const hs = slot(ho, 'scale.data')
  check('held, scale.data is there', Array.from(held.numbers(hs)).length, 8)

  const r = new RdsReader(new GzWindow(gz.length, (a, b) => gz.subarray(a, b), () => {}, 64))
  const o = r.read()
  // The slots either side must be intact — dropping the payload must not shift
  // the walk by a byte, which is the failure this whole file risks.
  check('streamed, counts is intact', Array.from(r.numbers(slot(o, 'counts'))), [1, 2, 3, 4])
  check('and data after it is intact',
    Array.from(r.numbers(slot(o, 'data'))), [0.5, 0.7, 1.1, 1.3])
  const ss = slot(o, 'scale.data')
  check('scale.data is marked dropped rather than wrong', ss.dropped === true, true)
  rejects('and asking for it says why', () => r.numbers(ss), 'scale.data')
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll RDS tests passed\n')
process.exit(failed ? 1 : 0)
