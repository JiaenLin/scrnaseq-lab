// Reading .h5ad, in both layouts AnnData has used.
//
// Modern (AnnData ≥ 0.7): obs and var are groups, one dataset per column, with
// categoricals as a sub-group of categories + codes, and `encoding-type`
// attributes naming everything.
//
// Legacy (< 0.7): obs, var and obsm are each a single compound dataset, and a
// categorical column is bare integer codes whose levels live off in
// uns/<column>_categories. Files in this layout are still everywhere — the
// scanpy PBMC tutorial object is one — and an app that only reads the modern
// one fails on exactly the files people already have.

import {
  categorical, fromFactor, internGeneSet, maybeDiscrete, numeric, pickGeneAlias,
  type Column, type EmbeddingRef, type GeneAlias, type MatrixRef, type RunField, type RunRecord, type RunSection, type Scan,
} from './scan.ts'
import { classify, dropZeros, stridedSample, type CellMajor } from './matrix.ts'
import {
  MAX_PLACE, MAX_VALUE, RECORD_KEYS, SEARCHED, isLine, legacyProvenance, notRecorded, recorded,
  type RecordedId, type Sources,
} from './facts.ts'

/* eslint-disable @typescript-eslint/no-explicit-any */
type H5 = any

/** h5wasm returns attributes through getters; some builds wrap them. */
function attr(o: H5, name: string): unknown {
  const a = o?.attrs?.[name]
  if (a == null) return undefined
  return typeof a === 'object' && a !== null && 'value' in a ? (a as H5).value : a
}

const text = (v: unknown): string =>
  Array.isArray(v) ? String(v[0]) : v instanceof Uint8Array
    ? new TextDecoder().decode(v) : String(v)

/**
 * Which kind of HDF5 node this is, told by what it can do rather than by what
 * it is called.
 *
 * These used to compare `constructor.name` against 'Group' and 'Dataset',
 * which holds only for as long as nothing renames those classes. h5wasm is
 * bundled and minified into this app, so the names surviving is luck rather
 * than a guarantee — and the day it stopped being true the failure would be
 * silent: every group would read as "not a group", and obs would come back
 * empty with no error to explain it.
 */
const isGroup = (o: H5): boolean =>
  !!o && typeof o.keys === 'function' && typeof o.get === 'function'
const isDataset = (o: H5): boolean =>
  !!o && !isGroup(o) && (Array.isArray(o.shape) || typeof o.slice === 'function')

// ---- one file: the input ----------------------------------------------------
//
// HDF5 opens other files by itself: the file an external link names, a
// virtual dataset's source file, a dataset's external raw storage. The node
// build of h5wasm opens them on the host's own filesystem and the browser build
// among the files the tab holds, and through one such link another file's
// content entered a bundle as the object's own, named nowhere: the same input
// gave two bundles when only the file beside it changed. So the converter reads
// one file, its input. `openInput` seals the filesystem HDF5 opens files
// through: while an input it opened is open, no other file can be opened.
// Every lookup of the scan goes through `child`, which names, once, what it
// does not follow, and one note in the bundle lists them.

/** The seal on one h5wasm module's filesystem: how many inputs are open, and whether the door is. */
interface Seal { inputs: number; door: boolean; plugins: Set<string> }
const SEALS = new WeakMap<object, Seal>()

/**
 * Opens the seal refused, counted. HDF5 does not fail every read it could not
 * finish: a dataset whose raw values are stored in another file came back from
 * h5wasm as whatever its buffer held, with no error. A refusal during a read is
 * the one sign that the values are not this file's.
 */
const REFUSED = { opens: 0, filters: 0 }

/** What a dataset is, said, when HDF5 looked for its filter as a plugin: not values in another file. */
const FILTERED = 'a dataset compressed with a filter this reader does not have, LZF for one: write the file again with gzip compression'

/**
 * This filesystem's seal, put on at its first use: from then on, while an input
 * opened by `openInput` is open, every open is refused (EACCES) except one made
 * through the door, which only the converter's own operations use.
 */
function sealOf(FS: H5): Seal {
  let seal = SEALS.get(FS)
  if (!seal) {
    const s: Seal = { inputs: 0, door: false, plugins: new Set() }
    const open = FS.open
    FS.open = function (this: unknown, ...args: unknown[]) {
      if (s.inputs > 0 && !s.door) {
        REFUSED.opens++
        // HDF5 looks for a filter it lacks in its own plugin directories, through this open.
        const at = String(args[0])
        if ([...s.plugins].some(d => at === d || at.startsWith(`${d}/`))) REFUSED.filters++
        throw new FS.ErrnoError(2)
      }
      return open.apply(this, args)
    }
    SEALS.set(FS, seal = s)
  }
  return seal
}

/** One of the converter's own operations on the filesystem: the door is open only while it runs. */
function throughDoor<T>(FS: H5, fn: () => T): T {
  const s = sealOf(FS)
  const was = s.door
  s.door = true
  try { return fn() } finally { s.door = was }
}

/** What the scan met and did not follow, per open input (by file id): path -> what it is. */
const UNFOLLOWED = new Map<string, Map<string, string>>()

/**
 * Open `path` read-only as the converter's input: the one way in the scan
 * accepts. Until the handle is closed, HDF5 can open no other file.
 */
export function openInput(h5wasm: H5Module, lib: H5Lib, path: string): H5 {
  guardReads(lib)
  const FS = lib.FS
  const s = sealOf(FS)
  for (const d of lib.get_plugin_search_paths?.() ?? []) s.plugins.add(String(d))
  const f = throughDoor(FS, () => new h5wasm.File(path, 'r'))
  if (!(Number(f?.file_id) >= 0)) return f
  const id = String(f.file_id)
  s.inputs++
  UNFOLLOWED.set(id, new Map())
  const close = f.close
  let closed = false
  f.close = function (this: unknown, ...args: unknown[]) {
    if (!closed) { closed = true; s.inputs--; UNFOLLOWED.delete(id) }
    return close.apply(f, args)
  }
  return f
}

/** What a dataset is when its values lie in another file, as its mapping names it; null for any other. */
function elsewhere(d: H5): string | null {
  try {
    const files = [...new Set(((d.metadata?.virtual_sources ?? []) as { file_name: string }[])
      .map(v => v.file_name).filter(name => name !== '.'))]
    return files.length ? 'a virtual dataset whose values lie in another file' : null
  } catch {
    return null
  }
}

/**
 * Read one value of an array dataset, and say whether that read tried to open
 * another file: external raw storage, or a virtual source no mapping names
 * ('file'); or only HDF5's own plugin directories, looking for a filter it
 * lacks ('filter'). One element is enough, since the first lies where the rest do.
 */
function storedElsewhere(d: H5): 'file' | 'filter' | null {
  const shape = (d.shape as number[] | null) ?? []
  if (!shape.length || shape.some(x => !(x > 0))) return null
  // h5wasm's slice of a variable-length type ran out of its memory, and the
  // scan reads those whole anyway. Raw values stored in another file have a
  // fixed size.
  try {
    if (d.metadata?.vlen || d.metadata?.vlen_type) return null
  } catch {
    return null
  }
  const before = REFUSED.opens, filters = REFUSED.filters
  try {
    d.slice(shape.map(() => [0, 1]))
  } catch {
    // A read that fails for another reason is the reader's to report, later.
  }
  const refused = REFUSED.opens - before
  return !refused ? null : refused === REFUSED.filters - filters ? 'filter' : 'file'
}

/** A scalar's value, read only from this file. */
function ownValue(d: H5): unknown {
  const before = REFUSED.opens
  const v = d.value
  if (REFUSED.opens !== before) throw new Error('its value is stored in another file, which the converter does not read')
  return v
}

/**
 * The object at `name` below `g`, never reached through a link out of the file.
 *
 * Looked up one name at a time: for a path that runs through a link it cannot
 * follow, h5wasm hands back an object shaped like a group whose every read
 * fails, while one name at a time it hands back the link itself. An external
 * link (whatever file it names: under the seal none can be opened), a soft link
 * that reaches nothing in this file, and a virtual dataset whose values lie in
 * another file are each remembered once, by their own place in this file, and
 * are absent to the scan. What a link points to is never carried: its target
 * text can be a path on the machine that wrote the file.
 */
function child(g: H5, name: string): H5 | null {
  let at: H5 = g
  for (const part of name.split('/')) {
    if (!part) continue
    if (!isGroup(at)) return null
    const o = at.get(part)
    const says = o?.type === 'ExternalLink' ? 'an external link, which points outside this file'
      : o?.type === 'BrokenSoftLink' ? 'a soft link that reaches nothing in this file'
        : isDataset(o) ? elsewhere(o) ?? { file: 'a dataset whose values are stored in another file', filter: FILTERED, none: null }[storedElsewhere(o) ?? 'none'] : null
    if (says) {
      UNFOLLOWED.get(String(at.file_id))?.set(`${at.path}/${part}`.replace(/^\/+/, ''), says)
      return null
    }
    if (!o) return null
    at = o
  }
  return at
}

/** What `part` below `g` is, when the scan did not follow it; null otherwise. */
function notFollowed(g: H5, part: string): string | null {
  return UNFOLLOWED.get(String(g?.file_id))?.get(`${g.path}/${part}`.replace(/^\/+/, '')) ?? null
}

/** An obs column the scan could not read, and why: absent, and named. */
interface Lost { name: string; why: string }

/** Why a part the scan needs is missing, when it is something the converter did not follow. */
function unfollowed(f: H5, path: string): string | null {
  const says = UNFOLLOWED.get(String(f.file_id))?.get(path)
  return says === FILTERED ? `/${path} is ${says}` : says ? `/${path} is ${says}, and the converter reads only the file it is given` : null
}

/** BigInt64Array and friends do not coerce, so every read goes through here. */
function toF32(v: unknown): Float32Array {
  if (v instanceof Float32Array) return v
  const a = v as ArrayLike<number | bigint>
  const out = new Float32Array(a.length)
  for (let i = 0; i < a.length; i++) out[i] = Number(a[i])
  return out
}

function toI32(v: unknown): Int32Array {
  if (v instanceof Int32Array) return v
  const a = v as ArrayLike<number | bigint>
  const out = new Int32Array(a.length)
  for (let i = 0; i < a.length; i++) out[i] = Number(a[i])
  return out
}

const toStrings = (v: unknown): string[] =>
  (v as unknown[]).map(x => (x instanceof Uint8Array ? new TextDecoder().decode(x) : String(x)))

export class H5adError extends Error {}

// ---- one read: every value the converter takes from the file ----------------
//
// h5wasm reads a dataset's or an attribute's values through one call into
// HDF5, and drops the status that call returns. A read HDF5 could not finish (a
// damaged compressed chunk, a filter or a conversion it lacks, values in a file
// the seal refused) came back as whatever its buffer held: an embedding offered
// with half its rows garbage, a matrix used as the expression, sample codes
// that moved cells between samples, all deterministic and named nowhere.
// `guardReads` wraps that call, once per module and the way the seal wraps the
// filesystem's open, so every read the converter makes goes through it: a read
// HDF5 reports as failed throws `UnreadableError`, naming what it read, and its
// buffer is never handed on as values. Each use then decides, by the rule the
// README states: an optional part is named in a note and is absent; a part the
// bundle depends on stops the conversion, by name.

/** A read HDF5 reported as failed: none of its values are used. */
export class UnreadableError extends H5adError {
  /** What the failed read was of: a dataset's path, or an attribute and its object's. */
  readonly what: string
  constructor(what: string) {
    super(`HDF5 could not read ${what} (a damaged chunk, or a filter or a type it cannot decode), so none of its values are used`)
    this.what = what
  }
}

/** What `h5wasm.ready` resolves to: the HDF5 module, with its filesystem and its read calls. */
export interface H5Lib {
  FS: H5
  get_dataset_data: (...a: any[]) => number
  get_attribute_data: (...a: any[]) => number
  reclaim_vlen_memory: (...a: any[]) => number
  get_plugin_search_paths?: () => string[]
}

const GUARDED = new WeakSet<object>()

const pathOf = (path: unknown) => String(path).replace(/^\/+/, '') || '/'

/** Every read of this module checked, from now on: see above. */
function guardReads(lib: H5Lib): void {
  if (GUARDED.has(lib)) return
  if (typeof lib?.get_dataset_data !== 'function' || typeof lib.get_attribute_data !== 'function'
    || typeof lib.reclaim_vlen_memory !== 'function') {
    throw new H5adError('openInput takes the module h5wasm.ready resolves to, so that every read the converter makes is checked')
  }
  const { get_dataset_data: dataset, get_attribute_data: attribute, reclaim_vlen_memory: reclaim } = lib
  // The buffer of the read that just failed. Its variable-length entries are
  // whatever HDF5 left there, so it is not freed as strings: that would free
  // memory nobody allocated. It is leaked instead.
  const failed = new Set<unknown>()
  const checked = (call: (...a: any[]) => number, ptr: number, what: (a: any[]) => string) =>
    function (this: unknown, ...a: any[]): number {
      failed.clear()
      let status: number
      try {
        status = call.apply(this, a)
      } catch {
        status = -1
      }
      if (!(Number(status) >= 0)) {
        failed.add(a[ptr])
        throw new UnreadableError(what(a))
      }
      return status
    }
  lib.get_dataset_data = checked(dataset, 5, a => pathOf(a[1]))
  lib.get_attribute_data = checked(attribute, 3, a => `the attribute ${String(a[2])} of ${pathOf(a[1])}`)
  lib.reclaim_vlen_memory = function (this: unknown, ...a: any[]): number {
    if (failed.delete(a[3])) return 0
    return reclaim.apply(this, a)
  }
  GUARDED.add(lib)
}

/** HDF5's classes of stored type, as a note names them. */
const STORED: Record<number, string> = {
  2: 'a time type', 3: 'strings', 4: 'a bitfield', 5: 'opaque bytes', 6: 'a compound type', 7: 'references',
  8: 'an enumeration', 9: 'variable-length sequences', 10: 'an array type',
}

/**
 * Why a dataset's values are not numbers the converter can use, from its
 * stored type alone; null when they are. Never reinterpreted: h5wasm hands an
 * opaque or unknown type back as its raw bytes, which read as numbers.
 */
function notNumbers(d: H5, integers = false): string | null {
  let m: H5
  try {
    m = d.metadata
  } catch {
    return 'a dataset whose stored type cannot be read'
  }
  if (!m?.shape) return 'a dataset with no values (an empty dataspace)'
  if (m.type === 0) return null
  if (m.type === 1) {
    if (integers) return 'stored as floats, not integers'
    return m.size === 4 || m.size === 8 ? null : `stored as ${m.size * 8}-bit floats, which this reader does not read`
  }
  return `stored as ${STORED[m.type as number] ?? 'a type HDF5 does not name'}, not ${integers ? 'integers' : 'numbers'}`
}

/** A read the bundle depends on: a failed one stops the conversion, saying what it was for. */
function needed<T>(what: string, read: () => T): T {
  const filters = REFUSED.filters
  try {
    return read()
  } catch (e) {
    // HDF5 looked for a filter it lacks during the read: say which way forward there is.
    const why = e instanceof UnreadableError && REFUSED.filters !== filters ? `; ${e.what} is ${FILTERED}` : ''
    if (e instanceof UnreadableError) throw new H5adError(`${what} cannot be read, so the conversion stops: ${e.message}${why}`)
    throw e
  }
}

/** The longest name from the input that a note or a list carries whole, in characters (code points). */
export const NAME_MAX = 200

/** A name read from the input, as a note or a list carries it: whole up to NAME_MAX characters, else cut there and marked. */
export function shown(name: string): string {
  if (!longer(name, NAME_MAX)) return name
  const points = Array.from(name)
  return `${points.slice(0, NAME_MAX).join('')}… (cut: a name of ${points.length.toLocaleString('en-US')} characters)`
}

/**
 * A string column in anndata's nullable encoding: a group of `values` + `mask`.
 *
 * This is what anndata (0.11 on) writes for pandas' own string dtype, and on
 * pandas 3 that dtype is the default — so an object written with default
 * settings stores its obs and var names this way, and every string column it
 * did not turn into a categorical. The family's writers opt out of it; anything
 * written with the defaults does not, and the reader used to fail on the first
 * one with a property read on undefined. A masked entry is missing: its slot in
 * `values` is a placeholder, never the value.
 *
 * Null for anything else, so callers fall through to the encodings they knew.
 */
function nullableStrings(o: H5): { values: string[]; mask: Uint8Array } | null {
  if (!isGroup(o) || text(attr(o, 'encoding-type') ?? '') !== 'nullable-string-array') return null
  const v = child(o, 'values')
  if (!isDataset(v)) return null
  const values = toStrings(v.value)
  const mask = new Uint8Array(values.length)
  const m = child(o, 'mask')
  if (isDataset(m)) {
    const raw = m.value as ArrayLike<number | bigint | boolean>
    for (let i = 0; i < mask.length; i++) mask[i] = raw[i] ? 1 : 0
  }
  return { values, mask }
}

/** Masked entries as empty strings — where every caller treats an empty name as no name. */
const unmasked = (s: { values: string[]; mask: Uint8Array }): string[] =>
  s.values.map((v, i) => (s.mask[i] ? '' : v))

/** A nullable string column as a grouping: masked cells are NA, never a level of their own. */
function maskedColumn(name: string, s: { values: string[]; mask: Uint8Array }): Column {
  const index = new Map<string, number>()
  const levels: string[] = []
  const codes = new Int32Array(s.values.length)
  for (let i = 0; i < s.values.length; i++) {
    if (s.mask[i]) { codes[i] = -1; continue }
    let at = index.get(s.values[i])
    if (at === undefined) {
      at = levels.length
      levels.push(s.values[i])
      index.set(s.values[i], at)
    }
    codes[i] = at
  }
  return fromFactor(name, codes, levels)
}

/**
 * The largest .h5ad the in-tab reader can accept.
 *
 * 2^31 - 1: h5wasm's emulated filesystem stores lengths as signed 32-bit.
 * Verified empirically against this build — a 2047 MB write succeeds and a
 * 2048 MB write throws.
 */
export const MAX_H5AD_BYTES = 2 ** 31 - 1

/** The compound-dataset layout: one array of rows, each row an array of fields. */
interface Compound { names: string[]; rows: unknown[][]; members: H5[] }

function compound(ds: H5): Compound | null {
  const members = ds?.metadata?.compound_type?.members
  if (!members) return null
  return { names: members.map((m: H5) => m.name as string), rows: ds.value as unknown[][], members }
}

const fieldOf = (c: Compound, name: string): unknown[] | null => {
  const i = c.names.indexOf(name)
  return i < 0 ? null : c.rows.map(r => r[i])
}

/**
 * Read obs into columns, whichever layout it is in.
 *
 * `uns` is passed because in the legacy layout that is the only place a
 * categorical column's level names exist — the column itself is bare integers,
 * and reading it without them turns "B cells" into 2.
 */
function readObs(obs: H5, uns: H5, nCells: number, note: (s: string) => void, lost: Lost[]): Column[] {
  const out: Column[] = []
  const unread = (name: string, why: string) => {
    note(`obs/${shown(name)} could not be read (${why}); skipped`)
    lost.push({ name, why })
  }

  const comp = needed('the obs table', () => compound(obs))
  if (comp) {
    for (let i = 0; i < comp.names.length; i++) {
      const name = comp.names[i]
      if (name === 'index' || name === '_index') continue
      const raw = comp.rows.map(r => r[i])
      const catsDs = uns ? child(uns, `${name}_categories`) : null
      let cats: string[] | null
      try {
        cats = isDataset(catsDs) ? toStrings(catsDs.value) : null
      } catch (e) {
        // Without its levels the column is bare integers: left out, never read as numbers.
        unread(name, (e as Error).message)
        continue
      }
      if (cats) {
        out.push(fromFactor(name, toI32(raw.map(Number)), cats))
      } else if (typeof raw[0] === 'string') {
        out.push(categorical(name, raw.map(String)))
      } else {
        out.push(maybeDiscrete(numeric(name, toF32(raw.map(Number))), nCells))
      }
    }
    return out
  }

  if (!isGroup(obs)) throw new H5adError('this file has no obs table, so there is nothing to group cells by')

  const indexName = text(attr(obs, '_index') ?? '_index')
  const order = attr(obs, 'column-order')
  const keys: string[] = order ? toStrings(order) : obs.keys().filter((k: string) => k !== indexName)

  for (const key of keys) {
    if (key === indexName) continue
    const o = child(obs, key)
    if (!o) {
      const says = notFollowed(obs, key)
      if (says) lost.push({ name: key, why: `obs/${shown(key)} is ${says}` })
      continue
    }
    try {
      if (isGroup(o)) {
        const strings = nullableStrings(o)
        if (strings) {
          if (strings.values.length !== nCells) {
            note(`obs/${key} has ${strings.values.length} values for ${nCells} cells; skipped`)
          } else {
            out.push(maskedColumn(key, strings))
          }
          continue
        }
        // A modern categorical: categories + codes, with -1 for NA.
        const cats = child(o, 'categories')
        const codes = child(o, 'codes')
        if (cats && codes) {
          const bad = notNumbers(codes, true)
          if (bad) throw new Error(`its codes are ${bad}`)
          out.push(fromFactor(key, toI32(codes.value), toStrings(cats.value)))
          continue
        }
        const gone = ['categories', 'codes'].find(p => notFollowed(o, p))
        if (gone) throw new Error(`obs/${shown(key)}/${gone} is ${notFollowed(o, gone)}`)
        note(`obs/${key} is a group this reader does not understand; skipped`)
        continue
      }
      if (o.metadata?.type !== 3 && o.metadata?.type !== 8) {
        const bad = notNumbers(o)
        if (bad) throw new Error(`it is ${bad}`)
      }
      const v = o.value
      if (v?.length !== nCells) {
        note(`obs/${key} has ${v?.length ?? 0} values for ${nCells} cells; skipped`)
        continue
      }
      if (typeof v[0] === 'string' || v[0] instanceof Uint8Array) {
        out.push(categorical(key, toStrings(v)))
      } else {
        out.push(maybeDiscrete(numeric(key, toF32(v)), nCells))
      }
    } catch (e) {
      unread(key, (e as Error).message)
    }
  }
  return out
}

/** Gene names from a var table, either layout. */
function readVar(v: H5, label: string): string[] {
  const comp = compound(v)
  if (comp) {
    const idx = fieldOf(comp, 'index') ?? fieldOf(comp, '_index')
    if (idx) return toStrings(idx)
  }
  if (isGroup(v)) {
    const name = text(attr(v, '_index') ?? '_index')
    const ds = child(v, name) ?? child(v, '_index') ?? child(v, 'index')
    const strings = nullableStrings(ds)
    if (strings) return unmasked(strings)
    if (isDataset(ds)) return toStrings(ds.value)
    if (ds) {
      throw new H5adError(`${label}/${name} is stored as "${text(attr(ds, 'encoding-type') ?? 'an unnamed group')}", which this reader cannot read as gene names`)
    }
  }
  throw new H5adError(`${label} has no gene names`)
}

/**
 * Every string-valued var column except the index.
 *
 * Read whole rather than sampled: a var table is one row per gene, so even an
 * atlas's is a few megabytes, and the alias is written into the bundle in full
 * anyway. Categoricals are expanded here — on this atlas the symbols are stored
 * as 31 017 categories over 31 053 codes, and reading only the categories would
 * misalign every gene after the first duplicate.
 */
function readVarStringColumns(v: H5, nGenes: number, note: (s: string) => void): {
  name: string; values: string[]
}[] {
  const out: { name: string; values: string[] }[] = []
  const comp = compound(v)
  if (comp) {
    for (const name of comp.names) {
      if (name === 'index' || name === '_index') continue
      const f = fieldOf(comp, name)
      if (!f || f.length !== nGenes) continue
      if (typeof f[0] !== 'string' && !(f[0] instanceof Uint8Array)) continue
      out.push({ name, values: toStrings(f) })
    }
    return out
  }
  if (!isGroup(v)) return out

  const indexName = text(attr(v, '_index') ?? '_index')
  for (const key of v.keys()) {
    if (key === indexName || key === '_index' || key === 'index') continue
    try {
      const o = child(v, key)
      if (!o) continue
      const strings = nullableStrings(o)
      if (strings) {
        if (strings.values.length === nGenes) out.push({ name: key, values: unmasked(strings) })
        continue
      }
      if (isGroup(o)) {
        const cats = child(o, 'categories')
        const codes = child(o, 'codes')
        if (!cats || !codes) continue
        const levels = toStrings(cats.value)
        const c = toI32(codes.value)
        if (c.length !== nGenes) continue
        out.push({ name: key, values: Array.from(c, i => (i >= 0 ? levels[i] ?? '' : '')) })
        continue
      }
      const val = o.value
      if (!val || val.length !== nGenes) continue
      if (typeof val[0] !== 'string' && !(val[0] instanceof Uint8Array)) continue
      out.push({ name: key, values: toStrings(val) })
    } catch (e) {
      note(`var/${shown(key)} could not be read (${(e as Error).message}); skipped`)
    }
  }
  return out
}

/** Symbols beside accessions, or the reverse — whichever the var table adds. */
function readVarAlias(
  v: H5, genes: string[], note: (s: string) => void,
): GeneAlias | null {
  try {
    return pickGeneAlias(genes, readVarStringColumns(v, genes.length, note))
  } catch (e) {
    note(`the var table's other gene names could not be read (${(e as Error).message})`)
    return null
  }
}

/** A matrix, dense or sparse, in either the modern or the h5sparse layout. */
function readMatrix(m: H5, nCells: number, nGenes: number, key: string): CellMajor {
  if (isDataset(m)) {
    // Dense, cells × genes. Materialized whole, then thinned to nonzeros.
    const flat = toF32(m.value)
    const shape = m.shape as number[]
    if (shape[0] !== nCells || shape[1] !== nGenes) {
      throw new H5adError(`${key} is ${shape.join('×')}, expected ${nCells}×${nGenes}`)
    }
    let nnz = 0
    for (let i = 0; i < flat.length; i++) if (flat[i] !== 0) nnz++
    const indptr = new Int32Array(nCells + 1)
    const indices = new Int32Array(nnz)
    const data = new Float32Array(nnz)
    let at = 0
    for (let c = 0; c < nCells; c++) {
      indptr[c] = at
      const base = c * nGenes
      for (let g = 0; g < nGenes; g++) {
        const v = flat[base + g]
        if (v === 0) continue
        indices[at] = g
        data[at] = v
        at++
      }
    }
    indptr[nCells] = at
    return { nCells, nGenes, indptr, indices, data }
  }

  const enc = text(attr(m, 'encoding-type') ?? attr(m, 'h5sparse_format') ?? '')
  const data = toF32(child(m, 'data').value)
  const indices = toI32(child(m, 'indices').value)
  const indptr = toI32(child(m, 'indptr').value)

  if (enc.includes('csc')) {
    // Stored gene-major already; flip it to cell-major so one path follows.
    if (indptr.length !== nGenes + 1) {
      throw new H5adError(`${key} says CSC but its indptr is ${indptr.length}, not ${nGenes + 1}`)
    }
    const counts = new Int32Array(nCells + 1)
    for (let k = 0; k < indices.length; k++) counts[indices[k] + 1]++
    for (let c = 0; c < nCells; c++) counts[c + 1] += counts[c]
    const cursor = counts.slice(0, nCells)
    const ci = new Int32Array(indices.length)
    const cd = new Float32Array(data.length)
    for (let g = 0; g < nGenes; g++) {
      for (let k = indptr[g]; k < indptr[g + 1]; k++) {
        const at = cursor[indices[k]]++
        ci[at] = g
        cd[at] = data[k]
      }
    }
    return dropZeros({ nCells, nGenes, indptr: counts, indices: ci, data: cd })
  }

  if (indptr.length !== nCells + 1) {
    throw new H5adError(`${key} has an indptr of ${indptr.length} for ${nCells} cells — is it transposed?`)
  }
  return dropZeros({ nCells, nGenes, indptr, indices, data })
}


/**
 * How much of the gene list a CSR matrix ever touches.
 *
 * Read the same way values are sampled — a few contiguous blocks spread across
 * the index vector — because a single block covers only a handful of cells and
 * would understate a matrix whose detected genes vary from cell to cell. That
 * is the whole distinction being measured: a real expression matrix touches
 * different genes in different cells and adds up to most of the list, while a
 * model fitted on selected features touches the same few in every cell however
 * far you read.
 */
function geneCoverage(m: H5, nGenes: number): number {
  const idx = child(m, 'indices')
  const n = (idx.shape as number[])[0]
  if (!n || !nGenes) return 1
  const blocks = 24
  const per = Math.max(1, Math.min(200_000, Math.floor(n / blocks)))
  const seen = new Set<number>()
  for (let b = 0; b < blocks; b++) {
    const at = Math.max(0, Math.min(n - per, Math.floor((n / blocks) * b)))
    const s = idx.slice([[at, at + per]]) as ArrayLike<number | bigint>
    for (let i = 0; i < s.length; i++) seen.add(Number(s[i]))
    if (seen.size >= nGenes) break
  }
  return seen.size / nGenes
}

/** Row offsets of a CSR matrix, as counts per cell. */
function rowLengths(m: H5, nCells: number): Int32Array {
  const indptr = toI32(child(m, 'indptr').value)
  const out = new Int32Array(nCells)
  for (let c = 0; c < nCells; c++) out[c] = indptr[c + 1] - indptr[c]
  return out
}

/** Values read per slice. Large enough that chunk decompression amortizes. */
const SCAN_BLOCK = 8_000_000

/**
 * Read disjoint cell subsets from a CSR matrix in one forward pass.
 *
 * A window slides forward over the stored values, never backwards, so each
 * compressed chunk is decompressed at most once. Cells not in any group are
 * skipped outright — the window jumps to the next selected row — so a split
 * that touches a quarter of the cells reads far less than the whole file.
 */
function readGroupsCSR(
  m: H5, groups: Int32Array[], nCells: number, nGenes: number,
  onProgress?: (frac: number) => void,
): CellMajor[] {
  const indptr = toI32(child(m, 'indptr').value)
  const dataDs = child(m, 'data')
  const idxDs = child(m, 'indices')
  const total = indptr[nCells]

  // Which group each cell belongs to, and how big each group's arrays must be.
  const groupOf = new Int32Array(nCells).fill(-1)
  groups.forEach((cells, g) => { for (const c of cells) groupOf[c] = g })
  const out = groups.map(cells => {
    let nnz = 0
    for (const c of cells) nnz += indptr[c + 1] - indptr[c]
    return {
      nCells: cells.length, nGenes,
      indptr: new Int32Array(cells.length + 1),
      indices: new Int32Array(nnz),
      data: new Float32Array(nnz),
    } as CellMajor
  })
  const at = new Int32Array(groups.length)      // write cursor per group
  const row = new Int32Array(groups.length)     // next row index per group

  let winStart = 0, winEnd = 0
  let winData: Float32Array = new Float32Array(0)
  let winIdx: Int32Array = new Int32Array(0)

  for (let c = 0; c < nCells; c++) {
    const g = groupOf[c]
    if (g < 0) continue
    const p0 = indptr[c], p1 = indptr[c + 1]
    if (p0 < winStart || p1 > winEnd) {
      winStart = p0
      winEnd = Math.min(total, Math.max(p0 + SCAN_BLOCK, p1))
      winData = toF32(dataDs.slice([[winStart, winEnd]]))
      winIdx = toI32(idxDs.slice([[winStart, winEnd]]))
      onProgress?.(p0 / total)
    }
    const o = out[g]
    o.indptr[row[g]] = at[g]
    for (let k = p0; k < p1; k++) {
      const v = winData[k - winStart]
      if (v === 0) continue
      o.indices[at[g]] = winIdx[k - winStart]
      o.data[at[g]] = v
      at[g]++
    }
    row[g]++
  }
  // Explicit zeros were skipped, so the arrays may be shorter than reserved.
  return out.map((o, g) => {
    o.indptr[row[g]] = at[g]
    return at[g] === o.data.length ? o : {
      ...o, indices: o.indices.subarray(0, at[g]), data: o.data.subarray(0, at[g]),
    }
  })
}

/**
 * Read an embedding once and keep it.
 *
 * A 292 k-cell coordinate array is 2.3 MB — nothing beside the matrix — and a
 * split asks for it once per part per embedding. At 43 parts and two
 * embeddings that is 86 traversals of the same dataset through a block cache
 * that the matrix pass keeps evicting.
 */
function once(e: EmbeddingRef): EmbeddingRef {
  let held: Float32Array | null = null
  return { ...e, load: () => (held ??= e.load()) }
}

/** Enough values to classify a matrix, without reading all of it. */
function sampleMatrix(m: H5, want = 20000): Float64Array {
  if (isDataset(m)) {
    const shape = m.shape as number[]
    const rows = Math.min(shape[0], Math.max(1, Math.ceil(want / shape[1])))
    const step = Math.max(1, Math.floor(shape[0] / rows))
    const vals: number[] = []
    for (let r = 0; r < rows; r++) {
      const at = Math.min(shape[0] - 1, r * step)
      const row = m.slice([[at, at + 1], [0, shape[1]]]) as ArrayLike<number>
      for (let i = 0; i < row.length; i++) vals.push(Number(row[i]))
    }
    return Float64Array.from(vals)
  }
  const d = child(m, 'data')
  const n = (d.shape as number[])[0]
  // One strided read would mean n round trips into WASM; a few contiguous
  // blocks spread across the vector cost far less and sample just as fairly.
  const blocks = 8
  const per = Math.max(1, Math.floor(Math.min(want, n) / blocks))
  const vals: number[] = []
  for (let b = 0; b < blocks && vals.length < want; b++) {
    const at = Math.min(n - per, Math.floor((n / blocks) * b))
    if (at < 0) break
    const block = d.slice([[at, at + per]]) as ArrayLike<number>
    for (let i = 0; i < block.length; i++) vals.push(Number(block[i]))
  }
  return vals.length ? Float64Array.from(vals) : stridedSample(n, want, i => Number(d.slice([[i, i + 1]])[0]))
}

/** Whole rows of a CSR or dense matrix, run by run: what the per-cell totals are measured on. */
function readRuns(m: H5, nGenes: number, runs: [number, number][]): CellMajor {
  const rows: CellMajor[] = []
  for (const [a, b] of runs) {
    if (b <= a) continue
    if (isDataset(m)) {
      const flat = toF32(m.slice([[a, b], [0, nGenes]]))
      const indptr = new Int32Array(b - a + 1), idx: number[] = [], val: number[] = []
      for (let r = 0; r < b - a; r++) {
        indptr[r] = idx.length
        for (let g = 0; g < nGenes; g++) {
          const v = flat[r * nGenes + g]
          if (v !== 0) { idx.push(g); val.push(v) }
        }
      }
      indptr[b - a] = idx.length
      rows.push({ nCells: b - a, nGenes, indptr, indices: Int32Array.from(idx), data: Float32Array.from(val) })
      continue
    }
    const ptr = toI32(child(m, 'indptr').slice([[a, b + 1]]))
    const s0 = ptr[0], s1 = ptr[b - a]
    const data = s1 > s0 ? toF32(child(m, 'data').slice([[s0, s1]])) : new Float32Array(0)
    const indices = s1 > s0 ? toI32(child(m, 'indices').slice([[s0, s1]])) : new Int32Array(0)
    rows.push({ nCells: b - a, nGenes, indptr: ptr.map(v => v - s0), indices, data })
  }
  const nCells = rows.reduce((n, r) => n + r.nCells, 0)
  const nnz = rows.reduce((n, r) => n + r.data.length, 0)
  const out: CellMajor = { nCells, nGenes, indptr: new Int32Array(nCells + 1), indices: new Int32Array(nnz), data: new Float32Array(nnz) }
  let c = 0, at = 0
  for (const r of rows) {
    for (let i = 0; i < r.nCells; i++) out.indptr[c + i] = at + r.indptr[i]
    out.indices.set(r.indices, at)
    out.data.set(r.data, at)
    c += r.nCells
    at += r.data.length
  }
  out.indptr[nCells] = at
  return out
}

// ---------------------------------------------------------------------------
// What the object records of what was run (facts.ts).
//
// Read by key at named places, and nothing found by a pattern over names: the
// integration method and its embedding as `default_method` and
// `default_embedding` in a dict directly under uns (the reader names no tool);
// a doublet or an ambient-RNA tool as `doublet_method` or `ambient_method`
// there; the normalization as scanpy's own record, `uns/log1p`. Each is carried
// only as one line of MAX_VALUE characters or fewer. Where nothing is at the
// place the fact is "not recorded", with where the converter looked; two places
// that both record one are both named in the notes and neither is stated.

const WHAT: Record<RecordedId, string> = {
  integration: 'an integration method', integration_embedding: "an integration method's embedding",
  normalization: 'a normalization', doublets: 'a doublet tool', ambient: 'an ambient-RNA tool',
}

/**
 * What one place of the object holds, read once: THE reader of every recorded
 * fact (readRecorded below). It tells a value from nothing there, and never
 * makes nothing a value: no entry, a null (anndata writes None as a dataset
 * with no dataspace, which a read gives as null and Number() as 0), an empty
 * array, an empty or blank string, a number that is not finite. A read that
 * fails is its own answer: what is there is not known, so no fact follows.
 */
export type Held =
  | { is: 'value'; value: string | number | boolean }
  | { is: 'nothing'; what: string }
  | { is: 'unreadable'; what: string }
  | { is: 'other'; what: string }

function held(parent: H5, key: string): Held {
  let o: H5 | null
  try { o = child(parent, key) } catch (e) { return { is: 'unreadable', what: `an entry that could not be read (${(e as Error).message})` } }
  if (!o) return { is: 'nothing', what: 'no entry' }
  if (isGroup(o)) return { is: 'other', what: 'a group' }
  if (!isDataset(o)) return { is: 'other', what: 'neither a dataset nor a group' }
  let info: H5
  try { info = o.metadata ?? {} } catch { return { is: 'unreadable', what: 'a value whose type cannot be read' } }
  if (info.shape === null) return { is: 'nothing', what: 'a null value' }
  const shape = (info.shape as number[] | undefined) ?? []
  const size = shape.reduce((a, b) => a * b, 1)
  if (shape.length && size === 0) return { is: 'nothing', what: 'an empty array' }
  if (size > 1) return { is: 'other', what: `an array of ${size.toLocaleString('en-US')} values` }
  let v: unknown
  try { v = ownValue(o) } catch (e) { return { is: 'unreadable', what: `a value that could not be read (${(e as Error).message})` } }
  if (Array.isArray(v) || ArrayBuffer.isView(v)) v = (v as unknown as ArrayLike<unknown>)[0]
  if (typeof v === 'bigint') v = Number(v)
  if (v === null || v === undefined) return { is: 'nothing', what: 'a null value' }
  if (typeof v === 'number') return Number.isFinite(v) ? { is: 'value', value: v } : { is: 'nothing', what: `a number that is not finite (${v})` }
  if (typeof v === 'boolean') return { is: 'value', value: v }
  if (typeof v === 'string') return v.trim() ? { is: 'value', value: v } : { is: 'nothing', what: v ? 'a blank string' : 'an empty string' }
  return { is: 'other', what: 'a value that is not a string, a number or a boolean' }
}

/** What a Held is, for a note: "a group", "a number", "nothing (a null value)". */
const heldAs = (h: Held): string =>
  h.is === 'value' ? `a ${typeof h.value}` : h.is === 'nothing' ? `nothing (${h.what})` : h.what

function readRecorded(uns: H5, obsmKeys: string[], note: (s: string) => void): NonNullable<Sources['recorded']> {
  const out: NonNullable<Sources['recorded']> = {}
  const ids = Object.keys(RECORD_KEYS) as (keyof typeof RECORD_KEYS)[]
  const hits = new Map<RecordedId, { value: string; at: string }[]>(ids.map(id => [id, []]))
  // Places that hold something this converter does not carry, or that could
  // not be read: then the fact is not stated (each is named), never "not
  // recorded" of an object that may hold one.
  const refused = new Set<RecordedId>()
  const keys: string[] = isGroup(uns) ? uns.keys() : []
  for (const k of keys) {
    let o: H5 | null
    // An entry whose type cannot be read (a float16 attribute) is not a record
    // here; the run reader names it.
    try { o = child(uns, k); if (o && !isDict(o)) o = null } catch { o = null }
    if (!o) continue
    for (const id of ids) {
      const field = RECORD_KEYS[id]
      const at = `uns/${shown(k)}/${field}`
      const h = held(o, field)
      if (h.is === 'nothing') {
        // Nothing at the place records nothing: it is not a value, and not a refusal.
        if (h.what !== 'no entry') note(`${at} holds nothing (${h.what}), so it records no ${WHAT[id]}`)
      } else if (h.is === 'unreadable') {
        refused.add(id)
        note(`${at} could not be read (${h.what}), so the bundle states nothing of the object's ${WHAT[id]}`)
      } else if (h.is === 'other' || typeof h.value !== 'string') {
        refused.add(id)
        note(`${at} is ${heldAs(h)}, so it is not read as the object's record of ${WHAT[id]}`)
      } else if (!isLine(h.value, MAX_VALUE) || !isLine(at, MAX_PLACE)) {
        refused.add(id)
        note(`${at} holds a value that is not one line of ${MAX_VALUE} characters or fewer without markup, so it is not carried as the object's record of ${WHAT[id]}; it stays in the object`)
      } else {
        hits.get(id)!.push({ value: h.value, at })
      }
    }
  }
  for (const id of ids) {
    const h = hits.get(id)!
    if (h.length === 1) out[id] = recorded(h[0].value, h[0].at)
    else if (!h.length && !refused.has(id)) out[id] = notRecorded(SEARCHED[id])
    else if (!h.length) continue
    else note(`${h.length} places record ${WHAT[id]}: ${h.map(x => `${x.at} (${x.value})`).join(', ')}; the bundle states none of them as the object's`)
  }
  const e = out.integration_embedding
  if (e && 'at' in e && !obsmKeys.includes(e.value)) note(`${e.at} names ${e.value}, which is not an embedding of this object`)

  // scanpy's record of log1p, a dict at uns/log1p: its base is stated only
  // from a real number at uns/log1p/base (scanpy's default is None).
  const top: Held = isGroup(uns) ? held(uns, 'log1p') : { is: 'nothing', what: 'no entry' }
  let n: H5 | null = null
  if (top.is === 'other' && top.what === 'a group') {
    try { n = child(uns, 'log1p'); if (n && !isDict(n)) n = null } catch { n = null }
  }
  if (n) {
    const b = held(n, 'base')
    const base = b.is === 'value' && typeof b.value === 'number' ? b.value : null
    if (base === null && !(b.is === 'nothing' && b.what === 'no entry')) {
      note(`uns/log1p/base is ${heldAs(b)}, not a number, so the normalization is stated with no base`)
    }
    out.normalization = recorded(base === null ? 'log1p' : `log1p, base ${base}`, 'uns/log1p')
  } else if (top.is === 'unreadable') {
    note(`uns/log1p could not be read (${top.what}), so the bundle states nothing of the object's normalization`)
  } else {
    if (!(top.is === 'nothing' && top.what === 'no entry')) {
      note(`uns/log1p is ${top.is === 'other' && top.what === 'a group' ? 'a group of another encoding than dict' : heldAs(top)}, not a dict as scanpy records log1p, so no normalization is read from it`)
    }
    out.normalization = notRecorded(SEARCHED.normalization)
  }
  return out
}

// ---------------------------------------------------------------------------
// The run's identity.
//
// An object says in uns what made it, and every tool says it in its own place:
// the integration stage under uns/<tool>/provenance, the annotator as
// uns/<column>_provenance, the QC stage as a group whose schema names
// provenance, and the harness's run record (when there is one) as
// uns/harness_run. The reader names no tool: a record is what calls itself
// provenance, by its key or by its own schema, at the top of uns or one level
// down.
//
// Of each record it carries the run's identity, not the record: six fields
// read by name, each only as a string of RUN_FIELD_MAX characters or fewer,
// under the converter's own keys. They are the fields Studio's Run card shows
// first. Everything else a record holds (a tool's arguments and input path,
// package versions, settings, label counts) stays in the object, and is counted,
// never named. An earlier design copied whatever a record held, so its input
// became all of HDF5 (links, cycles, other files, types JSON cannot hold, names
// that are code in JavaScript) for a card that shows six strings.

/** Where the reader looks for the run's records. Written into the bundle beside what it found. */
export const RUN_SEARCHED = 'uns/harness_run; every entry of uns, and every child of a dict directly under uns, named provenance or *_provenance; every dict directly under uns whose own schema names provenance'

/** The fields read from each record, by name and in this order: the run's identity. */
export const RUN_FIELDS: readonly RunField[] = ['tool', 'version', 'commit', 'stage', 'run_key', 'job']

/** The longest value carried, in characters (code points). */
export const RUN_FIELD_MAX = 256

/** The most records read. Past it, records are counted in leftOut and the first ten named. */
export const RUN_MAX_RECORDS = 64

/** `a`; `a and b`; `a, b and c`. */
function and(items: string[]): string {
  return items.length < 2 ? items.join('') : `${items.slice(0, -1).join(', ')} and ${items[items.length - 1]}`
}

/** Paths for a sentence: the first ten, then a count. */
function some(paths: string[]): string {
  const shown = paths.slice(0, 10).join(', ')
  return paths.length > 10 ? `${shown} and ${(paths.length - 10).toLocaleString('en-US')} more` : shown
}

/** Is `s` longer than `max` characters, counted as code points? */
function longer(s: string, max: number): boolean {
  if (s.length <= max) return false
  let n = 0
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i)
    if (c >= 0xd800 && c <= 0xdbff && i + 1 < s.length) {
      const d = s.charCodeAt(i + 1)
      if (d >= 0xdc00 && d <= 0xdfff) i++
    }
    if (++n > max) return true
  }
  return false
}

/**
 * What an HDF5 object is, from its type and shape alone, said of a field that
 * is not carried; null for a scalar string, which is read and then measured.
 */
function kindOf(o: H5): string | null {
  if (isGroup(o)) return 'a group'
  if (!isDataset(o)) return 'neither a dataset nor a group'
  let info: H5
  try {
    info = o.metadata ?? {}
  } catch {
    return 'a value whose type cannot be read'
  }
  // anndata writes None as a dataset with no dataspace: a null, not a number.
  if (info.shape === null) return 'a null value'
  const shape = (info.shape as number[] | undefined) ?? []
  if (shape.length) {
    const n = shape.reduce((a, b) => a * b, 1)
    return `an array of ${n.toLocaleString('en-US')} value${n === 1 ? '' : 's'}`
  }
  if (info.type === 3) return null
  if (info.type === 0 || info.type === 1) return 'a number'
  const members = info.enum_type?.members as Record<string, number> | undefined
  if (members && Object.keys(members).length === 2 && 'TRUE' in members && 'FALSE' in members) return 'a boolean'
  return 'a value that is not a string'
}

/**
 * A dict in anndata's layout, or a plain group from a writer that names no encoding.
 *
 * (Tailwind reads this directory for class names, so a variable here is never
 * called `t` and negated: the app has a `t` class, and the scan took the
 * negation for its important variant and shipped it as CSS.)
 */
const isDict = (g: H5): boolean => {
  const encoding = text(attr(g, 'encoding-type') ?? '')
  return isGroup(g) && (encoding === '' || encoding === 'dict')
}

/** Does this group name provenance in its own `schema` string? */
function schemaSaysProvenance(o: H5): boolean {
  if (!isDict(o)) return false
  try {
    const s = child(o, 'schema')
    return isDataset(s) && /provenance/i.test(text(ownValue(s)))
  } catch {
    return false
  }
}

/** One record, read by the closed rule. */
interface Identity {
  /** The fields carried, under the converter's own keys, in RUN_FIELDS order. */
  record: Partial<Record<RunField, string>>
  /** Each of the six it holds that is not carried, with what it is. */
  refused: string[]
  /** How many of its entries are not one of the six: counted, never named. */
  others: number
  /** Where a value longer than RUN_FIELD_MAX was met, which ends the section's reading. */
  long: string | null
}

/**
 * The six fields of one record, looked up by name. Nothing else of the record
 * is read: its entries are counted, never opened. A field that is a link out
 * of the file is absent here, and named by the scan's note on links.
 */
function readIdentity(rec: H5, path: string): Identity {
  const record: Partial<Record<RunField, string>> = {}
  const refused: string[] = []
  const six: readonly string[] = RUN_FIELDS
  const others = (rec.keys() as string[]).filter(k => !six.includes(k)).length
  for (const field of RUN_FIELDS) {
    const o = child(rec, field)
    if (!o) continue
    const kind = kindOf(o)
    if (kind) {
      refused.push(`${field} (${kind})`)
      continue
    }
    let v: unknown
    try {
      v = ownValue(o)
    } catch {
      refused.push(`${field} (a string that could not be read from this file)`)
      continue
    }
    const s = typeof v === 'string' ? v : Array.isArray(v) && v.length === 1 && typeof v[0] === 'string' ? v[0] : null
    if (s === null) {
      refused.push(`${field} (a value that is not a string)`)
    } else if (longer(s, RUN_FIELD_MAX)) {
      refused.push(`${field} (a string longer than ${RUN_FIELD_MAX} characters)`)
      // Nothing after it is read. The fields of the six this record holds
      // after it are named, by the converter's own names, and never opened.
      const holds = rec.keys() as string[]
      for (const later of RUN_FIELDS.slice(RUN_FIELDS.indexOf(field) + 1)) {
        if (holds.includes(later)) refused.push(`${later} (not read: it comes after that string)`)
      }
      return { record, refused, others, long: `${path}/${field}` }
    } else {
      record[field] = s
    }
  }
  return { record, refused, others, long: null }
}

/**
 * The identity of the run that made the object, record by record.
 *
 * Records are found in uns order by RUN_SEARCHED's rule, and the first
 * RUN_MAX_RECORDS of them are read by `readIdentity`. One that carries none of
 * the six is not written, and its note says so; one that is not a dict (a
 * string such as the harness's JSON, an array, a group of another encoding) is
 * named with what it is, and nothing of it is read. The first value longer than
 * RUN_FIELD_MAX ends the reading, so a long value is read once, however many
 * records reach it. The records past the cap and after that value are counted
 * in `leftOut`, in uns order, the first ten named. When no record is listed, so
 * is every other record found (none of the six, not a dict, unreadable): a
 * reader of the bundle that sees no record still learns that one was found.
 * Beside a listed record those are named in the notes only. Nothing here can
 * stop the conversion.
 */
function readRun(uns: H5, note: (s: string) => void): RunSection {
  const records: RunRecord[] = []
  if (!isGroup(uns)) return { searched: RUN_SEARCHED, records }
  const named = (k: string) => k === 'provenance' || k.endsWith('_provenance')
  const n = (x: number) => x.toLocaleString('en-US')

  // ---- found, in uns order; each entry inspected on its own -------------
  const found: { path: string; o: H5 }[] = []
  const uninspected: string[] = []
  let why = ''
  for (const k of uns.keys() as string[]) {
    const path = `uns/${shown(k)}`
    try {
      const o = child(uns, k)
      if (!o) continue
      if (k === 'harness_run' || named(k) || schemaSaysProvenance(o)) {
        found.push({ path, o })
      } else if (isDict(o)) {
        for (const c of o.keys() as string[]) {
          if (!named(c)) continue
          const oc = child(o, c)
          if (oc) found.push({ path: `${path}/${shown(c)}`, o: oc })
        }
      }
    } catch (e) {
      // An attribute h5wasm cannot read (float16): this entry cannot be
      // classified, and the conversion goes on without it.
      uninspected.push(path)
      why ||= (e as Error).message
    }
  }

  // ---- read, in uns order, to the cap -------------------------------------
  const pastCap: string[] = []
  const afterLong: string[] = []
  const failed: string[] = []
  const unlisted: string[] = []
  let long: string | null = null
  for (let i = 0; i < found.length; i++) {
    const { path, o } = found[i]
    if (long !== null) { afterLong.push(path); unlisted.push(path); continue }
    if (i >= RUN_MAX_RECORDS) { pastCap.push(path); unlisted.push(path); continue }
    try {
      if (!isDict(o)) {
        const what = isGroup(o) ? 'a group of another encoding than dict' : kindOf(o) ?? 'a string'
        note(`${path} is ${what}, not a dict of named fields, so run reads nothing from it`)
        unlisted.push(path)
        continue
      }
      const id = readIdentity(o, path)
      const carried = Object.keys(id.record)
      const not = [...id.refused, ...(id.others ? [`${n(id.others)} other entr${id.others === 1 ? 'y' : 'ies'}`] : [])]
      const rest = not.length ? `; not carried: ${and(not)}` : ''
      if (carried.length) {
        records.push({ path, record: id.record })
        note(`${path}: run carries its ${and(carried)}${rest}`)
      } else {
        note(`${path} holds none of ${and([...RUN_FIELDS])} as a string of ${RUN_FIELD_MAX} characters or fewer, so run does not list it${rest}`)
        unlisted.push(path)
      }
      long = id.long
    } catch (e) {
      failed.push(`${path} (${(e as Error).message})`)
      unlisted.push(path)
    }
  }

  // ---- said, once each -------------------------------------------------------
  if (uninspected.length) {
    note(`${some(uninspected)} could not be inspected (${why}), so the run section cannot tell whether ${uninspected.length > 1 ? 'they hold records' : 'it holds a record'}`)
  }
  if (failed.length) {
    note(`${some(failed)}: could not be read, so run carries nothing of ${failed.length > 1 ? 'them' : 'it'}`)
  }
  if (long !== null) {
    note(`${long} is a string longer than ${RUN_FIELD_MAX} characters, so run read no value after it${afterLong.length ? `: ${n(afterLong.length)} record${afterLong.length > 1 ? 's' : ''} found after it ${afterLong.length > 1 ? 'were' : 'was'} not read: ${some(afterLong)}` : ''}`)
  }
  if (pastCap.length) {
    note(`${n(found.length)} records were found and run reads the first ${RUN_MAX_RECORDS}: the ${n(pastCap.length)} after them ${pastCap.length > 1 ? 'were' : 'was'} not read: ${some(pastCap)}`)
  }
  // Beside a listed record, `leftOut` holds only what the budget left out (past
  // the cap, after the long value): a reader of the bundle says those were left
  // out over its budget, which is false of the others, and the notes name them.
  // When none is listed, it holds every record found, so that a reader that
  // sees no record still learns that one was found.
  const leftOut = records.length ? [...afterLong, ...pastCap] : unlisted
  return {
    searched: RUN_SEARCHED, records,
    ...(leftOut.length ? { leftOut: { count: leftOut.length, paths: leftOut.slice(0, 10) } } : {}),
  }
}

export interface H5Module {
  ready: Promise<H5Lib>
  File: new (name: string, mode: string) => H5
}

/** Read an .h5ad far enough to ask the user which columns to use. */
/**
 * Scan an already-open HDF5 handle.
 *
 * Split from the entry points below because how the file got opened — copied
 * into memory, or mounted lazily from a File inside a worker — changes nothing
 * about how it is read.
 */
export function scanOpenH5ad(f: H5, filename: string): Scan {
  if (!UNFOLLOWED.has(String(f?.file_id))) {
    throw new H5adError('the converter reads only a file it opened itself (openInput), so that no link can lead it to another')
  }
  const notes: string[] = []
  const note = (s: string) => notes.push(s)

  const obsRaw = child(f, 'obs')
  if (!obsRaw) throw new H5adError(unfollowed(f, 'obs') ?? 'this HDF5 file has no /obs — it is not an .h5ad')

  const X = child(f, 'X')
  if (!X) throw new H5adError(unfollowed(f, 'X') ?? 'this file has no /X matrix')
  const shape = isDataset(X)
    ? (X.shape as number[])
    : (toI32(attr(X, 'shape') ?? attr(X, 'h5sparse_shape') ?? []) as unknown as number[])
  const nCells = Number(shape[0])
  const nGenes = Number(shape[1])
  if (!nCells || !nGenes) throw new H5adError('the /X matrix has no readable shape')

  const legacy = needed('the obs table', () => compound(obsRaw)) != null
  note(legacy
    ? 'read in the legacy AnnData layout — obs and obsm are compound datasets and categorical levels come from uns'
    : 'read in the modern AnnData layout')

  const uns = child(f, 'uns')
  const lost: Lost[] = []
  const columns = readObs(obsRaw, uns, nCells, note, lost)

  // ---- gene lists -------------------------------------------------------
  const geneSets: Record<string, string[]> = {}
  const geneAliases: Record<string, GeneAlias> = {}
  const varDs = child(f, 'var')
  const varGenes = needed('the gene names (var)', () => readVar(varDs, 'var'))
  const varSet = internGeneSet(geneSets, varGenes, 'var')
  const varAlias = readVarAlias(varDs, varGenes, note)
  if (varAlias) geneAliases[varSet] = varAlias
  const rawVarDs = child(f, 'raw/var') ?? child(f, 'raw.var')
  let rawSet: string | null = null
  if (rawVarDs) {
    try {
      const rawGenes = readVar(rawVarDs, 'raw.var')
      rawSet = internGeneSet(geneSets, rawGenes, 'raw')
      // internGeneSet reuses an identical list, so raw and var can land on the
      // same key — do not overwrite an alias that is already there.
      if (!geneAliases[rawSet]) {
        const a = readVarAlias(rawVarDs, rawGenes, note)
        if (a) geneAliases[rawSet] = a
      }
    } catch (e) {
      // No raw genes. When HDF5 could not read them, that is said: raw.X is not used.
      if (e instanceof UnreadableError) note(`raw.var could not be read (${e.message}); raw.X is not used`)
    }
  }

  // ---- candidate matrices ----------------------------------------------
  const matrices: MatrixRef[] = []
  const unreadable: string[] = []
  const lostMatrices: { key: string; why: string }[] = []
  const add = (key: string, m: H5, geneSet: string, cells: number, genes: number) => {
    if (!m) return
    try {
      // Every member a load will read, looked up now and its stored type
      // checked: what is not in this file, or not numbers, is named by the
      // scan, never met later by a load.
      if (isGroup(m)) {
        for (const part of ['data', 'indices', 'indptr']) {
          const d = child(m, part)
          if (!isDataset(d)) throw new Error(`its ${part} is not a dataset in this file`)
          const bad = notNumbers(d, part !== 'data')
          if (bad) throw new Error(`its ${part} is ${bad}`)
        }
      } else {
        const bad = notNumbers(m)
        if (bad) throw new Error(`it is ${bad}`)
      }
      const csr = !isDataset(m) &&
        !text(attr(m, 'encoding-type') ?? attr(m, 'h5sparse_format') ?? '').includes('csc')
      // A read of a matrix after the scan, when the bundle is built from it,
      // that HDF5 reports as failed stops the conversion, by name: it never
      // falls to another matrix.
      const chosen = <A extends unknown[], T>(read: (...a: A) => T) => (...a: A): T =>
        needed(shown(key), () => read(...a))
      matrices.push({
        key, geneSet, nCells: cells, nGenes: genes,
        kind: classify(sampleMatrix(m)),
        load: chosen(() => readMatrix(m, cells, genes, key)),
        ...(csr || isDataset(m) ? { readRuns: chosen((runs: [number, number][]) => readRuns(m, genes, runs)) } : {}),
        ...(csr ? {
          nnzPerCell: chosen(() => rowLengths(m, cells)),
          // Only CSR: in a CSC matrix `indices` are cell ids, and counting
          // distinct cells would answer a different question entirely.
          coverage: () => geneCoverage(m, genes),
          loadGroups: chosen((groups: Int32Array[], onProgress?: (frac: number) => void) =>
            readGroupsCSR(m, groups, cells, genes, onProgress)),
        } : {}),
      })
    } catch (e) {
      const why = (e as Error).message
      note(`${shown(key)} could not be read (${why}); skipped`)
      unreadable.push(`${shown(key)} (${why})`)
      lostMatrices.push({ key, why })
    }
  }

  add('X', X, varSet, nCells, nGenes)
  const layers = child(f, 'layers')
  if (isGroup(layers)) {
    for (const k of layers.keys()) add(`layers["${k}"]`, child(layers, k), varSet, nCells, nGenes)
  }
  const rawX = child(f, 'raw/X') ?? child(f, 'raw.X')
  if (rawX && rawSet) {
    add('raw.X', rawX, rawSet, nCells, geneSets[rawSet].length)
  }
  if (!matrices.length) throw new H5adError(`no matrix in this object can be read, so the conversion stops: ${unreadable.join('; ')}`)

  // ---- embeddings -------------------------------------------------------
  const embeddings: EmbeddingRef[] = []
  const obsm = child(f, 'obsm')
  let comp: Compound | null = null
  try {
    comp = obsm ? compound(obsm) : null
  } catch (e) {
    note(`obsm could not be read (${(e as Error).message}); no embedding is offered from it`)
  }
  if (comp) {
    comp.members.forEach((m: H5, i: number) => {
      const width = m.array_type?.shape?.[0] ?? 0
      if (width < 2) return
      embeddings.push(once({
        key: m.name, nDims: width,
        load: () => {
          const xy = new Float32Array(nCells * 2)
          for (let c = 0; c < nCells; c++) {
            const row = comp.rows[c][i] as ArrayLike<number>
            xy[2 * c] = Number(row[0])
            xy[2 * c + 1] = Number(row[1])
          }
          return xy
        },
      }))
    })
  } else if (isGroup(obsm)) {
    for (const k of obsm.keys()) {
      const ds = child(obsm, k)
      if (!isDataset(ds)) continue
      // Its stored type checked, then read now: an embedding that is not
      // numbers, or that HDF5 cannot read, is named here and never offered.
      const bad = notNumbers(ds)
      if (bad) {
        note(`obsm/${shown(k)} is ${bad}; it is not offered as an embedding`)
        continue
      }
      const s = ds.shape as number[]
      if (s.length < 2 || s[1] < 2 || s[0] !== nCells) continue
      let xy: Float32Array
      try {
        const flat = toF32(ds.slice([[0, nCells], [0, 2]]))
        xy = flat.length === nCells * 2 ? flat : toF32(ds.value).filter((_v, i) => i % s[1] < 2)
      } catch (e) {
        note(`obsm/${shown(k)} could not be read (${(e as Error).message}); it is not offered as an embedding`)
        continue
      }
      embeddings.push({ key: k, nDims: s[1], load: () => xy })
    }
  }

  // ---- what the object records of what was run -------------------------
  // Read by key at named places, never from the look of a name (facts.ts).
  const recorded = readRecorded(uns, embeddings.map(e => e.key), note)
  const provenance = legacyProvenance(recorded)
  const run = readRun(uns, note)

  // ---- what was not followed, named once -----------------------------------
  const all = [...(UNFOLLOWED.get(String(f.file_id)) ?? [])]
  const outside = all.filter(([, says]) => says !== FILTERED)
  const filtered = all.filter(([, says]) => says === FILTERED).map(([path]) => shown(path))
  if (filtered.length) {
    note(`${some(filtered)} ${filtered.length > 1 ? 'are datasets' : 'is a dataset'} compressed with a filter this reader does not have, LZF for one, so ${filtered.length > 1 ? 'they are' : 'it is'} not read: write the file again with gzip compression`)
  }
  if (outside.length) {
    note(`the converter reads only the file it is given and follows no link out of it: ${some(outside.map(([path, says]) => `${shown(path)} (${says})`))} ${outside.length > 1 ? 'are' : 'is'} not read`)
  }

  return {
    format: 'h5ad', source: filename, nCells, columns, matrices, embeddings,
    geneSets, geneAliases, provenance, recorded, run, notes,
    ...(lost.length ? { unreadableColumns: lost } : {}),
    ...(lostMatrices.length ? { unreadableMatrices: lostMatrices } : {}),
  }
}

// ---------------------------------------------------------------------------
// Two ways in.

/**
 * Copy the whole file into h5wasm's in-memory filesystem, then scan it.
 *
 * Simple and synchronous once the bytes are in hand, but it holds the file
 * twice and cannot exceed MAX_H5AD_BYTES. Used by the Node tests; the app now
 * takes the lazy path below.
 */
export async function scanH5ad(
  h5wasm: H5Module, bytes: Uint8Array, filename: string,
): Promise<Scan> {
  const lib = await h5wasm.ready
  const FS = lib.FS

  // h5wasm is a 32-bit WebAssembly build of HDF5, and its emulated filesystem
  // holds the whole file. Lengths there are signed 32-bit, so 2 GB is a hard
  // ceiling — measured exactly: 2047 MB writes, 2048 MB does not.
  if (bytes.length >= MAX_H5AD_BYTES) {
    throw new H5adError(
      `this .h5ad is ${(bytes.length / 1e9).toFixed(1)} GB, which is past the 2 GB limit of the in-memory reader. Open it in the app, which mounts the file instead of copying it, or run tools/export_h5ad.py offline.`)
  }

  const path = `scan_${filename.replace(/[^\w.-]/g, '_')}`
  try { FS.unlink(path) } catch { /* first run */ }
  try {
    throughDoor(FS, () => FS.writeFile(path, bytes))
  } catch {
    throw new H5adError(
      `this .h5ad (${(bytes.length / 1e9).toFixed(1)} GB) would not fit in the reader's memory.`)
  }
  return scanOpenH5ad(openOrFail(h5wasm, lib, path), filename)
}

function openOrFail(h5wasm: H5Module, lib: H5Lib, path: string): H5 {
  try {
    return openInput(h5wasm, lib, path)
  } catch (e) {
    if (e instanceof H5adError) throw e
    throw new H5adError('this file is not readable as HDF5. An .h5ad is an HDF5 file — if this came from Seurat, save it as .rds instead.')
  }
}

/**
 * Mount a File lazily and scan it, at any size.
 *
 * WORKERFS is Emscripten's read-only filesystem over Blob/File objects: it
 * serves reads straight from the File with FileReaderSync rather than copying
 * it into memory, so nothing is ever held except the bytes HDF5 actually asks
 * for. That is what lifts the 2 GB ceiling — a 9.3 GB atlas opens in
 * milliseconds because opening it reads only the superblock and the metadata.
 *
 * FileReaderSync exists only inside a worker, so this must run in one.
 */
/** Bytes fetched per block. Big enough that the per-read cost disappears. */
const BLOCK = 4 << 20
/** Blocks kept resident: 48 x 4 MB, so the cache is bounded at ~192 MB. */
const MAX_BLOCKS = 48

/**
 * Present a File to HDF5 as a file it can read, without copying it anywhere.
 *
 * Emscripten ships WORKERFS for exactly this, and it is unusable here: it
 * answers every read with a fresh Blob slice through FileReaderSync, and that
 * call costs about 4.4 ms however few bytes it returns. HDF5 walks its metadata
 * in thousands of small reads, so a scan that takes 2 s against a local file
 * had still not finished after ten minutes through WORKERFS. Measured on the
 * 9.3 GB atlas: 2000 reads of 4 KB took 8.8 s - 0.9 MB/s - while one large read
 * ran at 20 MB/s. The cost is per call, not per byte.
 *
 * So reads are served out of 4 MB blocks with an LRU behind them. HDF5's small
 * reads then land in memory, and the 4.4 ms is paid once per block instead of
 * once per read.
 *
 * FileReaderSync exists only inside a worker, so this must run in one.
 */
function mountBlockCached(FS: H5, file: File, dir: string, name: string): string {
  const fr = new FileReaderSync()
  const cache = new Map<number, Uint8Array>()

  const blockAt = (i: number): Uint8Array => {
    const hit = cache.get(i)
    if (hit) {
      cache.delete(i)        // reinsert, so iteration order is least-recent first
      cache.set(i, hit)
      return hit
    }
    const start = i * BLOCK
    const slice = file.slice(start, Math.min(file.size, start + BLOCK))
    const block = new Uint8Array(fr.readAsArrayBuffer(slice) as ArrayBuffer)
    cache.set(i, block)
    if (cache.size > MAX_BLOCKS) cache.delete(cache.keys().next().value as number)
    return block
  }

  try { FS.mkdir(dir) } catch { /* left from a previous file */ }
  try { FS.unlink(`${dir}/${name}`) } catch { /* first time */ }

  const node = FS.createFile(dir, name, {}, true, false)
  // stat() reports `usedBytes`, not `size`. Setting only the latter leaves HDF5
  // looking at what it believes is an empty file — which it does not report as
  // an error, it just never finds a superblock and reads nothing, forever.
  Object.defineProperty(node, 'usedBytes', { get: () => file.size })
  node.size = file.size
  node.stream_ops = {
    llseek(stream: H5, offset: number, whence: number): number {
      const base = whence === 1 ? stream.position : whence === 2 ? node.size : 0
      const pos = base + offset
      if (pos < 0) throw new FS.ErrnoError(28)
      return pos
    },
    read(_s: H5, buffer: Uint8Array, offset: number, length: number, position: number): number {
      if (position >= node.size) return 0
      const end = Math.min(node.size, position + length)
      let pos = position, done = 0
      while (pos < end) {
        const i = Math.floor(pos / BLOCK)
        const block = blockAt(i)
        const from = pos - i * BLOCK
        const n = Math.min(end - pos, block.length - from)
        if (n <= 0) break
        buffer.set(block.subarray(from, from + n), offset + done)
        done += n
        pos += n
      }
      return done
    },
  }
  return `${dir}/${name}`
}

/**
 * Open a File lazily and scan it, at any size.
 *
 * Nothing is held but the blocks HDF5 actually asks for, so the 2 GB ceiling of
 * the copy-it-all path does not apply here.
 */
export async function scanH5adFile(h5wasm: H5Module, file: File): Promise<{
  scan: Scan; handle: H5
}> {
  const lib = await h5wasm.ready
  if (typeof FileReaderSync === 'undefined') {
    throw new H5adError('a file this large can only be read inside a worker, and this is not one.')
  }
  const safe = file.name.replace(/[^\w.-]/g, '_')
  const handle = openOrFail(h5wasm, lib, mountBlockCached(lib.FS, file, '/lazy', safe))
  return { scan: scanOpenH5ad(handle, file.name), handle }
}
