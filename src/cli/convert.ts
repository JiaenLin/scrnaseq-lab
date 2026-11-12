// The whole conversion, headless: the app's own modules in the app's own order.
//
// scan → roles → chooseSplit → buildBundle → collectionPieces, exactly as the
// worker runs them, driven from a shell instead of a page. Nothing here
// reimplements the format: every part comes out of the app's buildBundle and
// the file out of its collectionPieces. What it adds is what a shell needs and
// a page does not: every role settable by a flag, every guess said aloud, an
// exit status per kind of failure, an output that either appears whole or not
// at all, and three fields readers that do not know them ignore — `sentinels`
// in each part and in collection.json, and a sha256 for each part there.
//
// Parts are parked on disk as they are written rather than held. The app must
// hold them — a browser has nowhere else to put them — but on an atlas that is
// several gigabytes of live arrays, and a shell on a shared machine should not
// page everything else out to keep them.

import h5wasm from 'h5wasm/node'
import { createHash } from 'node:crypto'
import {
  closeSync, existsSync, mkdtempSync, openSync, readFileSync, readSync, realpathSync, renameSync, rmSync,
  statSync, writeFileSync, writeSync,
} from 'node:fs'
import { tmpdir } from 'node:os'
import { basename, dirname, join, resolve } from 'node:path'
import {
  BuildError, buildBundle, oneConditionPerSample, streamableMatrix, type BundleMeta, type BuildResult,
} from '../lib/build.ts'
import { COLLECTION_SCHEMA, collectionPieces, type CollectionMeta, type PartInfo } from '../lib/collection.ts'
import { H5adError, openInput, scanOpenH5ad } from '../lib/h5ad.ts'
import type { CellMajor } from '../lib/matrix.ts'
import type { Scan } from '../lib/scan.ts'
import { measureObject, plan, planInfo } from '../lib/plan.ts'
import { shardScan } from '../lib/shard.ts'
import type { ConvertArgs } from './args.ts'
import { CliError, EXIT } from './errors.ts'
import { resolveRoles } from './roles.ts'

/** collection@1, plus what this writer adds. Readers that do not know a field ignore it. */
type Index = Omit<CollectionMeta, 'parts'> & {
  /** The whole object's level order for the condition and each extra column. */
  condOrder?: string[]
  extraOrder: Record<string, string[]>
  /** Cluster levels that are sentinels, not cell types. Empty means none. */
  sentinels: string[]
  /** Each part's own sha256, so a part can be checked without trusting the rest. */
  parts: (PartInfo & { sha256: string })[]
}

export interface Summary {
  output: string
  /** sha256 of the file written, hex. */
  sha256: string
  bytes: number
  parts: number
  label: string
  choices: {
    cluster: string
    sample: string | null
    condition: string | null
    extras: string[]
    interactions: string[][]
    embedding: string
    sentinels: string[]
  }
  /** Flags that were guessed rather than given. */
  guessed: string[]
  /** The bundle's notes: every decision the conversion made, in order. */
  notes: string[]
}

const HDF5_SIGNATURE = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]

/** The HDF5 superblock signature: at 0, or after a user block of 512, 1024, … bytes. */
function looksLikeHdf5(path: string): boolean {
  const fd = openSync(path, 'r')
  try {
    const buf = Buffer.alloc(8)
    for (let at = 0; at <= 1 << 20; at = at ? at * 2 : 512) {
      if (readSync(fd, buf, 0, 8, at) < 8) return false
      if (HDF5_SIGNATURE.every((b, i) => buf[i] === b)) return true
    }
    return false
  } finally {
    closeSync(fd)
  }
}

const sha256 = (b: Uint8Array) => createHash('sha256').update(b).digest('hex')
const safe = (s: string) => (s || 'part').replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 80)

function writeAll(fd: number, b: Uint8Array): void {
  let at = 0
  while (at < b.length) at += writeSync(fd, b, at, b.length - at)
}

/**
 * Is `p` the input file, however it is reached?
 *
 * The same path spelled another way, a symlink to it, a path through a
 * symlinked directory and a hard link all share the input's device and inode,
 * so that is what is compared; where a file system reports no inode, the real
 * paths are. A path that does not exist is no file at all.
 */
function isInput(p: string, input: string): boolean {
  try {
    const a = statSync(p, { bigint: true })
    const b = statSync(input, { bigint: true })
    if (a.ino !== 0n && b.ino !== 0n) return a.dev === b.dev && a.ino === b.ino
    return realpathSync.native(p) === realpathSync.native(input)
  } catch {
    return false
  }
}

/** Input and output checked before any work: an hour's conversion should not end in a typo. */
function checkPaths(a: ConvertArgs): void {
  if (!existsSync(a.input)) throw new CliError(`${a.input}: no such file`, EXIT.input)
  if (statSync(a.input).isDirectory()) throw new CliError(`${a.input} is a directory, not an .h5ad`, EXIT.input)
  // The zip is written to OUT.partial and renamed onto OUT. If either is the
  // input, the object is replaced by its own conversion, and that is the one
  // mistake the command cannot undo.
  for (const target of [a.output, `${a.output}.partial`]) {
    if (isInput(target, a.input)) {
      throw new CliError(`${target} is the input file ${a.input}${target === a.output ? '' : ' (the command writes OUT.partial first, then renames it)'}: converting onto it would replace the object with its own zip. Name another output`, EXIT.usage)
    }
  }
  if (/\.rds$/i.test(a.input)) {
    throw new CliError(`${a.input}: the command reads .h5ad. Convert a Seurat .rds in the app, which reads it in the browser, or save the object as .h5ad`, EXIT.input)
  }
  if (!looksLikeHdf5(a.input)) {
    throw new CliError(`${a.input} is not an HDF5 file, so it is not an .h5ad`, EXIT.input)
  }
  const dir = dirname(resolve(a.output))
  if (!existsSync(dir) || !statSync(dir).isDirectory()) {
    throw new CliError(`${a.output}: the directory ${dir} does not exist`, EXIT.output)
  }
  if (existsSync(a.output) && statSync(a.output).isDirectory()) {
    throw new CliError(`${a.output} is a directory; name the .zip to write`, EXIT.output)
  }
}

/** Errors from the shared modules, sorted into the command's exit codes. */
function classify(e: unknown): unknown {
  if (e instanceof CliError) return e
  if (e instanceof H5adError) return new CliError(e.message, EXIT.input)
  if (e instanceof BuildError) return new CliError(e.message, EXIT.build)
  return e
}

export async function convert(a: ConvertArgs, say: (line: string) => void): Promise<Summary> {
  checkPaths(a)
  const t0 = Date.now()
  const step = (s: string) => say(`  [${((Date.now() - t0) / 1000).toFixed(1).padStart(5)} s] ${s}`)

  // The input is the one file HDF5 may open while it is read: an external
  // link, a virtual dataset's source or external raw storage cannot reach the
  // host's other files.
  const f = openInput(h5wasm, await h5wasm.ready, a.input)
  if (!(Number(f.file_id) >= 0)) throw new CliError(`${a.input} could not be opened as HDF5`, EXIT.input)

  const scratch = mkdtempSync(join(tmpdir(), 'scrnaseq-lab-'))
  try {
    let scan: Scan
    try {
      scan = scanOpenH5ad(f, basename(a.input))
    } catch (e) {
      throw classify(e)
    }
    say(`${basename(a.input)}: ${scan.nCells.toLocaleString('en-US')} cells · ${scan.columns.length} obs columns · ${scan.matrices.length} matrices · embeddings ${scan.embeddings.map(e => `${e.key} (${e.nDims}-D)`).join(', ') || 'none'}`)

    const r = resolveRoles(scan, a)
    for (const line of r.said) say(line)
    for (const line of r.info) say(line)
    for (const line of r.notes) say(`note: ${line}`)
    const choices = r.choices

    // The split is the app's decision on the object's own columns, made before
    // the interactions are added, so that it never cuts along a column the
    // object did not have.
    const measurable = scan.matrices.find(m => m.nnzPerCell)
    // Its row lengths are read here, whole: one HDF5 cannot read is the input's
    // error, named like any other, not the command's.
    const lengths = (() => {
      try { return measurable?.nnzPerCell ? measurable.nnzPerCell() : null } catch (e) { throw classify(e) }
    })()
    // One plan, the one the page shows (lib/plan.ts): split or not, the matrix
    // read, what travels. Made on the object's own columns, before the
    // interactions are added.
    const p = plan({ nCells: scan.nCells, columns: [...scan.columns], ...planInfo(scan, lengths) }, choices, a.splitAt)
    const split = p.split
    scan.columns.push(...r.added)
    scan.notes.push(...r.notes)

    // A sample whose cells hold more than one condition becomes one sample per
    // condition, decided here on the whole object so that every part of a
    // split names it the same way. Said on the terminal as well as in the notes.
    const said = scan.notes.length
    scan = oneConditionPerSample(scan, choices)
    for (const line of scan.notes.slice(said)) say(`note: ${line}`)

    // Measured on the whole object before any part is cut, so every part states the same.
    const ms = (() => {
      try { return measureObject(scan, p, choices) } catch (e) { throw classify(e) }
    })()
    for (const line of ms.notes) say(`note: ${line}`)
    const ctx = { plan: p, measured: ms }

    const parked: { file: string; tmp: string }[] = []
    const partInfo: Index['parts'] = []
    const used = new Set<string>()
    const keep = (key: string, res: BuildResult): BundleMeta => {
      let file = `parts/${safe(key)}.zip`
      for (let i = 2; used.has(file); i++) file = `parts/${safe(key)}-${i}.zip`
      used.add(file)
      const tmp = join(scratch, `${parked.length}.zip`)
      writeFileSync(tmp, res.zip)
      parked.push({ file, tmp })
      partInfo.push({
        key, file, nCells: res.meta.nCells, nnz: res.meta.nnz, bytes: res.zip.length, sha256: sha256(res.zip),
      })
      step(`${file}: ${res.meta.nCells.toLocaleString('en-US')} cells, ${(res.zip.length / 1e6).toFixed(1)} MB`)
      return res.meta
    }

    let first: BundleMeta | undefined
    try {
      if (!split) {
        first = keep(choices.label, buildBundle(scan, choices, s => step(s), ctx))
      } else {
        step(`split along ${split.column.name}: ${split.plan.shards.length} parts in ${split.plan.passes.length} passes`)
        const src = streamableMatrix(scan, s => { scan.notes.push(s); step(s) })
        if (!src?.loadGroups) throw new CliError('no matrix in this object can be read in parts', EXIT.build)
        if (src.key !== p.stream) throw new CliError('the plan and the build chose different matrices to read in parts: this is a bug in the converter', EXIT.build)
        for (let pi = 0; pi < split.plan.passes.length; pi++) {
          const group = split.plan.passes[pi]
          step(`pass ${pi + 1} of ${split.plan.passes.length}: ${group.length} parts`)
          let last = -1
          const mats: (CellMajor | null)[] = src.loadGroups(group.map(s => Int32Array.from(s.cells)), frac => {
            const tenth = Math.floor(frac * 10)
            if (tenth !== last) { last = tenth; step(`reading ${(frac * 100).toFixed(0)}%`) }
          })
          for (let gi = 0; gi < group.length; gi++) {
            const s = group[gi]
            const sub = shardScan(scan, { ...s, cells: Int32Array.from(s.cells) },
              new Map([[src.key, mats[gi]!]]), `${scan.source} · ${choices.label}`)
            const meta = keep(s.key, buildBundle(sub, { ...choices, label: `${choices.label} ${s.key}` }, undefined, ctx))
            first ??= meta
            // A pass can hold a gigabyte across its parts; drop each as it is used.
            mats[gi] = null
          }
        }
      }
    } catch (e) {
      throw classify(e)
    }
    if (!first) throw new CliError('the conversion produced nothing', EXIT.build)

    const levelsOf = (name: string | null) => (name ? scan.columns.find(c => c.name === name)?.levels : undefined)
    const index: Index = {
      schema: COLLECTION_SCHEMA,
      label: choices.label,
      source: first.source,
      splitBy: split ? split.column.name : null,
      reason: split ? split.reason : null,
      nCells: partInfo.reduce((n, p) => n + p.nCells, 0),
      nGenes: first.nGenes,
      clusterOrder: levelsOf(choices.cluster),
      // The whole object's order for every column a part carries only part of,
      // so the studio does not fall back to collation (e16.25 before e16.5).
      condOrder: levelsOf(choices.condition),
      extraOrder: Object.fromEntries((choices.extras ?? []).flatMap(name => {
        const levels = levelsOf(name)
        return levels ? [[name, levels]] : []
      })),
      embeddings: first.embeddings?.map(e => e.key),
      geneIdKind: first.geneIdKind,
      geneAlias: first.geneAlias ? { kind: first.geneAlias.kind, column: first.geneAlias.column } : null,
      sentinels: r.sentinels,
      parts: partInfo,
      notes: first.notes,
    }

    // Written beside the target and renamed into place: a conversion that
    // fails half-way leaves nothing that looks like a finished file.
    const partial = `${a.output}.partial`
    const hash = createHash('sha256')
    let bytes = 0
    let fd: number | null = null
    try {
      fd = openSync(partial, 'w')
      let held = -1
      let heldBytes: Uint8Array = new Uint8Array(0)
      // collectionPieces touches a part's bytes several times (CRC, sizes,
      // payload) and finishes one before the next, so one part in memory at a
      // time is enough.
      const parts = parked.map((p, i) => ({
        file: p.file,
        get bytes() {
          if (held !== i) { heldBytes = readFileSync(p.tmp); held = i }
          return heldBytes
        },
      }))
      for (const piece of collectionPieces(index as CollectionMeta, parts)) {
        hash.update(piece)
        writeAll(fd, piece)
        bytes += piece.length
      }
      closeSync(fd)
      fd = null
      renameSync(partial, a.output)
    } catch (e) {
      if (fd !== null) closeSync(fd)
      rmSync(partial, { force: true })
      throw new CliError(`writing ${a.output} failed: ${(e as Error).message}`, EXIT.output)
    }
    const digest = hash.digest('hex')
    step(`wrote ${a.output}: ${(bytes / 1e6).toFixed(2)} MB in ${partInfo.length} part${partInfo.length > 1 ? 's' : ''}, sha256 ${digest}`)

    const interactionNames = new Set(r.added.map(c => c.name))
    return {
      output: a.output,
      sha256: digest,
      bytes,
      parts: partInfo.length,
      label: choices.label,
      choices: {
        cluster: choices.cluster,
        sample: choices.sample,
        condition: choices.condition,
        extras: (choices.extras ?? []).filter(x => !interactionNames.has(x)),
        interactions: r.interactions,
        embedding: choices.embedding,
        sentinels: r.sentinels,
      },
      guessed: r.guessed,
      notes: first.notes,
    }
  } finally {
    f.close()
    rmSync(scratch, { recursive: true, force: true })
  }
}
