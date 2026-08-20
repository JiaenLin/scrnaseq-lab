// Reading a Seurat object out of the parsed .rds tree.
//
// A Seurat object is an S4 object, which in R's serialization is nothing but a
// bag of attributes — one per slot. So this file is mostly knowing the slot
// names and what shape each holds:
//
//   @assays$RNA@counts / @data     dgCMatrix, genes × cells, CSC by cell
//   @meta.data                     a data.frame — a list of columns
//   @reductions$umap@cell.embeddings   dense matrix, cells × dims, column-major
//   @active.ident                  a factor, often where the annotation ended up

import {
  attrOf, classOf, columnAsStrings, factorLabels, namesOf, RdsError, RdsReader,
  slot, strings, type RNum, type RValue,
} from './rds.ts'
import {
  categorical, fromFactor, internGeneSet, maybeDiscrete, numeric, pickGeneAlias,
  type Column, type EmbeddingRef, type GeneAlias, type MatrixRef, type Scan,
} from './scan.ts'
import { classify, dropZeros, stridedSample, type CellMajor } from './matrix.ts'

/** A dgCMatrix, read as cell-major sparse — which is how it is already stored. */
function dgc(r: RdsReader, m: RValue | undefined, key: string): CellMajor | null {
  if (!m) return null
  const cls = classOf(m)
  if (!cls.some(c => /^(dgCMatrix|dgTMatrix|dgRMatrix)$/.test(c))) return null
  if (cls.includes('dgTMatrix') || cls.includes('dgRMatrix')) {
    throw new RdsError(
      `${key} is a ${cls[0]}, not the dgCMatrix Seurat normally writes. Re-save after as(mat, "dgCMatrix").`)
  }
  const dim = slot(m, 'Dim')
  const i = slot(m, 'i'), p = slot(m, 'p'), x = slot(m, 'x')
  if (dim?.t !== 'num' || i?.t !== 'num' || p?.t !== 'num' || x?.t !== 'num') return null
  const d = r.numbers(dim) as Int32Array
  const nGenes = d[0], nCells = d[1]
  const indptr = r.numbers(p) as Int32Array
  if (indptr.length !== nCells + 1) {
    throw new RdsError(`${key} has ${indptr.length} column pointers for ${nCells} cells`)
  }
  return dropZeros({
    nCells, nGenes,
    indptr: Int32Array.from(indptr),
    indices: r.numbers(i) as Int32Array,
    data: r.floats(x),
  })
}


/**
 * Per-cell non-zero counts, straight off the column pointers.
 *
 * A dgCMatrix is CSC with one column per cell, so this is a difference of
 * `p` — a megabyte for a quarter-million cells — and a split can be costed
 * exactly before a single value is read.
 */
function dgcNnzPerCell(indptr: Int32Array, nCells: number): Int32Array {
  const out = new Int32Array(nCells)
  for (let c = 0; c < nCells; c++) out[c] = indptr[c + 1] - indptr[c]
  return out
}

/**
 * Several disjoint cell subsets, read in one forward pass.
 *
 * The reason the .rds path can convert an object it cannot hold. A dgCMatrix
 * stores cells as COLUMNS, so its values already arrive cell by cell: with `p`
 * in hand — which is small — every non-zero can be attributed to a cell as it
 * streams past, and written straight into the array of whichever shard that
 * cell belongs to. Nothing whole is ever materialised. The 618-million-nonzero
 * matrix that cannot be one array becomes twenty that can.
 *
 * `i` precedes `x` in the slot order R writes, and both are read from the same
 * fresh stream in that order, so one extra decompression serves both.
 */
function dgcGroups(
  r: RdsReader, value: RValue | undefined, key: string,
  groups: Int32Array[], onProgress?: (frac: number) => void,
): CellMajor[] {
  const dim = slot(value, 'Dim')
  const iv = slot(value, 'i') as RNum | undefined
  const pv = slot(value, 'p') as RNum | undefined
  const xv = slot(value, 'x') as RNum | undefined
  if (dim?.t !== 'num' || iv?.t !== 'num' || pv?.t !== 'num' || xv?.t !== 'num') {
    throw new RdsError(`${key} is not a readable dgCMatrix`)
  }
  const d = r.numbers(dim) as Int32Array
  const nGenes = d[0], nCells = d[1]
  const indptr = Int32Array.from(r.numbers(pv) as Int32Array)
  const nnz = indptr[nCells]

  const groupOf = new Int32Array(nCells).fill(-1)
  groups.forEach((cells, g) => { for (const c of cells) groupOf[c] = g })
  const out: CellMajor[] = groups.map(cells => {
    let k = 0
    for (const c of cells) k += indptr[c + 1] - indptr[c]
    return {
      nCells: cells.length,
      nGenes,
      indptr: new Int32Array(cells.length + 1),
      indices: new Int32Array(k),
      data: new Float32Array(k),
    }
  })
  // Each shard's cells keep the order the shard lists them in, which is the
  // order the split chose; the pass below visits cells in FILE order, so the
  // write position for a cell has to be looked up rather than counted.
  const slotOf = new Int32Array(nCells).fill(-1)
  groups.forEach((cells, g) => cells.forEach((c, j) => { slotOf[c] = j; void g }))
  const cursor = groups.map(cells => {
    const starts = new Int32Array(cells.length + 1)
    for (let j = 0; j < cells.length; j++) {
      starts[j + 1] = starts[j] + (indptr[cells[j] + 1] - indptr[cells[j]])
    }
    return starts
  })
  out.forEach((m, g) => m.indptr.set(cursor[g]))

  /**
   * Only a payload that was passed over needs the file opened again.
   *
   * A matrix small enough to hold is already in hand — `numbers()` returns it —
   * and re-decompressing sixteen gigabytes to re-read four megabytes would be
   * absurd. So the stream is created only if one of the two slots was deferred,
   * and both are then served from it.
   */
  const needsStream = iv.deferred || xv.deferred
  if (needsStream && !r.reopen) {
    throw new RdsError(`${key} was read past and there is no way to read it again`)
  }
  const st = needsStream ? r.reopen!() : null
  const BITE = 1 << 22
  const idx = new Int32Array(BITE)
  const val = new Float64Array(BITE)

  // Pass over `i`, then over `x`, in that order — which is the order R wrote
  // them, so one forward stream serves both.
  for (const which of ['i', 'x'] as const) {
    const v = which === 'i' ? iv : xv
    const wide = which === 'x'
    const held = v.deferred ? null : r.numbers(v)
    let k = 0
    let cell = 0
    while (k < nnz) {
      const take = Math.min(BITE, nnz - k)
      if (held) {
        for (let t = 0; t < take; t++) {
          if (wide) val[t] = held[k + t]
          else idx[t] = held[k + t]
        }
      } else if (wide) st!.read(val, v.off + k * 8, take, true)
      else st!.read(idx, v.off + k * 4, take, false)
      for (let t = 0; t < take; t++) {
        // The column this non-zero belongs to. k only ever increases, so the
        // pointer only ever walks forward — no search per value.
        while (cell < nCells && indptr[cell + 1] <= k + t) cell++
        const g = cell < nCells ? groupOf[cell] : -1
        if (g < 0) continue
        const at = cursor[g][slotOf[cell]] + (k + t - indptr[cell])
        if (wide) out[g].data[at] = val[t]
        else out[g].indices[at] = idx[t]
      }
      k += take
      if (onProgress) onProgress(((which === 'i' ? 0 : 1) + k / nnz) / 2)
    }
  }
  return out
}

/** Gene names off a dgCMatrix's Dimnames. */
function dgcGenes(m: RValue | undefined): string[] {
  const dn = slot(m, 'Dimnames')
  return dn?.t === 'list' ? strings(dn.v[0]) : []
}

const dgcDim = (r: RdsReader, m: RValue | undefined): [number, number] | null => {
  const dim = slot(m, 'Dim')
  if (dim?.t !== 'num') return null
  const d = r.numbers(dim) as Int32Array
  return [d[0], d[1]]
}

/** A dense R matrix as interleaved x,y — R stores column-major. */
function embedding2D(r: RdsReader, m: RValue | undefined, nCells: number): Float32Array {
  if (m?.t !== 'num') return new Float32Array(nCells * 2)
  const dim = attrOf(m, 'dim')
  const rows = dim?.t === 'num' ? (r.numbers(dim) as Int32Array)[0] : nCells
  const out = new Float32Array(nCells * 2)
  for (let c = 0; c < nCells; c++) {
    out[2 * c] = r.at(m, c)
    out[2 * c + 1] = r.at(m, rows + c)
  }
  return out
}

/**
 * The other naming of an assay's features, if the object carries one.
 *
 * Seurat has no dedicated place for it: `@meta.features` is a data.frame with
 * one row per feature, and what is in it depends entirely on what the analyst
 * put there. So every string column is offered to the same content-based
 * chooser the .h5ad path uses, and usually nothing qualifies — which is the
 * right answer for an object whose rows are already symbols.
 */
function assayAlias(
  r: RdsReader, assay: RValue | undefined, genes: string[],
): GeneAlias | null {
  const mf = slot(assay, 'meta.features')
  if (mf?.t !== 'list') return null
  const names = namesOf(mf)
  const columns: { name: string; values: string[] }[] = []
  mf.v.forEach((c, i) => {
    const lev = strings(attrOf(c, 'levels'))
    if (lev.length && c.t === 'num') {
      const codes = r.numbers(c) as Int32Array
      columns.push({
        name: names[i] ?? `V${i}`,
        values: Array.from(codes, k => (k >= 1 && k <= lev.length ? lev[k - 1] : '')),
      })
      return
    }
    const s = c.t === 'str' ? c.v.map(x => x ?? '') : columnAsStrings(r, c)
    if (s) columns.push({ name: names[i] ?? `V${i}`, values: s })
  })
  return pickGeneAlias(genes, columns)
}

/** Every assay's matrix slots, v3/v4 and v5 alike. */
function assayMatrices(
  assays: RValue | undefined, note: (s: string) => void,
): { key: string; value: RValue; genes: string[]; assay: RValue }[] {
  const out: { key: string; value: RValue; genes: string[]; assay: RValue }[] = []
  const names = namesOf(assays)
  if (assays?.t !== 'list') return out

  assays.v.forEach((a, ai) => {
    const an = names[ai] ?? `assay${ai}`
    const cls = classOf(a)

    // Seurat v5 keeps its matrices in a `layers` list instead of named slots.
    const layers = slot(a, 'layers')
    if (layers?.t === 'list') {
      const ln = namesOf(layers)
      layers.v.forEach((m, li) => {
        const genes = dgcGenes(m)
        out.push({ key: `${an}@layers$${ln[li] ?? li}`, value: m, genes, assay: a })
      })
      if (!out.length) note(`assay ${an} is ${cls[0] ?? 'v5'} with no readable layers`)
      return
    }

    for (const s of ['counts', 'data', 'scale.data']) {
      const m = slot(a, s)
      if (!m || m.t === 'null') continue
      out.push({ key: `${an}@${s}`, value: m, genes: dgcGenes(m), assay: a })
    }
  })
  return out
}

/** Read a Seurat .rds far enough to ask which columns to use. */
export function scanSeurat(r: RdsReader, obj: RValue, filename: string): Scan {
  const notes: string[] = []
  const note = (s: string) => notes.push(s)

  const cls = classOf(obj)
  if (!cls.includes('Seurat')) {
    throw new RdsError(
      cls.length
        ? `this .rds holds a ${cls.join('/')} object, not a Seurat object. If it is a SingleCellExperiment, convert it first: Seurat::as.Seurat(sce).`
        : 'this .rds does not hold an S4 object, so it is not a Seurat object.')
  }

  const version = slot(obj, 'version')
  const ver = version?.t === 'list' ? strings(version.v[0])[0] : strings(version)[0]
  if (ver) note(`written by Seurat ${ver}`)

  // ---- matrices ---------------------------------------------------------
  const found = assayMatrices(slot(obj, 'assays'), note)
  if (!found.length) throw new RdsError('this Seurat object has no assay matrices')

  const geneSets: Record<string, string[]> = {}
  const geneAliases: Record<string, GeneAlias> = {}
  const matrices: MatrixRef[] = []
  let nCells = 0

  for (const { key, value, genes, assay } of found) {
    const d = dgcDim(r, value)
    if (!d) {
      // scale.data is a plain dense matrix, not a dgCMatrix. It is scaled by
      // definition, so it is never the matrix we want — noting it is enough.
      note(`${key} is dense (scaled data); skipped`)
      continue
    }
    const [g, c] = d
    nCells = Math.max(nCells, c)
    if (genes.length !== g) {
      note(`${key} has no gene names on its Dimnames`)
      continue
    }
    const setKey = internGeneSet(geneSets, genes, key.replace(/@.*/, '') || 'features')
    if (!(setKey in geneAliases)) {
      const alias = assayAlias(r, assay, genes)
      if (alias) geneAliases[setKey] = alias
    }
    const x = slot(value, 'x') as RNum | undefined
    if (!x || x.t !== 'num') continue
    const pv = slot(value, 'p') as RNum | undefined
    matrices.push({
      key, geneSet: setKey, nGenes: g, nCells: c,
      kind: classify(stridedSample(x.n, 20000, i => r.at(x, i))),
      // The two the split path needs. With them an object far larger than the
      // tab can hold converts into a collection, exactly as an .h5ad atlas
      // does; without them the .rds path could only ever read a whole matrix
      // into one array, and 618 million non-zeros do not fit in one.
      nnzPerCell: pv && pv.t === 'num' && !pv.deferred
        ? () => dgcNnzPerCell(Int32Array.from(r.numbers(pv) as Int32Array), c)
        : undefined,
      loadGroups: r.reopen
        ? (groups, onProgress) => dgcGroups(r, value, key, groups, onProgress)
        : undefined,
      load: () => {
        const m = dgc(r, value, key)
        if (!m) throw new RdsError(`${key} could not be read`)
        return m
      },
    })
  }
  if (!matrices.length) {
    throw new RdsError('no assay in this object holds a sparse matrix with gene names')
  }

  // ---- meta.data --------------------------------------------------------
  const md = slot(obj, 'meta.data')
  const columns: Column[] = []
  const mdNames = namesOf(md)
  if (md?.t === 'list') {
    md.v.forEach((c, i) => {
      const name = mdNames[i] ?? `V${i}`
      const lev = strings(attrOf(c, 'levels'))
      if (lev.length && c.t === 'num') {
        // R factor codes are 1-based; ours are 0-based, with -1 for NA.
        const codes = r.numbers(c) as Int32Array
        const out = new Int32Array(codes.length)
        for (let k = 0; k < codes.length; k++) {
          out[k] = codes[k] >= 1 && codes[k] <= lev.length ? codes[k] - 1 : -1
        }
        columns.push(fromFactor(name, out, lev))
      } else if (c.t === 'str') {
        columns.push(categorical(name, c.v.map(s => s ?? 'NA')))
      } else if (c.t === 'num') {
        columns.push(maybeDiscrete(numeric(name, r.floats(c)), nCells))
      } else {
        const s = columnAsStrings(r, c)
        if (s) columns.push(categorical(name, s))
      }
    })
  }

  // active.ident is where RenameIdents leaves the annotation, and it is often
  // the only place it exists — a tutorial object relabelled but never written
  // back to meta.data would otherwise scan as having no cell types at all.
  const ident = factorLabels(r, slot(obj, 'active.ident'))
  if (ident && ident.length === nCells) {
    if (!columns.some(c => c.kind === 'categorical' &&
        c.levels.length === new Set(ident).size &&
        c.codes.every((v, i) => c.levels[v] === ident[i]))) {
      columns.push(categorical('active.ident', ident))
      note('active.ident is included as a column — RenameIdents does not write back to meta.data')
    }
  }

  // ---- reductions -------------------------------------------------------
  const red = slot(obj, 'reductions')
  const redNames = namesOf(red)
  const embeddings: EmbeddingRef[] = []
  if (red?.t === 'list') {
    red.v.forEach((d, i) => {
      const e = slot(d, 'cell.embeddings')
      if (e?.t !== 'num') return
      const dim = attrOf(e, 'dim')
      const shape = dim?.t === 'num' ? (r.numbers(dim) as Int32Array) : null
      if (!shape || shape[1] < 2) return
      embeddings.push({
        key: redNames[i] ?? `reduction${i}`, nDims: shape[1],
        load: () => embedding2D(r, e, nCells),
      })
    })
  }

  const names = columns.map(c => c.name)
  return {
    format: 'seurat', source: filename, nCells, columns, matrices, embeddings,
    geneSets, geneAliases,
    provenance: {
      normalization: null,
      integration: embeddings.some(e => /harmony|integrated|mnn/i.test(e.key)) ? 'present' : null,
      doublets: names.find(n => /doublet|scrublet/i.test(n)) ?? null,
      ambient: names.find(n => /soupx|cellbender|decontx/i.test(n)) ?? null,
    },
    notes,
  }
}
