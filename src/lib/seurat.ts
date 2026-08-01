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
  categorical, fromFactor, internGeneSet, maybeDiscrete, numeric,
  type Column, type EmbeddingRef, type MatrixRef, type Scan,
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

/** Every assay's matrix slots, v3/v4 and v5 alike. */
function assayMatrices(
  assays: RValue | undefined, note: (s: string) => void,
): { key: string; value: RValue; genes: string[] }[] {
  const out: { key: string; value: RValue; genes: string[] }[] = []
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
        out.push({ key: `${an}@layers$${ln[li] ?? li}`, value: m, genes })
      })
      if (!out.length) note(`assay ${an} is ${cls[0] ?? 'v5'} with no readable layers`)
      return
    }

    for (const s of ['counts', 'data', 'scale.data']) {
      const m = slot(a, s)
      if (!m || m.t === 'null') continue
      out.push({ key: `${an}@${s}`, value: m, genes: dgcGenes(m) })
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
  const matrices: MatrixRef[] = []
  let nCells = 0

  for (const { key, value, genes } of found) {
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
    const x = slot(value, 'x') as RNum | undefined
    if (!x || x.t !== 'num') continue
    matrices.push({
      key, geneSet: setKey, nGenes: g, nCells: c,
      kind: classify(stridedSample(x.n, 20000, i => r.at(x, i))),
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
    format: 'seurat', source: filename, nCells, columns, matrices, embeddings, geneSets,
    provenance: {
      normalization: null,
      integration: embeddings.some(e => /harmony|integrated|mnn/i.test(e.key)) ? 'present' : null,
      doublets: names.find(n => /doublet|scrublet/i.test(n)) ?? null,
      ambient: names.find(n => /soupx|cellbender|decontx/i.test(n)) ?? null,
    },
    notes,
  }
}
