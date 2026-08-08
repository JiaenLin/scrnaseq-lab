# scRNA-seq Lab

**One job: turn an already-annotated Scanpy `.h5ad` or Seurat `.rds` into the file
[scRNA-seq Studio](https://jiaenlin.github.io/scrnaseq-studio/) opens** — in the browser, with
no Python, no R, and no upload.

It does not cluster, annotate or analyse anything. Bring an object you have already processed;
you explore it in the Studio, not here.

**→ [jiaenlin.github.io/scrnaseq-lab](https://jiaenlin.github.io/scrnaseq-lab/)**

Drop the object in. The lab scans its metadata, shows every column it found, and asks the two
questions it cannot answer for you: **which column holds the cell type annotation**, and
**which — if any — is the condition to group by**. Then it writes `bundle.zip`.

## Why this exists

The studio reads a bundle, not the object itself, because resolving the two formats' quirks
in the viewer would mean shipping that ambiguity to every figure. Conversion used to mean
running `export_h5ad.py` or `export_seurat.R` from a terminal — fine if you have the
environment, a wall if you do not. This is the same conversion, in a tab.

Everything stays local. The file is read by JavaScript in your browser; there is no server to
send it to.

## What it reads

| Format | Layouts | Read by |
|---|---|---|
| `.h5ad` | modern AnnData (obs/var as groups) and legacy (< 0.7, compound datasets, `uns/<col>_categories`) | [h5wasm](https://github.com/usnistgov/h5wasm), a WASM build of HDF5 |
| `.rds` | Seurat v3/v4 (`@counts`, `@data`) and v5 (`@layers`), gzip or uncompressed | `src/lib/rds.ts` — a direct reader for R's serialization format |

The `.rds` path parses R's format itself rather than shipping webR. webR is a ~30 MB download
and a whole R runtime to read one file; the reader here is one source file, and it does not
need Seurat installed — it reads the object structurally, so it works on files written by a
Seurat version you can no longer install.

It is fast because it never materializes what it will not use. A Seurat object's `scale.data`
slot alone can be a gigabyte of doubles nobody asked for, so numeric vectors are recorded as
(offset, length) and skipped while walking. Measured on the 288 MB `pbmc3k_final.rds`:

| | |
|---|---|
| gunzip 288 MB → 351 MB | ~0.9 s |
| walk the whole object graph | ~60 ms |
| scan → questions on screen (in-browser, end to end) | ~2.6 s |
| convert → `bundle.zip` (8.9 MB, 2.2 M stored values) | ~0.6 s |

## What the object needs

| Needs | Which is | If missing |
|---|---|---|
| A cell annotation | any categorical column | **Required** — you pick which one after the scan |
| An embedding | UMAP, t-SNE or PCA | **Required** by the Cells and Feature plot tabs |
| Expression | log-normalized values, or raw counts to build them from | **Required** |
| A sample column | donor, animal, run | Optional — without it, composition cannot show between-animal spread |
| A condition column | the experimental group | Optional — without it the object opens single-condition |

## Decisions it makes, and shows

Slot names lie. An `.h5ad` routinely keeps scaled values in `X` and `log1p(counts)` in `.raw`;
picking by name picks wrong. So matrices are classified **by their numbers**:

| | |
|---|---|
| `scaled` | has negatives — refused, because the studio plots expression, not z-scores |
| `counts` | integer |
| `log-counts` | `expm1` is integer — logged, never depth-normalized; exponentiated back so pseudobulk works |
| `lognorm` | small positive non-integers — used as-is |
| `linear` | large positive non-integers — `log1p` applied |

Every such decision is written into `meta.notes` and appears again on the studio's Overview
tab, so whoever reads the figures sees the same caveats as whoever ran the conversion.

The lab also refuses to let one column play two roles: a column chosen as the cell annotation
is withdrawn from the condition and sample questions — as is any column that turns out to be
the *same partition under different labels*, which `seurat_clusters` and `seurat_annotations`
almost always are.

## Output

`bundle.zip`, schema `scrnaseq-studio/bundle@1` — the same format written by the studio's
`tools/export_h5ad.py` and `tools/export_seurat.R`, documented in
[tools/BUNDLE.md](https://github.com/JiaenLin/scrnaseq-studio/blob/main/tools/BUNDLE.md).
This is the third independent writer of it, so `scripts/test-build.mjs` reads a built bundle
back and compares it against a dense reference rather than trusting the writer.

## The family — one job each

| App | Takes | Produces |
|---|---|---|
| **rnaseq-service** | raw FASTQ | an analysis request + nf-core sample sheet |
| **rnaseq-lab** | a bulk counts matrix | `bundle.zip` for rnaseq-studio |
| **rnaseq-studio** | `bundle.zip` | the figures you read |
| **scrnaseq-lab** (here) | an annotated `.h5ad` / `.rds` | `bundle.zip` for scrnaseq-studio |
| **scrnaseq-studio** | `bundle.zip` | the figures you read |

## Development

```bash
npm install
npm run dev
npm test        # RDS reader, bundle builder, detection — no fixtures required
npm run build
```

The tests write R's serialization format themselves and read it back, so the reader is checked
against an independent implementation of the same spec rather than against a recorded file.

## Known limits

- **`.h5ad` files must be under 2 GB.** h5wasm is a 32-bit WebAssembly build of HDF5 and its
  emulated filesystem holds the whole file with signed-32-bit lengths, so 2^31-1 bytes is a hard
  ceiling (measured: 2047 MB writes, 2048 MB throws). The lab refuses larger files up front with
  that explanation rather than failing deep inside the reader. Subset the object first, or run
  `tools/export_h5ad.py` offline, which streams and has no such limit.
- **Very large objects are out of range regardless of format.** The conversion holds the matrix
  in memory twice — once cell-major, once gene-major — so it needs roughly `nnz x 16` bytes.
  A 700 M-nonzero atlas would want ~11 GB and produce a multi-GB bundle the studio could not
  open either. Convert a subset.
- **Dense `.h5ad` X** is materialized in full before being thinned to nonzeros. Sparse is the
  normal case and is streamed; a large dense matrix will be memory-hungry.
- **Seurat v5** objects are read through `@layers`. Tested against v3/v4 objects; a v5 object
  with a non-standard layer layout may need the assay named explicitly.
- **bzip2/xz `.rds`** is refused with instructions to re-save as gzip (R's default).
- Cell-level metadata beyond cluster/sample/condition is not carried into the bundle.
