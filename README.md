# scRNA-seq Lab

**One job: turn an already-annotated Scanpy `.h5ad` or Seurat `.rds` into the file
[scRNA-seq Studio](https://jiaenlin.github.io/scrnaseq-studio/) opens** — in the browser, with
no Python, no R, and no upload.

It does not cluster, annotate or analyse anything. Bring an object you have already processed;
you explore it in the Studio, not here.

**→ [jiaenlin.github.io/scrnaseq-lab](https://jiaenlin.github.io/scrnaseq-lab/)**

Drop the object in. The lab scans its metadata, shows every column it found, and asks the two
questions it cannot answer for you: **which column holds the cell type annotation**, and
**which — if any — is the condition to group by**. Then it writes one `.zip`.

Nothing else is asked. Any size is accepted, and an object too large to hold at once is split
automatically into that one file — a storage detail the studio hides completely.

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
| `.h5ad` | modern AnnData (obs/var as groups) and legacy (< 0.7, compound datasets, `uns/<col>_categories`) | [h5wasm](https://github.com/usnistgov/h5wasm) in a worker, over a block-cached reader — **no size limit** |
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

## Atlases

**There is no size limit, and nothing to decide.** Drop a 9.3 GB object in; you
get one `.zip` back.

The reader never copies the file. It is mounted behind a 4 MB block cache and
HDF5 pulls only the bytes it asks for, so opening a 9.3 GB object costs about
what opening a small one costs.

An object too large to hold at once is **split automatically**, into one file.
The lab costs every categorical column from the CSR row offsets — exactly, before
reading a single value — and takes the best cut itself. You are never asked. The
studio opens the result as one intact dataset with every cell, so the split is a
byte layout and never something you see.

Measured end to end in Chrome on the 9.3 GB, 292,495-cell, 736 M-nonzero
`developing_mouse_nervous_system.h5ad`:

| | |
|---|---|
| scan, questions on screen | **1.5–13 s** |
| cost every possible split, exactly | **12 ms** — from the row offsets, no values read |
| convert to one file | **~13 min** (43 parts, 7 forward passes, 5.83 GB, ZIP64) |
| the studio opens it | **4.0 s**, all 292,495 cells |

Two measurements shaped the design, both on that file:

- **Random row reads cost 20 ms each.** The payload is chunked and compressed, so
  a scattered row decompresses a whole chunk for ~2,500 values. One pass that way
  would take 100 minutes.
- **Sequential streaming sustains 11 M nonzeros/s.** So every read is
  forward-only, a window slides and never seeks backwards, and every part in a
  pass is filled together.

## What comes out

One `.zip`, always — whether it holds one piece or forty-three, so the studio has
exactly one shape to understand.

| | |
|---|---|
| `collection.json` | what it holds, what it was split along, and why |
| `parts/<key>.zip` | each a complete `scrnaseq-studio/bundle@1` |

Entries are **STORED, never deflated**: the parts are already compressed, and
storing them verbatim makes each one a contiguous byte range, so the studio reads
the index from the file's tail and pulls out only what it needs. Past 4 GB it
goes ZIP64 — a real atlas needs it, and without it the offsets wrap and the file
reads back as "not a collection".

Inside each part the expression matrix is also written **gene-chunked**: blocks of
64 genes, each deflated, with an offset index, so one gene can be read out of a
multi-gigabyte file without holding the matrix. Gap-encoding the cell ids made
that copy *smaller* than the flat one — −49.5% on real atlas data. One gene costs
**3.9 ms cold, 0.018 ms warm**.

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

- **`.rds` is still read whole.** The lazy path is `.h5ad` only: R's serialization is a stream
  with no index, so it cannot be read out of order. Very large Seurat objects remain limited by
  memory — save as `.h5ad` if yours is big.
- **Splitting needs a categorical column to split along.** An object with hundreds of millions
  of values and no usable grouping cannot be divided.
- **Dense `.h5ad` X** is materialized in full before being thinned to nonzeros. Sparse is the
  normal case and is streamed.
- **Writing a file larger than a few hundred MB needs `showSaveFilePicker`** (Chrome, Edge).
  A Chromium `Blob` past that size spills to disk and then fails every read, so elsewhere the
  fallback path is limited to what a Blob can carry.
