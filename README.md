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
| `.h5ad` | modern AnnData (obs/var as groups, including the nullable strings anndata writes by default on pandas ≥ 3) and legacy (< 0.7, compound datasets, `uns/<col>_categories`) | [h5wasm](https://github.com/usnistgov/h5wasm) in a worker, over a block-cached reader — **no size limit** |
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
| Gene symbols | a `var` column beside accession-named rows | Optional — without it, gene sets written as symbols will not match, and the lab says so |

## What it carries that the object had

**Every two-column embedding, each once, not only the one you picked.** An object usually
holds a UMAP and a t-SNE of the same cells; which to look at is the reader's question, not
the converter's. The one you choose opens first and is still `embed.f32`, so nothing that
reads bundles has to change; every other two-column embedding travels beside it at 8 bytes a
cell. Two kinds do not, and `meta.notes` names each one left out:

- **A latent space.** A PCA, a corrected PCA (`X_harmony`), or anything else wider than two
  columns: its first two columns are not a picture of the cells. You can still choose one
  to open first. It is then drawn by its first two columns, and the notes flag it.
- **A byte copy** of an embedding already offered. An integration stage writes its chosen
  method's UMAP twice, as `X_umap` and `X_umap_<method>`. It is offered once, under the name
  that opens first, or else the first in the object's order. Copies are compared on the
  whole object, so every part of a split object offers the same ones.

**The run's identity, not its record.** The lab carries the run's identity, not the run's
record. From each record it finds it reads six fields by name (`tool`, `version`, `commit`,
`stage`, `run_key`, `job`) and keeps each one stored as a string of 256 characters or fewer,
unchanged and under its own name. Nothing else in a record is read: its other fields are
counted and stay in the object. An absent field is not shown; a longer one, or one that is not
a string, is said so.

A record is anything that calls itself provenance. The lab names no tool:

- an entry named `provenance` or `*_provenance`, at the top of `uns` or one level down (an
  integration stage's `uns/<tool>/provenance`, an annotator's `uns/<column>_provenance`);
- a group directly under `uns` whose own `schema` string names provenance;
- the harness's run record, `uns/harness_run`.

The records go into `meta.run.records` in `uns` order, each as a path and a `record` that holds
only those of the six it carries. `meta.notes` says, for each record, which fields were carried,
each of the six that was not and what it is (a number, a boolean, an array, a group, a string
longer than 256 characters), and how many other entries the record holds: a count, never their
names. So a tool's arguments, input paths and host names do not travel. A record that holds none
of the six as such a string is not listed, and its note says so; one that is not a dict at all
(the harness's record is one JSON string) is named with what it is, and nothing of it is read.
An integration stage writes its `job` as a dict of scheduler, id and host, so that `job` is
named as a group and not carried.

At most 64 records are read, in `uns` order. Past that they are counted in `run.leftOut`, and
the first ten are named. The first string longer than 256 characters ends the reading, so that
one long value is read once, however many records reach it; the records after it are counted in
`run.leftOut` too. When the object carries no record, `run.records` is empty and a note says
so; `run.searched` says where the lab looked. A Seurat object has no `uns`, so its bundle has no
`run`. Readers that do not know the field ignore it. An `uns` entry that cannot even be
inspected (an attribute the reader cannot read) is named in a note and never stops the
conversion. What a record says is recorded, never vouched for: a file can name any tool and
commit.

**One file: the one it is given.** HDF5 by itself opens other files: the file an external link
names, a virtual dataset's source, a dataset's external raw storage. From the command those are
this machine's own files, and in the browser the files the tab holds. While the input is open,
the converter lets HDF5 open no other file. What the scan meets that leads out of the file is
read as absent and named once, in one note, with its target as the link states it: an external
link, whatever file it names (the input's own name included), a soft link that reaches nothing
in the file, and a virtual dataset whose values lie in another file. This holds for the
observations, the variables, the matrix and its layers, the embeddings and `uns` alike, so the
same input gives the same bundle, whatever sits beside it. An ordinary soft or hard link inside
the file still reads.

**Both namings of the genes.** `genes.txt` is whatever the matrix rows are indexed by —
`ENSMUSG00000038751` as often as `Sox2` — and `meta.geneIdKind` says which, judged from the
names themselves. Where the object also had the other naming in `var`, it is written line
for line beside it as `gene_alias.txt`. No mapping table ships with the app: one would be
species-specific and would go stale, and the file being converted already knows the answer.
Genes with no symbol keep their accession and are counted; several accessions sharing one
symbol keep their own rows and are counted — merging them would sum two genes under one
name.

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

**Every cell keeps its own condition.** The bundle records a condition per sample, not per
cell. Some samples hold cells of two conditions: a hashed lane carrying two treatments, or a
donor sampled before and after. Each such sample is split into one sample per condition,
named `sample|condition`, as the studio's own `tools/export_h5ad.py` does.

- Only those samples are renamed. A new name that another sample already has gets `#2`.
- The app says so before converting, and `meta.notes` names the samples that were split.
- With no sample column, the cells of each condition are one sample.
- The decision is made on the whole object before any split into parts, so a sample has one
  name in every part.

The lab used to give every cell its sample's first condition. That put whole groups of cells
in the wrong condition, with nothing on screen to show it.

## Output

`bundle.zip`, schema `scrnaseq-studio/bundle@1` — the same format written by the studio's
`tools/export_h5ad.py` and `tools/export_seurat.R`, documented in
[tools/BUNDLE.md](https://github.com/JiaenLin/scrnaseq-studio/blob/main/tools/BUNDLE.md).
This is the third independent writer of it, so `scripts/test-build.mjs` reads a built bundle
back and compares it against a dense reference rather than trusting the writer.

## From a shell: `scrnaseq-lab convert`

The same conversion as a command, for a workstation or a pipeline. It runs the app's own
modules, so each part is the bundle the app's builder writes, with two additions in its
`meta.json`: `sentinels`, and the command's own notes (a design factor left out, the sentinels'
cell counts). `collection.json` holds the app's fields plus `sentinels` and a `sha256` for each
part. Readers that do not know a field ignore it. The app itself is unchanged: it converts in
the browser and uploads nothing.

```bash
git clone https://github.com/JiaenLin/scrnaseq-lab && cd scrnaseq-lab
npm ci && npm link            # puts scrnaseq-lab on PATH; Node 22.18 or later

scrnaseq-lab convert final.h5ad final.zip \
  --cluster cell_type --sample sample \
  --condition factor_a --extra factor_b --interaction factor_a,factor_b \
  --embedding X_umap --sentinel EXCLUDED --sentinel UNRESOLVED --label my-object
```

Inside the clone, `npx scrnaseq-lab …` or `node bin/scrnaseq-lab.mjs …` work without
linking. The command runs the TypeScript sources directly, which Node does not do from inside
`node_modules`, so install it from a clone and not as a dependency. It reads `.h5ad`; convert
a Seurat `.rds` in the app.

| Flag | Role | When absent, the guess |
|---|---|---|
| `--cluster COL` | the cell annotation | the app's own pick: a named annotation, else the finest grouping |
| `--sample COL` | the sample column | a column named like a sample |
| `--condition COL` | the design factor the studio groups by first | a design factor, one named like a condition first |
| `--extra COL` (repeatable) | another column, carried under its own name | every other design factor |
| `--interaction A,B` (repeatable) | the interaction as a column `A:B`, levels `a1:b1` | the two carried factors' interaction, when they cross |
| `--embedding KEY` | the embedding that opens first | a UMAP, then t-SNE, a force layout, PCA |
| `--sentinel LEVEL` (repeatable) | a `--cluster` level that is not a cell type | `EXCLUDED` and `UNRESOLVED`, where present |
| `--label TEXT` | the dataset's name | the file name |

**A flag given is used as given, and a flag absent is guessed out loud.** Each guess is printed
on stderr as the flag that would produce it (`guessed --condition factor_a  (…)`). Pasting the
guesses back gives the same bytes. `none` opts out of `--sample`, `--condition`, `--extra`,
`--interaction` and `--sentinel`. A flag that names something the object lacks (a column, an
embedding, or a `--sentinel` level that no cell carries) is refused with exit 4 and the list of
what the object has. Nothing is skipped quietly. A `--sentinel` typo that only printed a warning
would leave the real sentinel counted as a cell type.

**Design factors are found by what they are, not what they are called.** A design factor is a
grouping that takes one value in every sample. Cell types vary inside a sample; a genotype, a
diet or an age does not. By default every factor travels: one as the condition and the rest as
extras, plus the interaction of a crossed pair. The studio offers each of them as the Group.
A factor that a flag leaves out is named in `meta.notes`, never dropped silently.

**Sentinels** stay in `clusters`, because the cells are real and are drawn. They are also
listed in `sentinels`, in `collection.json` and in each part, with a note counting their cells,
so a reader can leave them out of composition and tests. Without `--sentinel`, the guess marks
`EXCLUDED` and `UNRESOLVED` where cells carry them. With it, every level named must be carried
by some cell.

**The output is never the input.** The command refuses (exit 2) an output that is the input
file, however it is reached: the same path, a symlink to it, a path through a symlinked
directory, or a hard link. It also refuses an output whose `.partial` is the input.

**The same input gives the same bytes**, in any time zone and at any time. Zip entries carry a
fixed date, and nothing in the file records when it was made. stdout carries one JSON summary
with the file's `sha256`, and `collection.json` lists a `sha256` for each part. So
`shasum -a 256 final.zip` checks the whole file, and each part can be checked without
trusting the others.

stderr gets the guesses, the notes and the progress. The output is written beside its target
and renamed into place, so a failed conversion leaves no file behind. For a large object,
raise Node's heap: `NODE_OPTIONS=--max-old-space-size=8192 scrnaseq-lab convert …`.

| Exit | Meaning |
|---|---|
| 0 | converted |
| 1 | an internal error, or a Node that cannot run the command |
| 2 | usage: an unknown option, a missing value, a stray path, an output that is the input |
| 3 | input: missing, not HDF5, or not an `.h5ad` this reader understands |
| 4 | roles: a flag names a column or embedding the object lacks or cannot use, or a `--sentinel` level no cell carries |
| 5 | the object cannot become a bundle (scaled data only, no embedding, no annotation) |
| 6 | output: the file could not be written |

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
npm test        # RDS reader, bundle builder, detection, and the command run for real
npm run build
```

The tests write R's serialization format themselves and read it back, so the reader is checked
against an independent implementation of the same spec rather than against a recorded file.

The command's tests (`scripts/test-cli.mjs`) convert small synthetic objects in
`scripts/fixtures/`. The main one is shaped like the single-cell family's final object: two
crossed factors, the integration stage's `uns` record, and the annotator's sentinels. A copy of
it is written in anndata's default string encoding. `scripts/make-family-fixture.py` writes them from
made-up names, and Python is needed only to regenerate them:
`uv run --no-project --with anndata==0.13.2 --with pandas==3.0.5 python scripts/make-family-fixture.py`.
`scripts/test-integrity.mjs` converts copies of them that it edits itself through h5wasm. It
checks that every cell keeps its condition, that the run's identity travels and nothing else of
a record does, that no link leads it out of the file, and that only two-column embeddings are
offered. It also checks that nothing else in the fixture's bundle
moved. `uns-records.h5ad` holds run records h5wasm cannot write (a structured array with an
int64 field, float16, an HDF5 reference): they convert, and only their identity is read. `uns-unreadable.h5ad` holds an `uns` entry
whose `encoding-type` attribute is float16, which anndata itself refuses to read.

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
