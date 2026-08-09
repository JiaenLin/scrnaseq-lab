# Developer script — NOT part of `npm test`.
#
# A tiny .h5ad shaped like the real atlas, so the browser can be pointed at the
# two format additions without a nine-gigabyte file and a quarter of an hour:
# rows indexed by Ensembl accessions with the symbols in a categorical var
# column, and two 2D embeddings in obsm. Duplicate and empty symbols are in
# there on purpose — they are the cases the bundle has to state rather than
# quietly resolve.
#
#   python scripts/make-fixture.py <out.h5ad>

import sys
import h5py
import numpy as np

out = sys.argv[1] if len(sys.argv) > 1 else "scan_accession_fixture.h5ad"
rng = np.random.default_rng(7)
n_cells, n_genes = 300, 40

# Counts, sparse enough to look like real data.
dense = rng.poisson(0.6, size=(n_cells, n_genes)).astype(np.float32)
rows, cols = np.nonzero(dense)
order = np.argsort(rows, kind="stable")
rows, cols = rows[order], cols[order]
data = dense[rows, cols]
indptr = np.zeros(n_cells + 1, dtype=np.int32)
np.add.at(indptr, rows + 1, 1)
indptr = np.cumsum(indptr).astype(np.int32)

accessions = [f"ENSMUSG{i:011d}" for i in range(n_genes)]
# Two accessions carrying one symbol, and one with no symbol at all.
symbols = [f"Gene{i}" for i in range(n_genes)]
symbols[3] = symbols[2]
symbols[7] = ""
levels = sorted(set(symbols))
codes = np.array([levels.index(s) for s in symbols], dtype=np.int8)

cell_types = np.array([f"Type {i % 4}" for i in range(n_cells)])
type_levels = sorted(set(cell_types.tolist()))
type_codes = np.array([type_levels.index(t) for t in cell_types], dtype=np.int8)
samples = np.array([f"s{i % 3}" for i in range(n_cells)])
sample_levels = sorted(set(samples.tolist()))
sample_codes = np.array([sample_levels.index(s) for s in samples], dtype=np.int8)

vs = h5py.special_dtype(vlen=str)

with h5py.File(out, "w") as f:
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"

    X = f.create_group("X")
    X.attrs["encoding-type"] = "csr_matrix"
    X.attrs["encoding-version"] = "0.1.0"
    X.attrs["shape"] = np.array([n_cells, n_genes], dtype=np.int64)
    X.create_dataset("data", data=data)
    X.create_dataset("indices", data=cols.astype(np.int32))
    X.create_dataset("indptr", data=indptr)

    obs = f.create_group("obs")
    obs.attrs["encoding-type"] = "dataframe"
    obs.attrs["encoding-version"] = "0.2.0"
    obs.attrs["_index"] = "_index"
    obs.attrs["column-order"] = np.array(["cell_type", "sample"], dtype=object)
    obs.create_dataset("_index", data=np.array([f"c{i}" for i in range(n_cells)], dtype=object), dtype=vs)
    for name, cds, lev in (("cell_type", type_codes, type_levels), ("sample", sample_codes, sample_levels)):
        g = obs.create_group(name)
        g.attrs["encoding-type"] = "categorical"
        g.attrs["ordered"] = False
        g.create_dataset("categories", data=np.array(lev, dtype=object), dtype=vs)
        g.create_dataset("codes", data=cds)

    var = f.create_group("var")
    var.attrs["encoding-type"] = "dataframe"
    var.attrs["encoding-version"] = "0.2.0"
    var.attrs["_index"] = "Accession"
    var.attrs["column-order"] = np.array(["Gene", "Chromosome"], dtype=object)
    var.create_dataset("Accession", data=np.array(accessions, dtype=object), dtype=vs)
    # The symbols, stored the way the real atlas stores them: a categorical.
    g = var.create_group("Gene")
    g.attrs["encoding-type"] = "categorical"
    g.attrs["ordered"] = False
    g.create_dataset("categories", data=np.array(levels, dtype=object), dtype=vs)
    g.create_dataset("codes", data=codes)
    # A decoy: a plausible-looking string column that is not a naming.
    chrom = var.create_group("Chromosome")
    chrom.attrs["encoding-type"] = "categorical"
    chrom.attrs["ordered"] = False
    chrom.create_dataset("categories", data=np.array(["1", "2", "X"], dtype=object), dtype=vs)
    chrom.create_dataset("codes", data=(np.arange(n_genes) % 3).astype(np.int8))

    obsm = f.create_group("obsm")
    obsm.attrs["encoding-type"] = "dict"
    obsm.create_dataset("X_UMAP", data=rng.normal(size=(n_cells, 2)).astype(np.float32))
    obsm.create_dataset("X_tSNE", data=(rng.normal(size=(n_cells, 2)) * 10).astype(np.float64))

    f.create_group("layers").attrs["encoding-type"] = "dict"
    f.create_group("uns").attrs["encoding-type"] = "dict"

print(f"wrote {out}: {n_cells} cells x {n_genes} genes, X_UMAP + X_tSNE, symbols in var/Gene")
