"""Synthetic objects for the facts walk (scripts/test-facts.mjs), written by anndata.

    python scripts/make-facts-fixtures.py [--cells 240] [--out scripts/fixtures]

Each object copies a SHAPE of the family's final object and nothing else of it:
no value, name or label of any cohort. The shapes: X log-normalised to a
constant total per cell (CSR float32), the integer counts in layers["counts"],
the same values as X in layers["lognorm"]; two corrected latent spaces with
their two-column views, the method the object RECORDS being the second, and
X_umap a byte copy of that second view; a dict directly under uns holding
default_method and default_embedding beside a provenance record whose job is a
group among twenty other entries; four cell types and the annotator's two
sentinel labels, the cells of one of them with no coordinates on any embedding.

    facts-a.h5ad        the shape (the RIGHT object of most rows)
    facts-b-decoy.h5ad  the shape with every record removed and only NAMES left:
                        embeddings named like two methods, obs columns named like
                        a doublet tool and a mitochondrial share, uns keys named
                        like an ambient tool and a normalisation, a layer named
                        "counts" holding log values, a two-column X_pca
    facts-e-records.h5ad the shape plus records at the places the lab reads:
                        uns/log1p (scanpy's), doublet_method and ambient_method
                        in a dict under uns, per-cell QC columns in obs
    facts-f-two.h5ad    two dicts under uns each recording a method
    facts-g-lognorm.h5ad log values with no constant per-cell total, no counts,
                        and an uns string saying "CP10K"
    facts-h-hostile.h5ad default_method in dicts as a group, an array, a number,
                        10,000 characters, a line break, markup
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

p = argparse.ArgumentParser()
p.add_argument("--cells", type=int, default=240)
p.add_argument("--genes", type=int, default=40)
p.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "fixtures"))
a = p.parse_args()
n, g = a.cells, a.genes
rng = np.random.default_rng(20261112)

TYPES = ["type_a", "type_b", "type_c", "type_d"]
samples = np.array([f"s{i % 8 + 1}" for i in range(n)])
group = np.array(["g1" if int(s[1:]) <= 4 else "g2" for s in samples])
label = np.array([TYPES[i % 4] for i in range(n)], dtype=object)
label[np.arange(n) % 12 == 5] = "EXCLUDED"
label[np.arange(n) % 12 == 7] = "UNRESOLVED"
excluded = label == "EXCLUDED"

# genes: symbols, with a few accessions among them (the real object names some by accession)
genes = [f"GENE{i}" for i in range(g)]
for i in range(0, g, 10):
    genes[i] = f"ENSG{100000 + i:011d}"

base = rng.gamma(2.0, 2.0, size=g)
lam = np.outer(rng.uniform(0.5, 2.0, size=n), base)
for t, typ in enumerate(TYPES):
    lam[label == typ, t * 3:(t * 3) + 3] *= 6
counts = rng.poisson(lam).astype(np.int64)
counts[counts.sum(1) == 0, 1] = 1
C = sp.csr_matrix(counts)
tot = counts.sum(1, keepdims=True)
L = sp.csr_matrix(np.log1p(counts / tot * 1e4).astype(np.float32))


def emb(width):
    x = rng.normal(size=(n, width)).astype(np.float32)
    x[excluded] = np.nan
    return x


def provenance():
    rec = {"tool": "integrator", "version": "0.7.0", "commit": "0" * 40,
           "job": {"scheduler": "pbs", "id": "100.server", "host": "node1"}}
    for i in range(20):
        rec[f"entry_{i:02d}"] = f"value {i}"
    return rec


def obs_frame():
    return pd.DataFrame({
        "sample": pd.Categorical(samples),
        "group": pd.Categorical(group),
        "cell_type": pd.Categorical(label, categories=TYPES + ["EXCLUDED", "UNRESOLVED"]),
    }, index=[f"c{i}" for i in range(n)])


def shape():
    A = ad.AnnData(X=L.copy(), obs=obs_frame(), var=pd.DataFrame(index=genes))
    A.layers["counts"] = C.copy()
    A.layers["lognorm"] = L.copy()
    A.obsm["X_harmony"] = emb(10)
    A.obsm["X_umap_harmony"] = emb(2)
    A.obsm["X_scanvi"] = emb(10)
    A.obsm["X_umap_scanvi"] = emb(2)
    A.obsm["X_umap"] = A.obsm["X_umap_scanvi"].copy()
    A.uns["scintegrate"] = {"default_method": "scanvi", "default_embedding": "X_scanvi",
                            "provenance": provenance()}
    return A


def write(A, name):
    A.write_h5ad(os.path.join(a.out, name), compression="gzip")


os.makedirs(a.out, exist_ok=True)
write(shape(), "facts-a.h5ad")

B = shape()
del B.uns["scintegrate"]
B.uns["integrator"] = {"provenance": provenance()}
B.obs["predicted_doublet"] = pd.Categorical(np.where(np.arange(n) % 9 == 0, "True", "False"))
B.obs["scrublet_score"] = rng.uniform(size=n).astype(np.float32)
B.obs["percent_mito_estimate"] = rng.uniform(0, 5, size=n).astype(np.float32)
B.uns["cellbender"] = {"note": "a key named like a tool, recording nothing"}
B.uns["log1p_cp10k"] = "yes"
B.layers["counts"] = L.copy()
B.obsm["X_pca"] = emb(2)
B.var["gene_ids"] = [f"SYM{i}" for i in range(g)]
write(B, "facts-b-decoy.h5ad")

E = shape()
E.uns["log1p"] = {"base": 2}
E.uns["qc_stage"] = {"doublet_method": "scrublet", "ambient_method": "cellbender"}
E.obs["total_counts"] = counts.sum(1).astype(np.float32)
# A recorded value, carried as the object holds it (wide enough to be read as numbers, not levels).
E.obs["n_genes_by_counts"] = rng.integers(200, 5000, size=n).astype(np.float32)
E.obs["pct_counts_mt"] = rng.uniform(0, 5, size=n).astype(np.float32)
write(E, "facts-e-records.h5ad")

F = shape()
F.uns["other_stage"] = {"default_method": "harmony", "default_embedding": "X_harmony"}
write(F, "facts-f-two.h5ad")

G = shape()
del G.layers["counts"]
scale = rng.uniform(2e3, 2e4, size=(n, 1))
G.X = sp.csr_matrix(np.log1p(counts / tot * scale).astype(np.float32))
G.layers["lognorm"] = G.X.copy()
G.uns["normalization"] = "CP10K"
write(G, "facts-g-lognorm.h5ad")

H = shape()
del H.uns["scintegrate"]
H.uns["as_group"] = {"default_method": {"name": "scanvi"}}
H.uns["as_array"] = {"default_method": np.array(["scanvi", "harmony"])}
H.uns["as_number"] = {"default_method": 3}
H.uns["as_long"] = {"default_method": "s" * 10_000}
H.uns["as_lines"] = {"default_method": "scanvi\nharmony"}
H.uns["as_markup"] = {"default_method": "<b>scanvi</b>"}
write(H, "facts-h-hostile.h5ad")
print(f"wrote 6 objects of {n} cells x {g} genes to {a.out}")
