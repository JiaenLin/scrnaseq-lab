# Developer script — its OUTPUT is committed, so `npm test` needs no Python.
#
# Writes the synthetic objects in scripts/fixtures/ that `scripts/test-cli.mjs`
# converts through the real command. Everything in them is made up: neutral
# factor names (factor_a, factor_b), levels a1/a2 and b1/b2, samples s1..s8,
# cell types TypeA..TypeD and made-up gene names. Nothing comes from a real
# object, and nothing here may ever be replaced by one.
#
# The shape is the one the single-cell family delivers as its final object —
# the layout scIntegrate's emitter writes:
#
#   X                    log1p(CP10K), CSR
#   layers['counts']     integer counts;  layers['lognorm'] the same values as X
#   obs                  the batch key, the label columns and the design factors only
#   obsm                 X_pca, X_<method> (latent), X_umap_<method>, and X_umap
#                        (a copy of the chosen method's view)
#   uns['scintegrate']   the benchmark, the choice, and the run's provenance —
#                        nested dicts, string lists, scalars and an empty list
#
# plus the annotator's two sentinel labels, EXCLUDED and UNRESOLVED, among the
# cell types, because a real final object has them.
#
# Five files:
#
#   family-2x2.h5ad          two crossed factors, 2 samples per arm, written the
#                            way the family's emitters write (classic strings)
#   family-2x2-default.h5ad  the same object written with anndata's DEFAULT
#                            settings. On pandas >= 3 that stores obs/var names
#                            and unique string columns as nullable-string-array
#                            groups (values + mask), which the reader used to
#                            crash on.
#   nullable-mask.h5ad       a tiny object whose string columns stay
#                            nullable-string-arrays WITH missing values, so the
#                            mask itself is exercised (anndata only writes those
#                            when string-to-categorical conversion is off).
#   scaled-only.h5ad         a tiny object whose only matrix is z-scored, which
#                            no bundle can be built from: the command's "cannot
#                            convert" exit.
#   uns-records.h5ad         a tiny object whose run records hold what a JSON
#                            copy has to take care over: a structured array with
#                            an int64 field (anndata's rec-array), which once
#                            crashed the conversion, and leaves a JSON copy
#                            cannot hold as they are (float16, an HDF5
#                            reference, uint64 past 2**53, NaN and infinities,
#                            a dataframe). h5wasm writes none of these, so the
#                            test cannot make them itself.
#
# Pinned to what the family runs: anndata 0.13.2, pandas 3.0.5. With uv:
#
#   uv run --no-project --with anndata==0.13.2 --with pandas==3.0.5 \
#       python scripts/make-family-fixture.py [OUTDIR]
#
# OUTDIR defaults to scripts/fixtures next to this file.

import sys
import warnings
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.simplefilter("ignore", FutureWarning)

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "fixtures"
out_dir.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(20261109)

# ---- the design: factor_a x factor_b, crossed, two samples per arm ----------
FACTOR_A = ["a1", "a2"]
FACTOR_B = ["b1", "b2"]
REPS = 2
CELLS_PER_SAMPLE = 30
samples = []  # (sample id, level of factor_a, level of factor_b)
for a in FACTOR_A:
    for b in FACTOR_B:
        for _ in range(REPS):
            samples.append((f"s{len(samples) + 1}", a, b))

n = CELLS_PER_SAMPLE * len(samples)
sample_of = np.repeat(np.arange(len(samples)), CELLS_PER_SAMPLE)

TYPES = ["TypeA", "TypeB", "TypeC", "TypeD"]
type_of = rng.integers(0, len(TYPES), size=n)
forced = np.array([TYPES[t] for t in type_of], dtype=object)
label = forced.copy()
# The annotator's sentinels, at fixed places so every regeneration has them.
label[np.arange(n) % 25 == 3] = "EXCLUDED"
label[np.arange(n) % 40 == 17] = "UNRESOLVED"

obs = pd.DataFrame(index=[f"cell{i:04d}-{samples[s][0]}" for i, s in enumerate(sample_of)])
obs["sample"] = pd.Categorical([samples[s][0] for s in sample_of], categories=[s[0] for s in samples])
obs["cell_type"] = pd.Categorical(label, categories=TYPES + ["EXCLUDED", "UNRESOLVED"])
obs["cell_type_forced"] = pd.Categorical(forced, categories=TYPES)
obs["factor_a"] = pd.Categorical([samples[s][1] for s in sample_of], categories=FACTOR_A)
obs["factor_b"] = pd.Categorical([samples[s][2] for s in sample_of], categories=FACTOR_B)

# ---- expression: counts with a type signal, then log1p(CP10K) ----------------
G = 48
genes = [f"Gene{j + 1:03d}" for j in range(G)]
mean = rng.gamma(0.4, 2.0, size=(len(TYPES), G))
counts = rng.poisson(mean[type_of] * rng.uniform(0.6, 1.4, size=(n, 1))).astype(np.float32)
counts[:, :4] += rng.poisson(1.0, size=(n, 4))
counts[counts.sum(axis=1) == 0, 0] = 1  # no empty cell
C = sp.csr_matrix(counts)
lib = np.asarray(C.sum(axis=1)).ravel()
L = sp.csr_matrix(sp.diags(1e4 / lib) @ C, dtype=np.float32)
L.data = np.log1p(L.data)

var = pd.DataFrame(index=genes)
var["gene_id"] = [f"ENSXXXG{j + 1:011d}" for j in range(G)]

A = ad.AnnData(X=L, obs=obs, var=var)
A.layers["counts"] = C
A.layers["lognorm"] = L.copy()

# ---- embeddings, as the integration stage lays them out --------------------
pca = rng.normal(size=(n, 10)).astype(np.float32)
A.obsm["X_pca"] = pca
A.obsm["X_none"] = pca.copy()
A.obsm["X_harmony"] = rng.normal(size=(n, 10)).astype(np.float32)
A.obsm["X_umap_none"] = rng.normal(size=(n, 2)).astype(np.float32)
A.obsm["X_umap_harmony"] = rng.normal(size=(n, 2)).astype(np.float32)
A.obsm["X_umap"] = A.obsm["X_umap_harmony"].copy()

# ---- uns: the stage's record, provenance included --------------------------
A.uns["scintegrate"] = {
    "default_embedding": "X_harmony",
    "default_method": "harmony",
    "embeddings": {"X_none": "synthetic none", "X_harmony": "synthetic harmony"},
    "uncorrected_baseline": "X_pca",
    "X_is": "log1p of library-size-normalised counts, normalised over all samples together",
    "counts_layer": "layers['counts'] - raw integer counts",
    "benchmark": {
        "aggregate": {"none": 0.5, "harmony": 0.6},
        "chosen": {"default": "harmony",
                   "ranked": '[{"method": "harmony", "total": 0.6}, {"method": "none", "total": 0.5}]'},
    },
    "assessment": {"summary": {"needed": True}},
    "constraint_on_use": "synthetic fixture",
    "provenance": {
        "tool": "scintegrate", "version": "0.0.0+fixture", "commit": "0" * 40,
        "state_version": 1, "argv": ["scintegrate", "integrate", "--h5ad", "joint.h5ad"],
        "python": "3.12", "job": "fixture", "wrapped_versions": {"anndata": "0.13.2"},
        "sees": {"none": [], "harmony": ["sample"]}, "input": "joint.h5ad",
        "batch_key": "sample", "label_key": "cell_type", "seed": 0,
    },
}


def write(adata, path, *, classic, **kw):
    prev = ad.settings.allow_write_nullable_strings
    if classic:
        ad.settings.allow_write_nullable_strings = False
    try:
        adata.write_h5ad(path, **kw)
    finally:
        ad.settings.allow_write_nullable_strings = prev
    print(f"wrote {path.name}: {adata.n_obs} cells x {adata.n_vars} genes, {path.stat().st_size} bytes")


write(A, out_dir / "family-2x2.h5ad", classic=True)
write(A, out_dir / "family-2x2-default.h5ad", classic=False)

# ---- the mask: string columns with missing values that stay strings --------
m = 8
tiny_obs = pd.DataFrame(index=[f"c{i}" for i in range(m)])
tiny_obs["group"] = pd.Categorical(["g1", "g2"] * (m // 2))
tiny_obs["note"] = pd.array(["n0", "n1", None, "n3", "n4", None, "n6", "n7"], dtype="string")
tiny_var = pd.DataFrame(index=[f"Gene{j + 1:03d}" for j in range(4)])
tiny_var["gene_id"] = pd.array(["ENSXXXG00000000001", None, "ENSXXXG00000000003",
                                "ENSXXXG00000000004"], dtype="string")
T = ad.AnnData(X=sp.csr_matrix(rng.poisson(2.0, size=(m, 4)).astype(np.float32)),
               obs=tiny_obs, var=tiny_var)
T.obsm["X_umap"] = rng.normal(size=(m, 2)).astype(np.float32)
write(T, out_dir / "nullable-mask.h5ad", classic=False, convert_strings_to_categoricals=False)

# ---- nothing to build from: the only matrix is scaled -----------------------
z = rng.normal(size=(m, 4)).astype(np.float32)
z -= z.mean(axis=0)
S = ad.AnnData(X=z, obs=pd.DataFrame({"cell_type": pd.Categorical(["TypeA", "TypeB"] * (m // 2))},
                                     index=[f"c{i}" for i in range(m)]),
               var=pd.DataFrame(index=[f"Gene{j + 1:03d}" for j in range(4)]))
S.obsm["X_umap"] = rng.normal(size=(m, 2)).astype(np.float32)
write(S, out_dir / "scaled-only.h5ad", classic=True)

# ---- run records holding what a JSON copy has to take care over --------------
# Written last, so the four files above keep their bytes.
R = ad.AnnData(X=sp.csr_matrix(rng.poisson(2.0, size=(m, 4)).astype(np.float32)),
               obs=pd.DataFrame({"cell_type": pd.Categorical(["TypeA", "TypeB"] * (m // 2))},
                                index=[f"c{i}" for i in range(m)]),
               var=pd.DataFrame(index=[f"Gene{j + 1:03d}" for j in range(4)]))
R.obsm["X_umap"] = rng.normal(size=(m, 2)).astype(np.float32)
R.uns["scintegrate"] = {"provenance": {
    "tool": "synthetic-tool", "version": "0.0.0+fixture",
    # anndata writes a structured array as a rec-array: one HDF5 compound
    # dataset, read back as rows whose int64 field is a BigInt.
    "ranked": np.array([(b"harmony", 3, 0.6), (b"none", 1, 0.5)],
                       dtype=[("method", "S8"), ("rank", "<i8"), ("score", "<f8")]),
}}
R.uns["odd_provenance"] = {
    "tool": "synthetic-tool",
    "half": np.arange(4, dtype=np.float16),
    "big": np.uint64(2**63 + 5),
    "special": np.array([np.nan, np.inf, -np.inf, 1.0]),
    "flag": np.bool_(False),
    "frame": pd.DataFrame({"a": [1, 2]}, index=["r1", "r2"]),
}
write(R, out_dir / "uns-records.h5ad", classic=True)
# anndata writes no HDF5 references; h5py adds one, to uns itself.
with h5py.File(out_dir / "uns-records.h5ad", "r+") as f:
    f["uns/odd_provenance"].create_dataset("ref", data=f["uns"].ref, dtype=h5py.ref_dtype, track_times=False)
print(f"added uns/odd_provenance/ref: uns-records.h5ad is {(out_dir / 'uns-records.h5ad').stat().st_size} bytes")
