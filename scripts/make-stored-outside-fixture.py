"""Write scripts/fixtures/stored-outside.h5ad: nullable-mask.h5ad with one more embedding,
obsm/X_stored_outside, whose raw values HDF5 keeps in another file (external storage,
stored-outside.bin beside it). h5wasm cannot write external storage, so this one is h5py's.
Every value is synthetic.

    uv run --no-project --with h5py==3.14.0 --with numpy python scripts/make-stored-outside-fixture.py
"""
import os
import shutil

import h5py
import numpy as np

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
SRC, OUT, RAW = (os.path.join(HERE, n) for n in ('nullable-mask.h5ad', 'stored-outside.h5ad', 'stored-outside.bin'))

shutil.copyfile(SRC, OUT)
with h5py.File(OUT, 'r+') as f:
    n = next(iter(f['obsm'].values())).shape[0]
    np.full((n, 2), 1.5, dtype=np.float32).tofile(RAW)
    f['obsm'].create_dataset('X_stored_outside', shape=(n, 2), dtype=np.float32,
                             external=[('stored-outside.bin', 0, n * 2 * 4)])
print('wrote', OUT, 'and', RAW)
