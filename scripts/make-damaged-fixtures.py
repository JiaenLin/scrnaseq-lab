"""Write three fixtures whose reads fail inside HDF5, made by h5py (the real writer). Every value is synthetic.

- scripts/fixtures/damaged-reads.h5ad: family-2x2.h5ad with obsm/X_umap, X/data and obs/sample/codes rewritten
  gzip-compressed in small chunks, 8 bytes inside one chunk of each overwritten, and three embeddings whose stored
  type is not numbers: obsm/X_opaque (X_umap_none's values as opaque 4-byte items), obsm/X_compound (a compound of
  two floats) and obsm/X_null (an empty dataspace).
- scripts/fixtures/damaged-half.h5ad: nullable-mask.h5ad with obsm/X_umap, the embedding a guess takes first, stored
  as 16-bit floats, and obsm/X_tsne holding the same values as 32-bit floats.
- scripts/fixtures/damaged-only-matrix.h5ad: nullable-mask.h5ad with X/data, its only matrix, damaged the same way.

h5wasm hands back its buffer when such a read fails, so these are what the converter's one read is tested on.
Each damaged dataset is checked here to be unreadable by h5py itself.

    uv run --no-project --with h5py==3.14.0 --with numpy python scripts/make-damaged-fixtures.py
"""
import os
import shutil

import h5py
import numpy as np

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')


def at(name):
    return os.path.join(HERE, name)


def rewrite(f, path, chunks):
    """The dataset at `path`, written again gzip-compressed in chunks of `chunks`, its attributes kept."""
    old = f[path]
    values, attrs = old[()], dict(old.attrs)
    del f[path]
    d = f.create_dataset(path, data=values, compression='gzip', chunks=chunks)
    for k, v in attrs.items():
        d.attrs[k] = v


def damage(file, path, chunk):
    """Overwrite 8 bytes in the middle of one stored chunk, then check that h5py cannot read the dataset."""
    with h5py.File(file, 'r') as f:
        info = f[path].id.get_chunk_info(chunk)
    with open(file, 'r+b') as fh:
        fh.seek(info.byte_offset + info.size // 2 - 4)
        fh.write(b'\xff' * 8)
    with h5py.File(file, 'r') as f:
        try:
            f[path][()]
        except OSError:
            return
    raise SystemExit(f'{file}: {path} still reads after the damage')


out = at('damaged-reads.h5ad')
shutil.copyfile(at('family-2x2.h5ad'), out)
with h5py.File(out, 'r+') as f:
    rewrite(f, 'obsm/X_umap', (60, 2))
    rewrite(f, 'X/data', (500,))
    rewrite(f, 'obs/sample/codes', (60,))
    none = f['obsm/X_umap_none'][()]
    f['obsm'].create_dataset('X_opaque', data=np.ascontiguousarray(none).view('V4'))
    f['obsm'].create_dataset('X_compound', data=np.array([tuple(r) for r in none], dtype=[('a', 'f4'), ('b', 'f4')]))
    f['obsm'].create_dataset('X_null', data=h5py.Empty('f4'))
damage(out, 'obsm/X_umap', 2)
damage(out, 'X/data', 4)
damage(out, 'obs/sample/codes', 1)

out = at('damaged-half.h5ad')
shutil.copyfile(at('nullable-mask.h5ad'), out)
with h5py.File(out, 'r+') as f:
    umap = f['obsm/X_umap'][()]
    del f['obsm/X_umap']
    f['obsm'].create_dataset('X_umap', data=umap.astype('f2'))
    f['obsm'].create_dataset('X_tsne', data=umap)

out = at('damaged-only-matrix.h5ad')
shutil.copyfile(at('nullable-mask.h5ad'), out)
with h5py.File(out, 'r+') as f:
    rewrite(f, 'X/data', (10,))
damage(out, 'X/data', 1)
print('wrote damaged-reads.h5ad, damaged-half.h5ad and damaged-only-matrix.h5ad')
