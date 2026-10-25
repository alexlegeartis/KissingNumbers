"""Independently verify a compact Gamma72 symmetry certificate with NumPy.

python cpu_certificate.py certificate.npz [report.json]

The seed is checked in int64. Each stored integral isometry is checked exactly
using binary32 products after checking all partial-sum bounds. Every retained
line is hashed using exact binary64 factored products, then all keys are sorted
to prove disjointness. A collision causes rejection, never a false PASS.
The published Gamma72 lattice minimum and shell size are parent inputs.
"""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np


def verify(path, report_path=None):
    if not __debug__:
        raise RuntimeError('Run the verifier with assertions enabled (without -O).')
    started = time.perf_counter()
    G = np.load(Path(__file__).with_name('gamma72_gram.npy')).astype(np.int64)
    with np.load(path, allow_pickle=False) as z:
        C, matrices, masks, sizes, r = (z[name] for name in ('base','matrices','keep','sizes','hash_r'))
    assert C.dtype == matrices.dtype == np.dtype('int8')
    assert C.ndim == 2 and C.shape[1] == 72 and len(C) > 0
    assert matrices.shape == (len(sizes), 72, 72)
    assert sizes.ndim == 1 and sizes.dtype.kind in 'iu' and np.all(sizes > 0)
    assert r.shape == (72, 2) and r.dtype == np.dtype('int64')
    assert np.all(r >= -(2**30)) and np.all(r <= 2**30)
    assert np.abs(C.astype(np.int64)).max() <= 24
    assert np.abs(matrices.astype(np.int16)).max() <= 24
    assert np.array_equal(G, G.T) and np.abs(G).max() <= 8
    Gi = np.rint(np.linalg.inv(G.astype(np.float64))).astype(np.int64)
    assert np.array_equal(Gi @ G, np.eye(72, dtype=np.int64)) and Gi.diagonal().max() == 72
    C64 = C.astype(np.int64)
    products = C64 @ G @ C64.T
    assert np.all(products.diagonal() == 8), 'seed norm is not 8'
    np.fill_diagonal(products, 0)
    assert np.abs(products).max() <= 2, 'seed angle bound violated'
    del products
    keep = np.unpackbits(masks, axis=1)[:, :len(C)].astype(bool)
    assert keep.shape == (len(sizes), len(C))
    assert np.array_equal(keep.sum(axis=1), sizes)
    Gf = G.astype(np.float32)
    Cf = C.astype(np.float64)
    Rf = r.astype(np.float64)
    retained = np.empty(int(sizes.sum()), dtype=np.int64)
    cursor = 0
    for start in range(0, len(sizes), 256):
        A = matrices[start:start + 256].astype(np.float32)
        AG = A @ Gf
        # With |A|<=24 and |G|<=8 the first product has integer partial
        # sums <=13,824. Checking |AG|<=8 gives the same bound on the second.
        assert np.abs(AG).max() <= 8
        assert np.all(AG @ A.transpose(0, 2, 1) == Gf), 'matrix is not a lattice isometry'
        AR = A.astype(np.float64) @ Rf
        assert np.abs(AR).max() <= 72 * 24 * 2**30
        # Compute many images as one matrix product, avoiding tiny GEMMs.
        H = Cf @ AR.transpose(1, 0, 2).reshape(72, -1)
        H = H.reshape(len(C), len(A), 2).transpose(1, 0, 2).astype(np.int64)
        flip = (H[:, :, 0] < 0) | ((H[:, :, 0] == 0) & (H[:, :, 1] < 0))
        H[flip] *= -1
        keys = ((H[:, :, 0] & 0xffffffff) << 32) | (H[:, :, 1] & 0xffffffff)
        selected = keys[keep[start:start + len(A)]]
        retained[cursor:cursor + len(selected)] = selected
        cursor += len(selected)
    assert cursor == len(retained)
    geometry_seconds = time.perf_counter() - started
    retained.sort()
    assert np.all(retained[1:] != retained[:-1]), 'repeated line or hash collision'
    report = dict(status='PASS', n_classes=len(sizes), base_lines=len(C), total_lines=cursor,
                  sizes_min=int(sizes.min()), sizes_max=int(sizes.max()),
                  geometry_seconds=geometry_seconds, seconds=time.perf_counter()-started,
                  certificate_bytes=Path(path).stat().st_size,
                  certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                  scope='exact seed norms and every angle; every integral isometry; global retained-line disjointness',
                  numpy=np.__version__)
    srt = np.sort(sizes)[::-1].astype(np.int64)
    report['largest_class_sums'] = {str(n):int(srt[:n].sum()) for n in (1782,2466,3556,5800,9240,16632) if n<=len(sizes)}
    if report_path:
        Path(report_path).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    verify(sys.argv[1], sys.argv[2] if len(sys.argv)>2 else None)
