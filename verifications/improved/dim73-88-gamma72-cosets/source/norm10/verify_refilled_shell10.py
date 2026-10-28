"""Verify the enlarged explicit 73D norm-ten layer in two exact arithmetics.

All head/head and mixed products are checked in int64, then recomputed by
binary32 GEMM after proving that every term and partial sum is integral
and smaller than 2**24. No search pool or conflict graph is trusted.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent


def certificate():
    parent = HERE.parent/'dim73_surgery'
    with np.load(parent/'lattice/certificate.npz', allow_pickle=False) as data:
        g = data['gamma_gram'].astype(np.int64)
        x = data['heads'].astype(np.int64)
    u = np.load(parent/'surgery_class.npy', allow_pickle=False).astype(np.int64)
    path = HERE/'refilled_shell10.npz'
    with np.load(path, allow_pickle=False) as data:
        raw = data['heads']
        assert raw.dtype == np.int16 and raw.ndim == 2 and raw.shape[1] == 72
        y = raw.astype(np.int64)
    assert len(y) > len(np.load(HERE/'shell10_selected.npy'))
    assert np.max(np.abs(y)) <= 64
    assert len(np.unique(y, axis=0)) == len(y)
    yg = y@g
    assert np.all(np.sum(yg*y, axis=1) == 10)
    assert max(int(np.max(np.sum(np.abs(yg), axis=1)))*int(np.max(np.abs(v)))
               for v in (y, x, u)) < 2**24
    limits = []
    for source in (y, x, u):
        source32 = source.T.astype(np.float32)
        largest, smallest = -100, 100
        for i in range(0, len(y), 128):
            integer = yg[i:i+128]@source.T
            binary = yg[i:i+128].astype(np.float32)@source32
            assert np.array_equal(integer, binary)
            if source is y:
                integer[np.arange(len(integer)), i+np.arange(len(integer))] = 0
            largest = max(largest, int(integer.max()))
            smallest = min(smallest, int(integer.min()))
        limits.append((smallest, largest))
    assert -6 <= limits[0][0] and limits[0][1] <= 3
    assert max(map(abs, limits[1])) <= 5 and max(map(abs, limits[2])) <= 3
    # The same norm-eight points as in the preceding explicit layer:
    # +/-(4y/5,sqrt(8/5)). The inequalities below include their antipodes.
    assert -100 <= 16*limits[0][0]+40 <= 16*limits[0][1]+40 <= 100
    assert 2*max(map(abs, limits[1]))+10 <= 20
    assert 400-108-80 > 0 and 48**2*15 < (400-108-80)**2
    previous = 6218175600+2*len(x)+2*len(u)+2*len(np.load(HERE/'shell10_selected.npy'))
    return {
        'status': 'PASS', 'dimension': 73,
        'oriented_heads': len(y), 'added_layer_points': 2*len(y),
        'previous_bound': previous,
        'bound': 6218175600+2*len(x)+2*len(u)+2*len(y),
        'gain_over_previous_bound': 2*(len(y)-len(np.load(HERE/'shell10_selected.npy'))),
        'head_product_interval': list(limits[0]),
        'coset_product_interval': list(limits[1]),
        'cap_product_interval': list(limits[2]),
        'products_checked_in_each_arithmetic': len(y)*(len(y)+len(x)+len(u)),
        'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'scope': 'Explicit layer; parent surgery/lattice certificate and published Gamma72 minimum inherited.',
    }


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    report = certificate()
    path = HERE/'refilled_shell10_verification.json'
    if '--check' in sys.argv:
        assert json.loads(path.read_text()) == report, 'stale certificate'
    else:
        path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
