"""Independent direct-integer cross-check of every inserted point.

This supplements the complete geometry replay: it does not establish the
angles between two unchanged points by itself. No BLAS products are used.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from verify_geometry import read_points, require, canonical_hash

ROOT = Path(__file__).resolve().parent


def replay(record):
    x = read_points(ROOT/record['certificate'])
    with np.load(ROOT/record['delta'], allow_pickle=False) as z:
        removed, added = z['removed'], z['added']
    d = x.shape[1]
    for a in (removed, added):
        require(a.dtype.kind in 'iu' and a.ndim == 2 and a.shape[1] == d, 'Delta format')
        require(np.all(a >= -16) and np.all(a <= 16), 'Delta coefficient bound')
        require(len(np.unique(a, axis=0)) == len(a), 'Distinct delta points')
        require(np.all(np.sum(a.astype(np.int64)**2, axis=1) == 128), 'Delta norms')
    lookup = {row.tobytes(): i for i, row in enumerate(x)}
    require(all(row.astype(np.int16).tobytes() in lookup for row in added), 'Added membership')
    require(all(row.astype(np.int16).tobytes() not in lookup for row in removed), 'Removed membership')
    added_set = {row.astype(np.int16).tobytes() for row in added}
    retained = np.array([row for row in x if row.tobytes() not in added_set], dtype=np.int16)
    baseline = np.unique(np.concatenate((retained, removed.astype(np.int16))), axis=0)
    require(len(baseline) == record['packingstar'], 'Reconstructed baseline count')
    require(canonical_hash(baseline) == record['source_canonical_sha256'], 'Baseline binding')
    require(len(x) == len(baseline)-len(removed)+len(added) == record['total'], 'Delta count')
    xi = x.astype(np.int64)
    checks = 0
    maximum = -128
    for row in added.astype(np.int64):
        dots = np.zeros(len(x), dtype=np.int64)
        for j in np.flatnonzero(row):
            dots += xi[:, j] * int(row[j])
        dots[lookup[row.astype(np.int16).tobytes()]] = -128
        require(int(dots.max()) <= 64, 'Inserted-point angle')
        maximum = max(maximum, int(dots.max()))
        checks += len(x)-1
    return dict(status='PASS', dimension=d, removed=len(removed), added=len(added),
                gain=len(added)-len(removed), direct_integer_checks=checks,
                maximum_inserted_product=maximum,
                canonical_point_sha256=canonical_hash(x))


if __name__ == '__main__':
    wanted = set(map(int, sys.argv[1:]))
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        if not wanted or record['dimension'] in wanted:
            print(json.dumps(replay(record)))
