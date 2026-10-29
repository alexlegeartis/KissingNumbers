"""Complete exact geometry of integer kissing configurations, with no input assumption.

Small integer dot products are evaluated exactly in binary32 or binary64.
Every product and every partial sum is bounded below the exact-integer limit.
The two signs of a line are grouped only when both actually occur.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical_hash(points):
    return hashlib.sha256(np.asarray(points, dtype='<i2').tobytes()).hexdigest()


def read_points(path):
    with np.load(path, allow_pickle=False) as z:
        require(set(z.files) == {'points', 'squared_norm', 'inner_product_limit',
                                 'coordinate_denominator'}, 'Archive fields')
        require(int(z['squared_norm']) == 128 and int(z['inner_product_limit']) == 64,
                'Normalization metadata')
        require(int(z['coordinate_denominator']) == 2, 'Coordinate denominator')
        points = z['points']
        require(points.dtype.kind in 'iu' and points.ndim == 2, 'Integer matrix')
        require(27 <= points.shape[1] <= 31 and len(points) >= 2, 'Dimensions')
        require(np.all(points >= -16) and np.all(points <= 16), 'Coefficient bounds')
    return points.astype(np.int16)


def verify(path, precision=32, block_size=256):
    started = time.monotonic()
    x = read_points(path)
    n, d = x.shape
    require(np.all(np.sum(x.astype(np.int64)**2, axis=1) == 128), 'Squared norms')
    require(len(np.unique(x, axis=0)) == n, 'Distinct points')
    # Products, intermediate sums and final sums are integers of magnitude
    # at most d*16^2 < 2^24, so binary32 and binary64 represent them exactly.
    require(d * 16**2 < 2**24, 'Exact BLAS arithmetic bound')
    require(precision in (32, 64), 'Arithmetic choice')
    dtype = np.float32 if precision == 32 else np.float64
    canonical = x.copy()
    first = np.argmax(canonical != 0, axis=1)
    sign = np.where(canonical[np.arange(n), first] < 0, -1, 1)
    canonical *= sign[:, None]
    representatives, inverse, counts = np.unique(
        canonical, axis=0, return_inverse=True, return_counts=True)
    require(np.all((counts == 1) | (counts == 2)), 'Line multiplicities')
    paired = representatives[counts == 2].astype(dtype)
    unpaired = x[counts[inverse] == 1].astype(dtype)
    p, u = len(paired), len(unpaired)
    require(2*p + u == n, 'Actual opposite-pair decomposition')
    line_checks = 0
    maximum = -128
    for at in range(0, p, block_size):
        block = paired[at:at+block_size]
        g = block @ paired[at:].T
        i, j = np.tril_indices(len(block))
        g[i, j] = 0
        require(np.all(np.abs(g) <= 64), 'Opposite-pair angle')
        maximum = max(maximum, int(np.abs(g).max(initial=0)))
        line_checks += len(block)*(p-at) - len(block)*(len(block)+1)//2
    for at in range(0, p, 1024):
        g = paired[at:at+1024] @ unpaired.T
        require(np.all(np.abs(g) <= 64), 'Opposite pair versus unpaired point')
        maximum = max(maximum, int(np.abs(g).max(initial=0)))
    for at in range(0, u, block_size):
        block = unpaired[at:at+block_size]
        g = block @ unpaired[at:].T
        i, j = np.tril_indices(len(block))
        g[i, j] = -128
        require(np.all(g <= 64), 'Unpaired-point angle')
        maximum = max(maximum, int(g.max(initial=-128)))
    require(line_checks == p*(p-1)//2, 'Line-pair coverage')
    covered = 4*line_checks + p + 2*p*u + u*(u-1)//2
    require(covered == n*(n-1)//2, 'Complete point-pair coverage')
    return dict(status='PASS', dimension=d, points=n, squared_norm=128,
                inner_product_limit=64, antipodal_pairs=p, unpaired_points=u,
                distinct_line_pairs_checked=line_checks,
                paired_to_unpaired_checks=p*u, unpaired_pairs_checked=u*(u-1)//2,
                unordered_point_pairs_covered=covered, maximum_checked_product=maximum,
                canonical_point_sha256=canonical_hash(x),
                file_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                arithmetic='Exact bounded integer products in binary%d; no tolerance' % precision,
                seconds=time.monotonic()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coordinates', type=Path)
    parser.add_argument('--precision', type=int, choices=(32, 64), default=32)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = verify(args.coordinates, args.precision)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
