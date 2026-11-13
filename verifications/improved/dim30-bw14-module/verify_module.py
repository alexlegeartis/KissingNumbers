"""Independently reconstruct the BW14 module and the exact source exchange.

This checker uses integer coordinate sums for all changed pairs. The complete
configuration is checked separately by verify_geometry.py, without assuming
the validity of the external source.
"""
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
NORMAL = (0, 2, 4, 5, 6, 10, 11, 23)


def require(value, label):
    if not value:
        raise ValueError(label)


def point_keys(x):
    return {tuple(map(int, row)) for row in x}


def reconstruct():
    points = []
    for i, j in combinations(range(14), 2):
        for a, b in product((-4, 4), repeat=2):
            row = [0] * 14
            row[i], row[j] = a, b
            points.append(row)
    supports = []
    for normal in range(1, 16):
        for parity in (0, 1):
            support = [i for i in range(16)
                       if (normal & i).bit_count() % 2 == parity]
            if 14 in support or 15 in support:
                continue
            supports.append(support)
            for signs in product((-2, 2), repeat=8):
                if sum(s < 0 for s in signs) % 2:
                    continue
                row = [0] * 14
                for j, s in zip(support, signs):
                    row[j] = s
                points.append(row)
    require(len(supports) == 7 and len(points) == 1260, 'BW14 construction count')
    return np.asarray(points, dtype=np.int64)


def replay():
    provenance = json.loads((ROOT / 'data/provenance.json').read_text())
    full = reconstruct()
    require(len(point_keys(full)) == 1260, 'Module distinctness')
    require(np.all(np.sum(full * full, axis=1) == 32), 'Module norms')
    dot = full @ full.T
    np.fill_diagonal(dot, -32)
    require(np.max(dot) == 16, 'Complete direct-integer BW14 geometry')
    tail_norm = np.sum(full[:, 8:] ** 2, axis=1)
    base = full[tail_norm == 16]
    require(base.shape == (960, 14), 'Module head count')
    with np.load(ROOT / 'data/module.npz', allow_pickle=False) as z:
        r = z['transformation_numerator'].astype(np.int64)
        denominator = int(z['transformation_denominator'])
        heads = z['heads'].astype(np.int64)
        require(np.array_equal(z['base_heads'], base), 'Stored base versus reconstruction')
    require(denominator == 4 and r.shape == (14, 14), 'Transformation format')
    require(np.array_equal(r @ r.T, 16 * np.eye(14, dtype=np.int64)), 'Exact isometry')
    numerator = base @ r.T
    require(not np.any(numerator % 4), 'Integral transformed module')
    require(np.array_equal(numerator // 4, heads), 'Transformed heads')
    with np.load(ROOT / 'data/coordinates30_248020.npz', allow_pickle=False) as z:
        final128 = z['points'].astype(np.int64)
    require(final128.shape == (248020, 30) and not np.any(final128 % 2), 'Final format')
    final = final128 // 2
    with np.load(ROOT / 'data/exchange.npz', allow_pickle=False) as z:
        added = z['added'].astype(np.int64)
        removed = z['removed'].astype(np.int64)
    expected = np.zeros((960, 30), dtype=np.int64)
    expected[:, NORMAL], expected[:, 24:] = heads[:, :8], heads[:, 8:]
    require(point_keys(expected) == point_keys(added), 'Exact inserted module')
    require(added.shape == (960, 30) and removed.shape == (932, 30), 'Exchange counts')
    fk, ak, rk = point_keys(final), point_keys(added), point_keys(removed)
    require(len(fk) == len(final) and len(ak) == 960 and len(rk) == 932, 'Distinct exchange')
    require(ak <= fk and not (rk & fk), 'Exchange membership')
    source = np.asarray(sorted((fk - ak) | rk), dtype='<i2')
    require(len(source) == 247992, 'Reconstructed source count')
    source_hash = hashlib.sha256(source.tobytes()).hexdigest()
    require(source_hash == provenance['source_canonical_norm32_sha256'], 'Pinned source identity')
    outside = [j for j in range(24) if j not in NORMAL]
    require(not np.any(removed[:, outside]), 'Only normal module points removed')
    require(np.all(np.sum(removed[:, 24:] ** 2, axis=1) == 16), 'Only old heads removed')
    require(np.all(np.sum(final * final, axis=1) == 32), 'All final norms by integer sums')
    lookup = {row: i for i, row in enumerate(map(tuple, final))}
    maximum = -32
    comparisons = 0
    for point in added:
        dots = np.zeros(len(final), dtype=np.int64)
        for j in np.flatnonzero(point):
            dots += int(point[j]) * final[:, j]
        own = lookup[tuple(point)]
        require(dots[own] == 32, 'Inserted self product')
        dots[own] = -32
        require(np.max(dots) <= 16, 'Every changed pair by direct integer sums')
        maximum = max(maximum, int(np.max(dots)))
        comparisons += len(final) - 1
    return dict(status='PASS', dimension=30, source_count=247992, total=len(final),
                removed=len(removed), added=len(added), gain=28,
                source_canonical_norm32_sha256=source_hash,
                reconstructed_module_count=len(full), reconstructed_heads=len(heads),
                changed_point_comparisons=comparisons, maximum_changed_product=maximum,
                arithmetic='Direct int64 coordinate sums; no tolerance or floating matrix products')


if __name__ == '__main__':
    print(json.dumps(replay(), indent=2))
