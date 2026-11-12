"""Independent direct-int64 changed-pair replay and source reconstruction."""
from pathlib import Path
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parent


def require(value, label):
    if not value:
        raise ValueError(label)


def replay(record):
    with np.load(ROOT/record['certificate'], allow_pickle=False) as z:
        final = z['points'].astype(np.int64)
    with np.load(ROOT/record['exchange'], allow_pickle=False) as z:
        added = z['added'].astype(np.int64)
        removed = z['removed'].astype(np.int64)
    require(len(added) == record['added'] and len(removed) == record['removed'], 'Delta counts')
    fkeys = {tuple(map(int, row)): i for i, row in enumerate(final)}
    akeys = {tuple(map(int, row)) for row in added}
    rkeys = {tuple(map(int, row)) for row in removed}
    require(len(fkeys) == len(final) and len(akeys) == len(added)
            and len(rkeys) == len(removed), 'Distinct delta and final rows')
    require(akeys <= fkeys.keys() and not (rkeys & fkeys.keys()), 'Delta membership')
    baseline = np.asarray(sorted((fkeys.keys() - akeys) | rkeys), dtype='<i2')
    require(len(baseline) == record['packingstar'], 'Source count')
    require(hashlib.sha256(baseline.tobytes()).hexdigest()
            == record['source_canonical_sha256'], 'Reconstructed pinned source hash')
    require(len(final) == len(baseline) + len(added) - len(removed), 'Exchange arithmetic')
    require(np.all(np.sum(final*final, axis=1) == 128), 'Final norms by integer sums')
    max_product = -128
    comparisons = 0
    for point in added:
        dots = np.zeros(len(final), dtype=np.int64)
        for j in range(final.shape[1]):
            dots += int(point[j])*final[:, j]
        own = fkeys[tuple(map(int, point))]
        require(dots[own] == 128, 'Inserted self product')
        dots[own] = -128
        require(np.all(dots <= 64), 'Direct-integer inserted angle')
        max_product = max(max_product, int(dots.max()))
        comparisons += len(final)-1
    return dict(status='PASS', dimension=record['dimension'], added=len(added),
                removed=len(removed), total=len(final),
                changed_point_comparisons=comparisons,
                maximum_changed_product=max_product,
                arithmetic='Direct int64 coordinate sums; no BLAS or tolerance',
                source_canonical_sha256=record['source_canonical_sha256'])
