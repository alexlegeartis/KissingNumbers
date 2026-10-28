"""Verify the complete additional weight-eight layer in dimension 61."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dim49plus'))
from latitude_compression import quad_leq


def packed_check(masks):
    assert all(0 <= m < 2**61 and m.bit_count() == 8 for m in masks)
    assert len(set(masks)) == len(masks)
    assert all(sum(bool(m & (15 << (4*r))) for r in range(15)) <= 5 for m in masks)
    maximum = 0
    for i, a in enumerate(masks):
        for b in masks[:i]:
            overlap = (a & b).bit_count()
            assert overlap <= 4
            maximum = max(maximum, overlap)
    return maximum


def coordinate_check(masks):
    assert all(0 <= m < 2**61 for m in masks)
    incidence = np.array([[(m >> c) & 1 for c in range(61)] for m in masks], dtype=np.float32)
    assert len(np.unique(incidence, axis=0)) == len(masks)
    assert np.all(incidence.sum(axis=1) == 8)
    assert np.all(np.count_nonzero(incidence[:, :60].reshape(-1, 15, 4).sum(axis=2), axis=1) <= 5)
    # Every partial sum is an integer between zero and 61, exactly represented.
    with threadpool_limits(limits=3):
        for start in range(0, len(masks), 256):
            end = min(start+256, len(masks))
            dots = incidence[start:end] @ incidence[:end].T
            dots[np.arange(end-start), start+np.arange(end-start)] = 0
            assert np.all(dots <= 4)


def geometry():
    signs = [x for x in range(256) if x.bit_count() % 2 == 0]
    assert len(signs) == 128
    assert min((a ^ b).bit_count() for i, a in enumerate(signs) for b in signs[:i]) == 2
    assert Q(8-2*2, 8) == Q(1, 2)
    assert Q(4, 8) == Q(1, 2)
    assert Q(25, 8*15) <= Q(1, 4)  # at most five transversal rows
    assert Q(4, 8*2) == Q(1, 4)    # all existing weight-two points
    assert Q(8, 64) <= Q(1, 4) and Q(8*17, 1024) <= Q(1, 4)
    assert Q(81, 512) <= Q(1, 4)   # j=0 fibre, support contains coordinate 60
    assert quad_leq(Q(837, 8192), Q(28, 8192), 17, Q(1, 4))
    return len(signs)


def negative_controls():
    controls = [[255, 255], [511], [255, 383],
                [sum(1 << c for c in [0, 1, 4, 5, 8, 12, 16, 20])]]
    for example in controls:
        for checker in (packed_check, coordinate_check):
            try:
                checker(example)
            except AssertionError:
                continue
            raise AssertionError('A malformed support certificate was accepted')
    return len(controls)*2


def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are required.')
    path = HERE / 'bridge61.npz'
    with np.load(path, allow_pickle=False) as data:
        supports = data['supports']
    assert supports.dtype == np.uint64 and supports.ndim == 1
    masks = list(map(int, supports))
    maximum = packed_check(masks)
    coordinate_check(masks)
    sign_count = geometry()
    controls = negative_controls()
    parent = json.loads((HERE.parent / 'dim49plus/latitude_compression.json').read_text())['counts']
    assert parent['total'] == parent['dense']+parent['weight15']+parent['weight2']
    report = dict(status='Verified: exact packed and coordinate support replays and all cross-layer inequalities',
                  dimension=61, supports=len(masks), signs_per_support=sign_count,
                  support_pairs_checked_twice=len(masks)*(len(masks)-1)//2,
                  maximum_intersection=maximum, maximum_occupied_rows=5,
                  negative_controls_rejected=controls, previous_bound=parent['total'],
                  added_points=len(masks)*sign_count,
                  bound=parent['total']+len(masks)*sign_count,
                  certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  geometry='Unit vectors with eight nonzero coordinates +/-1/sqrt(8); even sign parity. The old families have respectively 61, 15 and 2 nonzero coordinates, proving disjointness.')
    (HERE / 'bridge61.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
