"""Irrational E8/cube module in dimension 31; discovery uses exact products.

At squared norm 32, h is an E8 norm-16 vector and the seven tail coordinates
are (4/sqrt(7))*s.  All comparisons with the integer incumbent reduce to
integer square roots.  Candidate/candidate products satisfy the kissing
condition iff 7*dot(h,h') <= 32*Hamming(s,s').
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

NORMAL = np.array([0, 2, 4, 5, 6, 10, 11, 23])


def e8_shell16():
    rows = []
    for i in range(8):
        for sign in [-1, 1]:
            h = np.zeros(8, dtype=np.int16)
            h[i] = 4*sign
            rows.append(h)
    for support in itertools.combinations(range(8), 4):
        for signs in itertools.product([-2, 2], repeat=4):
            h = np.zeros(8, dtype=np.int16)
            h[list(support)] = signs
            rows.append(h)
    for i in range(8):
        for signs in itertools.product([-1, 1], repeat=8):
            if sum(s < 0 for s in signs) % 2 != 1:
                continue
            h = np.array(signs, dtype=np.int16)
            h[i] *= 3
            rows.append(h)
    h = np.array(rows, dtype=np.int16)
    assert len(h) == len(np.unique(h, axis=0)) == 2160
    assert np.all(np.sum(h.astype(np.int32)**2, axis=1) == 16)
    assert np.all(h.sum(axis=1) % 4 == 0)
    return h


def integer_threshold(radical_numerator):
    # The largest integer a with a + b/sqrt(7) <= 16.
    b = int(radical_numerator)
    floor = math.isqrt(b*b//7)
    if b > 0:
        ceiling = floor + int(7*floor*floor != b*b)
        return 16 - ceiling
    return 16 + floor


def prepare(source):
    z = np.load(source, allow_pickle=False)
    assert int(z['squared_norm']) == 128
    assert np.all(z['points'] % 2 == 0)
    x = z['points'].astype(np.int16)//2
    assert len(x) == 251254 and x.shape[1] == 31
    assert np.all(np.sum(x.astype(np.int32)**2, axis=1) == 32)
    h = e8_shell16()
    s = np.array(list(itertools.product([-1, 1], repeat=7)), dtype=np.int16)
    q = np.sum(x[:, 24:].astype(np.int32)**2, axis=1)
    p = np.sum(x[:, NORMAL].astype(np.int32)**2, axis=1)
    relevant = ((q == 8) & (p > 8)) | ((q == 16) & (p > 0))
    ids = np.flatnonzero(relevant)
    assert len(ids) == 2136
    # Every omitted old point is covered by one of these direct exact guards.
    eq = np.unique(x[q == 0][:, NORMAL], axis=0)
    assert np.max(h.astype(np.int32) @ eq.T.astype(np.int32)) <= 16
    mids = x[(q == 8) & ~relevant]
    assert np.max(h.astype(np.int32) @ np.unique(mids[:, NORMAL], axis=0).T.astype(np.int32)) <= 8
    assert np.max(np.abs(s.astype(np.int32) @ np.unique(mids[:, 24:], axis=0).T.astype(np.int32))) <= 4
    # (4/sqrt7)*4 < 8, so midpoint sums are at most 16.
    assert 16*16 < 7*8*8
    pure = x[q == 32, 24:]
    assert np.max(np.abs(s.astype(np.int32) @ pure.T.astype(np.int32))) <= 8
    assert 32*32 < 7*16*16
    heads0 = x[(q == 16) & ~relevant, 24:]
    assert np.max(np.abs(s.astype(np.int32) @ heads0.T.astype(np.int32))) <= 8
    assert np.all(np.isin(q, [0, 8, 16, 32]))
    a = h.astype(np.int32) @ x[ids][:, NORMAL].T.astype(np.int32)
    b = 4*s.astype(np.int32) @ x[ids, 24:].T.astype(np.int32)
    thresholds = np.array([integer_threshold(v) for v in range(-128, 129)], dtype=np.int16)
    assert b.min() >= -128 and b.max() <= 128
    cost = np.empty((len(s), len(h)), dtype=np.int16)
    for j in range(len(s)):
        cost[j] = np.sum(a > thresholds[b[j]+128], axis=1)
    return x, h, s, ids, a, b, thresholds, cost


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cutoff', type=int, default=8)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    x, h, s, ids, a, b, thresholds, cost = prepare(args.source)
    ts, hs = np.where(cost <= args.cutoff)
    offsets = [0]
    adjacency = []
    for t, i in zip(ts, hs):
        blockers = ids[np.flatnonzero(a[i] > thresholds[b[t]+128])]
        adjacency.extend(blockers.tolist())
        offsets.append(len(adjacency))
    np.savez_compressed(args.out/'cube7_pool.npz', body=h[hs], signs=s[ts],
                        body_index=hs, tail_index=ts, blocker_offsets=np.array(offsets),
                        blockers=np.array(adjacency), all_costs=cost)
    report = {
        'source': str(args.source),
        'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'baseline': len(x), 'dimension': 31,
        'normal_octad': NORMAL.tolist(), 'core_vectors': len(h),
        'tail_vectors': len(s), 'candidate_count': int(cost.size),
        'relevant_old_points': len(ids),
        'minimum_blockers': int(cost.min()),
        'cost_histogram': {str(int(k)):int(v) for k,v in zip(*np.unique(cost,return_counts=True))},
        'cutoff': args.cutoff, 'pool_size':len(ts),
        'incidences':len(adjacency),
        'seconds':time.monotonic()-started,
        'scope':'All 276480 E8 norm16/cube7 norm16 candidates are priced exactly; no selection result follows from pricing.'
    }
    (args.out/'pricing.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
