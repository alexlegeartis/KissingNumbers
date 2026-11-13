"""Exact null-axis completions of the 28-point residual bulk fibers.

The source is split into eight normal, sixteen complementary, and seven tail
coordinates. A fixed norm-16 complementary vector leaves a norm-16 fiber in
R15 with pairwise nonpositive products. This script checks the actual rank-14
orthoplexes and prices their uniquely determined missing axes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research/reclaim_26_31_20261112/modules'))
from contact_refill import select


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cutoff', type=int, default=8)
    ap.add_argument('--seconds', type=float, default=90)
    ap.add_argument('--skip-selection', action='store_true',
                    help='Save the exact pool without the older unstrengthened MILP.')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    z = np.load(args.source, allow_pickle=False)
    X = z['points'].astype(np.int16)
    if int(z['squared_norm']) == 128:
        if np.any(X % 2):
            raise ValueError('Norm32 integrality')
        X //= 2
    if not np.all(np.sum(X.astype(np.int64)**2, axis=1) == 32):
        raise ValueError('Source norms')
    I = np.array([0, 2, 4, 5, 6, 10, 11, 23])
    J = np.array([j for j in range(24) if j not in I])
    A = np.r_[I, np.arange(24, 31)]
    ids = np.flatnonzero(np.sum(X[:, J].astype(np.int64)**2, axis=1) == 16)
    B, inverse, counts = np.unique(X[ids][:, J], axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind='stable')
    ptr = np.r_[0, np.cumsum(counts)]
    candidates = []
    radical = []
    shapes = {}
    for index in np.flatnonzero(counts == 28):
        f = X[ids[order[ptr[index]:ptr[index+1]]]][:, A].astype(np.int64)
        gram = f @ f.T
        if set(np.unique(gram)) != {-16, 0, 16} or np.sum(gram == -16) != 28:
            raise ValueError('Expected actual fourteen antipodal orthogonal pairs')
        numerator = 32 * np.eye(15, dtype=np.int64) - f.T @ f
        if np.any(numerator % 2):
            raise ValueError('Integral rank-one projector numerator')
        Q = numerator // 2
        if np.trace(Q) != 16 or not np.array_equal(Q @ Q, 16 * Q):
            raise ValueError('Exact rank-one complement')
        j = int(np.argmax(np.diag(Q)))
        q = int(Q[j, j])
        sq = math.isqrt(q)
        if sq * sq != q or np.any(Q[:, j] % sq):
            radical.append(dict(body=B[index].tolist(), numerator=Q[:, j].tolist(), squared_denominator=q))
            continue
        h = Q[:, j] // sq
        if np.dot(h, h) != 16 or np.any(f @ h):
            raise ValueError('Exact missing axis')
        key = ','.join(map(str, (int(np.dot(h[:8], h[:8])), int(np.dot(h[8:], h[8:])))))
        shapes[key] = shapes.get(key, 0) + 1
        for sign in (-1, 1):
            y = np.zeros(31, dtype=np.int16)
            y[J] = B[index]
            y[A] = sign * h
            candidates.append(y)
    C = np.unique(np.asarray(candidates, dtype=np.int16), axis=0)
    known = {row.tobytes() for row in X}
    C = np.asarray([row for row in C if row.tobytes() not in known], dtype=np.int16)
    print(json.dumps(dict(fibers=len(B), size28_fibers=int(np.sum(counts == 28)), candidates=len(C),
                          radical_axes=len(radical), missing_axis_shapes=shapes)), flush=True)
    F = X.astype(np.float32)
    blocker_counts = np.zeros(len(C), dtype=np.int32)
    active = np.arange(len(C))
    rng = np.random.default_rng(20261113)
    ordering = rng.permutation(len(X))
    for at in range(0, len(X), 512):
        if not len(active):
            break
        products = C[active].astype(np.float32) @ F[ordering[at:at+512]].T
        blocker_counts[active] += np.sum(products > 16, axis=1)
        active = active[blocker_counts[active] <= args.cutoff]
    admitted = C[active]
    blockers = []
    for y in admitted:
        dot = np.zeros(len(X), dtype=np.int64)
        for j in np.flatnonzero(y):
            dot += int(y[j]) * X[:, j].astype(np.int64)
        blockers.append(np.flatnonzero(dot > 16).tolist())
    ptr = np.r_[0, np.cumsum([len(b) for b in blockers])]
    flat = np.asarray([i for b in blockers for i in b], dtype=np.int64)
    np.savez_compressed(args.out / 'pool.npz', candidates=2*admitted,
                        blocker_indptr=ptr, blocker_indices=flat)
    selection = (dict(gain=0, status='not_run') if args.skip_selection else
                 select(2*admitted, blockers, args.seconds) if len(admitted) else dict(gain=0))
    if selection['gain'] > 0:
        keep = np.ones(len(X), dtype=bool)
        keep[selection['removed']] = False
        add = admitted[selection['added']]
        final = np.r_[X[keep], add]
        if len(np.unique(final, axis=0)) != len(final):
            raise ValueError('Final uniqueness')
        for i, y in enumerate(add):
            dot = final.astype(np.int64) @ y.astype(np.int64)
            dot[len(X)-len(selection['removed'])+i] = -32
            if np.max(dot) > 16:
                raise ValueError('Exact changed-pair replay')
        np.savez_compressed(args.out / f'coordinates31_{len(final)}.npz', points=2*final,
                            squared_norm=128, inner_product_limit=64, coordinate_denominator=2)
        np.savez_compressed(args.out / 'exchange.npz', added=2*add,
                            removed=2*X[selection['removed']])
        selection['changed_pairs'] = 'PASS; complete final proof still required'
    report = dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                  baseline=len(X), complementary_norm=16, residual_norm=16,
                  fiber_histogram={str(int(k)):int(v) for k,v in zip(*np.unique(counts, return_counts=True))},
                  exact_size28_orthoplexes=int(np.sum(counts == 28)), missing_axis_shapes=shapes,
                  new_integer_candidates=len(C), radical_axes=radical, cutoff=args.cutoff,
                  admitted=len(admitted), free_candidates=sum(not b for b in blockers),
                  selection=selection, seconds=time.monotonic()-started,
                  scope='Complete size28 null-axis generation; blocker-restricted selection only. No global capacity claim.')
    (args.out / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'radical_axes'}), flush=True)


if __name__ == '__main__':
    main()
