"""Bounded joint selection from the 759-octad 26D candidate pool.

Saved blocker indices are independently checked against the exact quadratic
source before optimization. A numerical optimizer only proposes a subset;
all changes in any positive witness are checked using integer comparisons.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[2]


def above16(a, b):
    c = a - 16
    return ((b == 0) & (c > 0)) | ((b > 0) & ((c >= 0) | (3*b*b > c*c))) | ((b < 0) & (c > 0) & (c*c > 3*b*b))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seconds', type=float, default=90)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source = ROOT.parent / 'PackingStar/25D-31D/new_bounds_20261111/D26_201180_coordinates.npz'
    pool = Path(__file__).resolve().parent / 'continuous26/transplant759.npz'
    z = np.load(source, allow_pickle=False)
    assert int(z['rational_denominator']) == 1
    assert np.all(z['quadratic_denominators'] == 1) and np.all(z['quadratic_radicands'] == 3)
    A = np.r_[z['rational_points'], z['quadratic_A']].astype(np.int64)
    B = np.r_[np.zeros_like(z['rational_points']), z['quadratic_B']].astype(np.int64)
    z = np.load(pool, allow_pickle=False)
    C = z['candidate_A'].astype(np.int64)
    assert not np.any(z['candidate_B']) and np.all(np.sum(C*C, axis=1) == 32)
    assert len(C) == len(np.unique(C, axis=0))
    ptr, idx = z['blocker_indptr'], z['blocker_indices']
    blockers = [idx[ptr[i]:ptr[i+1]] for i in range(len(C))]
    for i, point in enumerate(C):
        da = np.zeros(len(A), dtype=np.int64)
        db = np.zeros(len(A), dtype=np.int64)
        for j in np.flatnonzero(point):
            da += point[j]*A[:, j]
            db += point[j]*B[:, j]
        actual = np.flatnonzero(above16(da, db))
        assert np.array_equal(np.sort(blockers[i]), actual)
    print(json.dumps({'exact_blocker_replay': 'PASS', 'candidates': len(C)}), flush=True)
    old = np.unique(idx)
    mapping = {int(x): len(C)+i for i, x in enumerate(old)}
    rr, cc, vv, upper = [], [], [], []
    for i, bs in enumerate(blockers):
        for j in bs:
            row = len(upper)
            rr.extend((row, row)); cc.extend((i, mapping[int(j)])); vv.extend((1, -1)); upper.append(0)
    dots = C @ C.T
    ii, jj = np.where(np.triu(dots > 16, 1))
    for i, j in zip(ii, jj):
        row = len(upper)
        rr.extend((row, row)); cc.extend((int(i), int(j))); vv.extend((1, 1)); upper.append(1)
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(upper), len(C)+len(old))).tocsc()
    result = milp(np.r_[-np.ones(len(C)), np.ones(len(old))], integrality=np.ones(len(C)+len(old)),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, -np.inf, upper),
                  options={'time_limit': args.seconds, 'presolve': False, 'mip_rel_gap': 0, 'threads': 2})
    chosen = np.array([], dtype=np.int64) if result.x is None else np.flatnonzero(result.x[:len(C)] > .5)
    deleted = np.unique(np.concatenate([blockers[i] for i in chosen])) if len(chosen) else np.array([], dtype=np.int64)
    gain = len(chosen)-len(deleted)
    report = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  pool_sha256=hashlib.sha256(pool.read_bytes()).hexdigest(),
                  candidates=len(C), unique_blockers=len(old), candidate_conflicts=len(ii),
                  exact_blocker_replay='PASS', solver_status=int(result.status), message=str(result.message),
                  numerical_gain_bound=None if getattr(result, 'mip_dual_bound', None) is None else -float(result.mip_dual_bound),
                  added=len(chosen), removed=len(deleted), gain=gain,
                  scope='Numerical finite-pool search. Positive changes are checked exactly; a timeout gives no exclusion.')
    if gain > 0:
        keep = np.ones(len(A), dtype=bool); keep[deleted] = False
        FA = np.r_[A[keep], C[chosen]]; FB = np.r_[B[keep], np.zeros_like(C[chosen])]
        assert len(np.unique(np.c_[FA, FB], axis=0)) == len(FA)
        for i, point in enumerate(C[chosen]):
            da = np.zeros(len(FA), dtype=np.int64); db = np.zeros(len(FA), dtype=np.int64)
            for j in np.flatnonzero(point):
                da += point[j]*FA[:, j]; db += point[j]*FB[:, j]
            da[len(A)-len(deleted)+i] = -32
            assert not np.any(above16(da, db))
        np.savez_compressed(args.out/'proposal.npz', rational=FA, radical=FB, radicand=3,
                            squared_norm=32, inner_product_limit=16)
        np.savez_compressed(args.out/'exchange.npz', added_A=C[chosen], removed_A=A[deleted], removed_B=B[deleted])
        report['changed_geometry'] = 'PASS; full final replay still required'
    report['seconds'] = time.monotonic()-started
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
