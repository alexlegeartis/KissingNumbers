"""Search exact BW14 normal-module transplants into the current 30D source.

The 14-coordinate BW16 section has 1260 points. Its E8 equator and D6 tail
are shared with the incumbent, leaving 960 heads versus 932 incumbent heads.
Orientations use E8 Weyl reflections and signed permutations. Every positive
exchange is checked directly against the complete input, not just projections.
The validity of unchanged source pairs is an explicit inherited premise.
"""
import argparse
import hashlib
import json
import math
import time
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

NORMAL = np.array([0, 2, 4, 5, 6, 10, 11, 23])


def bw14():
    rows = []
    for i, j in combinations(range(14), 2):
        for a, b in product((-4, 4), repeat=2):
            x = np.zeros(14, np.int16)
            x[i], x[j] = a, b
            rows.append(x)
    supports = []
    for a in range(1, 16):
        for value in (0, 1):
            s = [i for i in range(16) if (a & i).bit_count() % 2 == value]
            if 14 in s or 15 in s:
                continue
            supports.append(s)
            for signs in product((-2, 2), repeat=8):
                if sum(t < 0 for t in signs) % 2:
                    continue
                x = np.zeros(14, np.int16)
                x[s] = signs
                rows.append(x)
    X = np.array(rows, np.int16)
    assert len(supports) == 7 and len(X) == 1260
    assert np.all(np.sum(X.astype(np.int64) ** 2, axis=1) == 32)
    D = X.astype(np.int32) @ X.astype(np.int32).T
    np.fill_diagonal(D, -32)
    assert D.max() == 16
    q = np.sum(X[:, 8:].astype(np.int64) ** 2, axis=1)
    assert dict(zip(*np.unique(q, return_counts=True))) == {0: 240, 16: 960, 32: 60}
    return X[q == 16]


def mutate(A, kind, choice):
    Z = A.copy()
    if kind == "root":
        dots = Z[:, :8] @ choice
        assert not np.any(dots % 4)
        Z[:, :8] -= (dots // 4)[:, None] * choice
    elif kind == "core_swap":
        i, j = choice
        Z[:, [i, j]] = Z[:, [j, i]]
    elif kind == "tail_swap":
        i, j = choice
        Z[:, [8 + i, 8 + j]] = Z[:, [8 + j, 8 + i]]
    else:
        Z[:, 8 + choice] *= -1
    return Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--iterations", type=int, default=200000)
    ap.add_argument("--seconds", type=float, default=250)
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--temperature", type=float, default=3)
    ap.add_argument("--restart", type=int, default=2500)
    ap.add_argument("--initial", type=Path)
    ap.add_argument("--normal", default=",".join(map(str, NORMAL)))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source = np.load(args.source, allow_pickle=False)
    X = source["points"].astype(np.int16)
    if int(source["squared_norm"]) == 128:
        assert not np.any(X % 2)
        X //= 2
    assert X.ndim == 2 and X.shape[1] == 30
    assert np.all(np.sum(X.astype(np.int64) ** 2, axis=1) == 32)
    normal = np.array([int(t) for t in args.normal.split(",")])
    assert len(normal) == len(set(normal)) == 8 and np.all((normal >= 0) & (normal < 24))
    cols = np.r_[normal, np.arange(24, 30)]
    other = [j for j in range(24) if j not in normal]
    projection = X[:, cols]
    q = np.sum(X[:, 24:].astype(np.int64) ** 2, axis=1)
    h = np.sum(X[:, normal].astype(np.int64) ** 2, axis=1)
    module = ~np.any(X[:, other], axis=1)
    relevant = ((q == 8) & (h > 8)) | ((q == 16) & (h > 0))
    ids = np.flatnonzero(relevant)
    external = relevant & ~module
    assert np.sum(module & (q == 0)) == 240
    assert np.sum(module & (q == 8)) == 0
    old_heads = int(np.sum(module & (q == 16)))
    nominal_gain = 960 - old_heads
    A0 = projection[external]
    base = bw14() if args.initial is None else np.load(args.initial)["heads"]
    assert base.shape == (960, 14)
    assert len(np.unique(base, axis=0)) == 960
    assert np.all(np.sum(base[:, :8].astype(np.int64)**2, axis=1) == 16)
    assert np.all(np.sum(base[:, 8:].astype(np.int64)**2, axis=1) == 16)
    base_products = base.astype(np.int32) @ base.astype(np.int32).T
    np.fill_diagonal(base_products, -32)
    assert base_products.max() <= 16
    basef = base.astype(np.float32)
    A = A0.copy()
    R = 4 * np.eye(14, dtype=np.int16)
    oldset = {row.tobytes() for row in projection}
    records = []
    rng = np.random.default_rng(args.seed)
    best_union = 0

    def cost(a):
        # Integer dot products remain exact in binary32 for these small entries.
        return int(np.any(basef @ a.astype(np.float32).T > 16, axis=0).sum())

    score = cost(A)
    best = score
    bestA, bestR = A.copy(), R.copy()

    def record(iteration):
        nonlocal best_union
        num = base.astype(np.int32) @ R.astype(np.int32).T
        assert not np.any(num % 4)
        Y = (num // 4).astype(np.int16)
        assert np.array_equal(A.astype(np.int32) * 4, A0.astype(np.int32) @ R)
        assert np.array_equal(R.astype(np.int32) @ R.astype(np.int32).T, 16 * np.eye(14, dtype=np.int32))
        new = np.array([row for row in Y if row.tobytes() not in oldset], np.int16).reshape((-1, 14))
        cross = new.astype(np.float32) @ projection[ids].astype(np.float32).T > 16
        graph = csr_matrix(cross)
        matching = maximum_bipartite_matching(graph, perm_type="column")
        gain = len(new) - int(np.sum(matching >= 0))
        rec = dict(iteration=iteration, external_deletions=score, whole_module_gain=nominal_gain-score,
                   union_gain=gain, new_points=len(new), cross_conflicts=graph.nnz)
        np.savez_compressed(args.out / "best_module.npz", heads=Y,
                            transformation_numerator=R, transformation_denominator=4,
                            base_heads=base)
        if gain > 0:
            reverse = np.full(len(ids), -1, np.int64)
            for i, j in enumerate(matching):
                if j >= 0:
                    reverse[j] = i
            left = set(map(int, np.flatnonzero(matching < 0)))
            right, queue = set(), list(left)
            for i in queue:
                for j0 in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
                    j = int(j0)
                    if matching[i] == j or j in right:
                        continue
                    right.add(j)
                    if reverse[j] >= 0 and int(reverse[j]) not in left:
                        left.add(int(reverse[j]))
                        queue.append(int(reverse[j]))
            removed = ids[sorted(right)]
            p = new[sorted(left)]
            added = np.zeros((len(p), 30), np.int16)
            added[:, normal], added[:, 24:] = p[:, :8], p[:, 8:]
            keep = np.ones(len(X), bool)
            keep[removed] = False
            F = np.r_[X[keep], added]
            assert len(F) == len(X) + gain and len(np.unique(F, axis=0)) == len(F)
            for i, row in enumerate(added):
                dots = np.zeros(len(F), np.int32)
                for j in np.flatnonzero(row):
                    dots += F[:, j].astype(np.int32) * int(row[j])
                dots[len(X)-len(removed)+i] = -32
                assert dots.max() <= 16
            np.savez_compressed(args.out / f"exchange30_gain{gain}.npz", removed=X[removed], added=added)
            np.savez_compressed(args.out / f"proposal30_{len(F)}.npz", points=F, squared_norm=32, inner_product_limit=16)
            rec.update(exact_changed_pairs="PASS", added=len(added), removed=len(removed))
            best_union = max(best_union, gain)
        records.append(rec)
        print(json.dumps(rec), flush=True)

    print(json.dumps(dict(old_module=int(module.sum()), external_points=len(A0), relevant_points=len(ids))), flush=True)
    record(0)
    completed = 0
    for it in range(1, args.iterations + 1):
        if time.monotonic() - started > args.seconds:
            break
        u = rng.random()
        if u < .5:
            kind, choice = "root", rng.choice([-1, 1], 8)
            choice[-1] = np.prod(choice[:-1])
        elif u < .65:
            kind, choice = "core_swap", rng.choice(8, 2, replace=False)
        elif u < .85:
            kind, choice = "tail_swap", rng.choice(6, 2, replace=False)
        else:
            kind, choice = "tail_sign", int(rng.integers(6))
        candidate = mutate(A, kind, choice)
        value = cost(candidate)
        temperature = args.temperature * (1 - (it % args.restart) / args.restart) + .1
        if value <= score or rng.random() < math.exp((score-value) / temperature):
            A, R, score = candidate, mutate(R, kind, choice), value
            if score < best:
                best, bestA, bestR = score, A.copy(), R.copy()
                record(it)
        if it % args.restart == 0:
            A, R, score = bestA.copy(), bestR.copy(), best
        completed = it
        if best == 0:
            break  # The full-module objective cannot have fewer than zero victims.
    report = dict(source=str(args.source), source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                  source_count=len(X), old_module=int(module.sum()), external_potential_blockers=len(A0),
                  module_head_count=len(base), best_external_cost=best, nominal_full_module_gain=nominal_gain,
                  normal_octad=normal.tolist(),
                  best_union_gain=best_union, iterations=completed, seconds=time.monotonic()-started,
                  seed=args.seed, records=records,
                  scope="Bounded E8-Weyl/signed-tail BW14 orientation search; not an exhaustive exclusion.")
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
