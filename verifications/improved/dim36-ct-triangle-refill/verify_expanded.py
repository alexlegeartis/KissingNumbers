"""Expand the new 36D witness and check every pair touching a changed point."""

import json
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

from verify import pairing, read_geometry


def blocks(lanes):
    a, b, c = (len(x) for x in lanes)
    return np.column_stack((
        np.repeat(lanes[0], b * c, axis=0),
        np.tile(np.repeat(lanes[1], c, axis=0), (a, 1)),
        np.tile(lanes[2], (a * b, 1)),
    ))


def expansion(g):
    kept = []
    changed = []
    zeros = np.zeros_like(g.axial)
    for k in range(3):
        lanes = [zeros] * 3
        lanes[k] = g.axial
        kept.append(np.concatenate(lanes, axis=1))
    for i, row in enumerate(g.cells[1:], start=1):
        if np.all(row != 0):
            continue
        if i not in g.switches:
            kept.append(blocks([g.cosets[j] for j in row]))
            continue
        small, large, empty = pairing(row, g)
        for v in g.cosets[row[large]]:
            lanes = [np.zeros((1, 12), dtype=np.int64) for _ in range(3)]
            lanes[large] = v[None, :]
            if tuple(v) in g.switches[i]:
                lanes[empty] = -g.cosets[row[small]]
                changed.append(blocks(lanes))
            else:
                lanes[small] = g.cosets[row[small]]
                kept.append(blocks(lanes))
    for row in sorted(g.final):
        a, b, c, p, q, r = row
        v = np.concatenate((g.roots[a, p], g.roots[b, q], g.roots[c, r]))
        (changed if row in g.inserted else kept).append(v[None, :])
    return np.concatenate(kept), np.concatenate(changed)


def main():
    assert __debug__, "run without -O"
    assert np.finfo(np.float32).nmant == 23
    g = read_geometry()
    kept, changed = expansion(g)
    X = np.concatenate((kept, changed))
    assert X.shape == (g.points, 36)
    assert len(np.unique(X, axis=0)) == len(X)
    v = X.reshape(-1, 3, 12)
    w = (v @ g.gram).reshape(-1, 36)
    assert np.all((X * w).sum(axis=1) == 12)
    largest = 36 * int(np.max(np.abs(w))) * int(np.max(np.abs(X)))
    assert largest < 2**24
    # These binary32 matrix products are exact integer products and
    # accumulations. The bound above applies to every partial sum.
    Y = X.astype(np.float32).T.copy(order="C")
    most = -12
    checks = 0
    k = len(kept)
    for start in range(k, len(X), 64):
        stop = min(start + 64, len(X))
        products = w[start:stop].astype(np.float32) @ Y[:, :stop]
        for j, row in enumerate(products):
            i = start + j
            most = max(most, int(row[:i].max()))
            checks += i
    assert checks == len(changed) * k + len(changed) * (len(changed) - 1) // 2
    assert most <= 6
    print(json.dumps({
        "status": "PASS",
        "points": len(X),
        "changed_points": len(changed),
        "changed_pairs": checks,
        "maximum_distinct_product": most,
        "binary32_integer_partial_sum_bound": largest,
    }, indent=2))


if __name__ == "__main__":
    main()
