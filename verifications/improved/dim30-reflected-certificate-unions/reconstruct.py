"""Reproduce a 248416-point 30D kissing configuration from the pinned source.

Usage: python reconstruct.py SOURCE OUTPUT
Requires NumPy and SciPy. Use a local copy of the public input identified below.
"""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

OCTAD = [0, 2, 4, 5, 6, 10, 11, 23]
COMPLEMENT = [j for j in range(24) if j not in OCTAD]
SOURCE_SHA256 = "ece0feb7b87fcfa835338f56688785b7386352893d34bb81c98cc453b62d93c1"
TARGET_COORDINATE_SHA256 = (
    "5c9d9b8978bc2aa0b1d20d34edb30d1ff13ac5cecb4454c106aee29949be92c3"
)


def keys(x):
    x = np.ascontiguousarray(x, dtype="<i2")
    return x.view(np.dtype((np.void, x.shape[1] * 2))).ravel()


def splice(old, donor):
    _, i, j = np.intersect1d(
        keys(old), keys(donor), assume_unique=True, return_indices=True
    )
    a = np.ones(len(old), dtype=bool)
    b = np.ones(len(donor), dtype=bool)
    a[i] = False
    b[j] = False
    ia, ib = np.flatnonzero(a), np.flatnonzero(b)
    # Here the symmetric differences are small. int16 accumulates these
    # products exactly because 30 times the squared maximum entry is 1920.
    graph = csr_matrix((donor[ib] @ old[ia].T) > 64)
    match = maximum_bipartite_matching(graph, perm_type="column")
    chosen = match < 0
    reverse = np.full(len(ia), -1)
    reverse[match[~chosen]] = np.flatnonzero(~chosen)
    removed = np.zeros(len(ia), dtype=bool)
    frontier = np.flatnonzero(chosen)
    while len(frontier):
        seen = np.unique(graph[frontier].indices)
        seen = seen[~removed[seen]]
        removed[seen] = True
        next_vertices = reverse[seen]
        if np.any(next_vertices < 0):
            raise ValueError("The matching admits an augmenting path.")
        frontier = next_vertices[~chosen[next_vertices]]
        chosen[frontier] = True
    keep = np.ones(len(old), dtype=bool)
    keep[ia[removed]] = False
    return np.vstack((old[keep], donor[ib[chosen]])), [
        int(chosen.sum()),
        int(removed.sum()),
    ]


def construct(source_path, target_path):
    if hashlib.sha256(source_path.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("Wrong source archive.")
    with np.load(source_path, allow_pickle=False) as src:
        raw = src["points"]
    if raw.shape != (248351, 30) or raw.dtype.kind != "i" or np.max(np.abs(raw)) > 8:
        raise ValueError("Source array is outside the supported range.")
    src = raw.astype(np.int16)
    current = src.copy()
    trades = []
    for columns in [COMPLEMENT, OCTAD, list(range(24, 30)), list(range(30))]:
        donor = src.copy()
        donor[:, columns] *= -1
        current, trade = splice(current, donor)
        trades.append(trade)
    extra = np.zeros((4, 30), dtype=np.int16)
    head = np.array([-2, 2, 2, -6, -2, 2, -2, -2], dtype=np.int16)
    for row, (a, b) in enumerate([(1, 1), (1, -1), (-1, 1), (-1, -1)]):
        extra[row, OCTAD] = a * head
        extra[row, 28] = 8 * b
    result = np.vstack((current, extra))
    result = result[np.lexsort(result.T[::-1])]
    digest = hashlib.sha256(result.astype("<i2").tobytes()).hexdigest()
    if len(result) != 248416 or digest != TARGET_COORDINATE_SHA256:
        raise ValueError("Final cardinality or coordinate hash mismatch.")
    np.savez_compressed(
        target_path,
        points=result,
        squared_norm=128,
        inner_product_limit=64,
        coordinate_denominator=2,
    )
    print(
        json.dumps(
            dict(trades=trades, final_insertions=4, count=len(result), sha256=digest)
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    p = parser.parse_args()
    construct(p.source, p.output)
