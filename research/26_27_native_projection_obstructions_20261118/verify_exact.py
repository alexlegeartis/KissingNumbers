"""Exhaustive certificate replay for an integral 60-degree spherical code."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np


def verify(points, block=256):
    if points.ndim != 2 or points.dtype.kind != "i" or not points.size:
        raise ValueError("Expected a nonempty signed-integer coordinate matrix.")
    n, dimension = points.shape
    if int(points.min()) < -128 or int(points.max()) > 128:
        raise ValueError("Coordinate range exceeds this verifier.")
    a = points.astype(np.int32)
    if not np.all(np.einsum("ij,ij->i", a, a) == 128):
        raise ValueError("Squared norm differs from 128.")
    signs = np.sign(a[np.arange(n), np.argmax(a != 0, axis=1)])
    orient = (a * signs[:, None]).astype(np.int16)
    _, idx, inv, counts = np.unique(
        orient, axis=0, return_index=True, return_inverse=True, return_counts=True
    )
    sums = np.bincount(inv, weights=signs).astype(int)
    if np.any(counts > 2) or np.any((counts == 2) & (sums != 0)):
        raise ValueError("Duplicate point.")
    paired = orient[idx[counts == 2]].astype(np.float32)
    single = a[idx[counts == 1]].astype(np.float32)
    budget = dimension * max(abs(int(a.min())), abs(int(a.max()))) ** 2
    if budget >= 2**24:
        raise ValueError("Single-precision integral exactness bound exceeded.")
    maximum = -128
    for left, right, tri, absolute in [
        (paired, paired, True, True),
        (paired, single, False, True),
        (single, single, True, False),
    ]:
        for start in range(0, len(left), block):
            first = start if tri else 0
            if first >= len(right):
                continue
            gram = left[start : start + block] @ right[first:].T
            if absolute:
                np.abs(gram, out=gram)
            if tri:
                np.fill_diagonal(gram, -128)
            mx = int(gram.max())
            if mx > 64:
                raise ValueError(f"Inner product {mx} exceeds 64.")
            maximum = max(maximum, mx)
    p, u = len(paired), len(single)
    pairs = 4 * math.comb(p, 2) + p + 2 * p * u + math.comb(u, 2)
    if pairs != math.comb(n, 2):
        raise ValueError("Pair coverage mismatch.")
    return dict(
        dimension=dimension,
        count=n,
        squared_norm=128,
        maximum_dot=maximum,
        antipodal_pairs=p,
        singletons=u,
        unordered_pairs=pairs,
        absolute_sum_bound=budget,
        coordinate_sha256=hashlib.sha256(a.astype("<i2").tobytes()).hexdigest(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("npz", type=Path)
    args = parser.parse_args()
    with np.load(args.npz, allow_pickle=False) as f:
        print(json.dumps(verify(f["points"]), indent=2))
