"""Exact bounds for large integer products using one short integer product.

Write A = 2**s H + R, 0 <= R < 2**s. Choose s so every partial
sum in H @ B is an integer below 2**24, then compute that product in
binary32. The remainder contributes at most 2**s * max_column_l1(B).
This certifies most entries immediately. Every remaining entry is checked
with Python integers. There is no numerical tolerance or sampling.
"""
import numpy as np

_last_direction_proof = None


def remember_direction_proof(W, c, bound):
    """Share a freshly established angle bound within one verifier process."""
    global _last_direction_proof
    _last_direction_proof = W.copy(), c.copy(), int(bound)


def direction_proof(W, c):
    """Use full array equality, so reusing the proof requires identical inputs."""
    if _last_direction_proof is not None:
        X, metric, bound = _last_direction_proof
        if np.array_equal(W, X) and np.array_equal(c, metric):
            return bound
    return None


def _maximum_abs(A):
    return max(abs(int(A.min())), abs(int(A.max()))) if A.size else 0


def integer_matmul(A, B):
    """Return an int64 product after proving all intermediate sums fit."""
    assert A.dtype == B.dtype == np.dtype('int64')
    bound = _maximum_abs(A) * _maximum_abs(B) * A.shape[1]
    if bound < 2**63:
        return A @ B
    result = A.astype(object) @ B.astype(object)
    assert int(result.min()) >= -(2**63) and int(result.max()) < 2**63
    return result.astype(np.int64)


def abs_product_bound(A, B, limits, block_entries=8_000_000):
    """Prove |(A @ B)[i,j]| <= limits[i]; return violating pairs and fallbacks."""
    assert A.dtype == B.dtype == np.dtype('int64')
    assert A.ndim == B.ndim == 2 and A.shape[1] == B.shape[0]
    limits = np.asarray(limits, dtype=np.int64)
    if limits.ndim == 0:
        limits = np.full(len(A), int(limits), dtype=np.int64)
    assert limits.shape == (len(A),) and np.all(limits >= 0)
    assert _maximum_abs(B) * B.shape[0] < 2**24
    l1 = int(np.abs(B).sum(axis=0).max())
    shift = 0
    while True:
        high = A >> shift
        maximum = _maximum_abs(high)
        if maximum < 2**24 and maximum * l1 < 2**24:
            break
        shift += 1
        assert shift <= 62
    scale = 2**shift
    remainder = A - (high << shift)
    assert np.all(remainder >= 0) and _maximum_abs(remainder) < scale
    Hf, Bf = high.astype(np.float32), B.astype(np.float32)
    size = max(1, block_entries // max(1, B.shape[1]))
    bad = fallbacks = 0
    for start in range(0, len(A), size):
        product = Hf[start:start + size] @ Bf
        # Clamp only to the proved range of this integer product, so the
        # comparison threshold also has an exact binary32 representation.
        cutoff = limits[start:start + len(product)] // scale - l1
        cutoff = np.clip(cutoff, -1, 2**24).astype(np.float32)
        rows, cols = np.nonzero(np.abs(product) > cutoff[:, None])
        fallbacks += len(rows)
        for first in range(0, len(rows), 4096):
            rr, cc = rows[first:first + 4096], cols[first:first + 4096]
            low = (remainder[start + rr].astype(object) * B.T[cc].astype(object)).sum(axis=1)
            exact = product[rr, cc].astype(np.int64).astype(object) * scale + low
            bad += int((np.abs(exact) > limits[start + rr]).sum())
    return bad, fallbacks
