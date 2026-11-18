"""Check the absence of one-point extensions at each native BW16 projection."""

import argparse
from collections import Counter
from itertools import product
import hashlib
from fractions import Fraction
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import qr

I = [0, 2, 4, 5, 6, 10, 11, 23]
J = [j for j in range(24) if j not in I]


def full_rank_modulo(A):
    a = np.asarray(A, dtype=np.int64).copy() % 257
    for i in range(len(a)):
        idx = np.flatnonzero(a[i:, i])
        if not len(idx):
            return False
        j = i + idx[0]
        a[[i, j]] = a[[j, i]]
        a[i] = (a[i] * pow(int(a[i, i]), -1, 257)) % 257
        a[i + 1 :] = (a[i + 1 :] - a[i + 1 :, i, None] * a[i]) % 257
    return True


def nullbasis(E):
    # Approximate factorizations suggest a minor and rational columns.
    # Exact finite-field rank and integer multiplication certify the result.
    _, R, ix = qr(E.T.astype(float), pivoting=True, mode="economic")
    r = int(np.count_nonzero(np.abs(np.diag(R)) > 1e-7))
    A = E[ix[:r]]
    _, _, iy = qr(A.astype(float), pivoting=True, mode="economic")
    C, F = iy[:r], iy[r:]
    if not full_rank_modulo(A[:, C]):
        raise ValueError("Rank witness failed.")
    B = np.zeros((E.shape[1], len(F)), dtype=object)
    B[F] = np.eye(len(F), dtype=int)
    v = np.linalg.solve(A[:, C].astype(float), -A[:, F].astype(float))
    for i, row in enumerate(C):
        for j in range(len(F)):
            B[row, j] = Fraction(float(v[i, j])).limit_denominator(128)
    den = math.lcm(*(Fraction(b).denominator for b in B.ravel()))
    N = (B * den).astype(np.int64)
    if den > 10000 or np.max(np.abs(N), initial=0) > 10000 or np.any(E @ N):
        raise ValueError("Nullspace witness failed.")
    return N, den


def poly_max(G, den, rows, rhs):
    # z spans a box because each free physical coordinate equals den*z_i.
    F = Fraction
    bound = F(8, den)
    if len(G) == 1:
        p = rows[:, 0]
        hi = bound
        lo = -bound
        for sign in [1, -1]:
            ids = np.flatnonzero(sign * p > 0)
            if not len(ids):
                continue
            j = ids[np.argmin(rhs[ids] / (sign * p[ids]))]
            c = F(int(rhs[j]), sign * int(p[j]))
            if sign == 1:
                hi = min(hi, c)
            else:
                lo = max(lo, -c)
        if lo > hi:
            return F(0)
        return max(lo * lo, hi * hi) * int(G[0, 0])
    if len(G) != 2:
        raise ValueError("Unexpected free dimension.")
    poly = [(-bound, -bound), (bound, -bound), (bound, bound), (-bound, bound)]
    lower = {}
    for (a, b), r in zip(rows.tolist(), rhs.tolist()):
        lower[a, b] = min(r, lower.get((a, b), r))
    for (a, b), r in lower.items():
        val = [a * x + b * y - r for x, y in poly]
        new = []
        for i, P in enumerate(poly):
            j = (i + 1) % len(poly)
            Q = poly[j]
            if val[i] <= 0:
                new.append(P)
            if (val[i] < 0 < val[j]) or (val[j] < 0 < val[i]):
                t = val[i] / (val[i] - val[j])
                new.append(tuple(u + t * (v - u) for u, v in zip(P, Q)))
        poly = new
        if not poly:
            return F(0)
    return max(
        x * x * int(G[0, 0]) + 2 * x * y * int(G[0, 1]) + y * y * int(G[1, 1])
        for x, y in poly
    )


def exclude(N, den, H, rhs):
    # Work in cube coordinates u with y=8*N*u/den. Every candidate of norm 8
    # occurs inside [-1,1]^q, since N has a den*identity minor.
    q = N.shape[1]
    G = N.T @ N
    C = H.astype(np.int64) @ N
    C, iv = np.unique(C, axis=0, return_inverse=True)
    b = np.full(len(C), 1000, np.int64)
    np.minimum.at(b, iv, rhs)
    keep = 8 * np.abs(C).sum(axis=1) > b * den
    C = 8 * C[keep]
    b = b[keep] * den
    radius = np.abs(C).sum(1)
    if q > 7:
        return False, 0
    corners = np.array(list(product([-1, 1], repeat=q)), np.int64)
    Q = np.zeros((1, q), np.int64)
    for dep in range(11):
        L = 2**dep
        # Exact bound for all products with the current cube coordinates.
        if (L + 1) * np.abs(C).sum(axis=1).max(initial=0) >= 2**60 or (
            L + 1
        ) ** 2 * abs(G).sum() >= 2**60:
            raise ValueError("Range")
        live = []
        for at in range(0, len(Q), 128):
            M = Q[at : at + 128]
            M = M[(M @ C.T <= L * b + radius).all(1)]
            if not len(M):
                continue
            MG = M @ G
            possible = (M * MG).sum(1) + 2 * np.abs(MG).sum(1) + np.abs(
                G
            ).sum() >= L**2 * den**2
            live.extend(M[possible])
        if not live:
            return True, dep
        if len(live) > 18000:
            return False, dep
        Q = (2 * np.array(live)[:, None, :] + corners).reshape(-1, q)
    return False, 10


def native_projections():
    # Project the standard Leech shell to the complement of the named Golay octad.
    # Only the (64,64) layer is needed. All sign patterns on its tetrads occur.
    g = sum(1 << i for i in (0, 2, 4, 5, 6, 10, 11))
    words = []
    for m in range(4096):
        word = 0
        for i in range(12):
            if m & (1 << i):
                word ^= g << i
        words.append(word | ((word.bit_count() % 2) << 23))
    tetrads = set()
    for word in words:
        if word.bit_count() == 8:
            support = tuple(j for j, k in enumerate(J) if word & (1 << k))
            if len(support) == 4:
                tetrads.add(support)
    signs = {tuple(2 if word & (1 << k) else -2 for k in J) for word in words}
    axes = 8 * np.eye(16, dtype=np.int16)
    pieces = [axes, -axes, np.array(sorted(signs), dtype=np.int16)]
    for block in sorted(tetrads):
        values = np.zeros((16, 16), dtype=np.int16)
        values[:, block] = np.array(list(product([-4, 4], repeat=4)), dtype=np.int16)
        pieces.append(values)
    W = np.unique(np.vstack(pieces), axis=0)
    if (
        len(W) != 4320
        or len(tetrads) != 140
        or len(signs) != 2048
        or np.any((W * W).sum(1) != 64)
    ):
        raise ValueError("Native projection enumeration failed.")
    return W


def certify(path):
    with np.load(path, allow_pickle=False) as z:
        raw = z["points"] if "points" in z else z["integer_points"]
        rational_rows = (
            len(z["fraction_denominators"]) if "fraction_denominators" in z else 0
        )
    if (
        raw.ndim != 2
        or raw.dtype.kind != "i"
        or raw.size == 0
        or raw.min() < -8
        or raw.max() > 8
    ):
        raise ValueError("Unexpected integer coordinate matrix.")
    X = raw.astype(np.int16)
    if np.any((X * X).sum(1) != 128):
        raise ValueError("Incorrect norm.")
    d = X.shape[1]
    if not 26 <= d <= 31:
        raise ValueError("Unsupported dimension.")
    Y = X[:, J]
    A = X[:, I + list(range(24, d))]
    W = native_projections()
    H, inv = np.unique(A, axis=0, return_inverse=True)
    order = np.argsort(inv)
    cuts = np.r_[0, np.flatnonzero(np.diff(inv[order])) + 1]
    right = np.ascontiguousarray(Y[order], dtype=np.float32)
    id_by_row = {v.tobytes(): i for i, v in enumerate(H)}
    opposite = np.array([id_by_row.get((-v).tobytes(), len(H)) for v in H])
    populations = Counter(v.tobytes() for v in Y[(Y * Y).sum(1) == 64])
    hist = Counter()
    cache = {}
    low_dim_max = Fraction(0)
    maxdepth = 0
    for start in range(0, len(W), 32):
        group = W[start : start + 32]
        # Every absolute partial sum is <= 16*8*8 < 2**24.
        dots = group.astype(np.float32) @ right.T
        costs = 64 - np.maximum.reduceat(dots, cuts, axis=1).astype(np.int16)
        for b, w in zip(costs, group):
            if populations[w.tobytes()] == 2 * (d - 16):
                hist["orthoplex_capacity"] += 1
                continue
            opposing = np.r_[b, 129][opposite]
            if np.any(b + opposing < 0):
                hist["incompatible_halfspaces"] += 1
                continue
            equality = (b == 0) & (opposing == 0)
            key = np.packbits(equality).tobytes()
            if key not in cache:
                cache[key] = nullbasis(H[equality])
            N, den = cache[key]
            dim = N.shape[1]
            if dim == 0 or np.any(((H @ N) == 0).all(1) & (b < 0)):
                hist["inconsistent_with_norm"] += 1
                continue
            if dim <= 2:
                bound = poly_max(N.T @ N, den, H @ N, b)
                if bound >= 64:
                    raise ValueError(
                        f"Inconclusive exact polytope test for {w.tolist()}."
                    )
                low_dim_max = max(low_dim_max, bound)
                hist[f"exact_polytope_dim_{dim}"] += 1
            else:
                passed, depth = exclude(N, den, H, b)
                if not passed:
                    raise ValueError(
                        f"Inconclusive rational subdivision for {w.tolist()}."
                    )
                maxdepth = max(depth, maxdepth)
                hist["dyadic_subdivision"] += 1
    if sum(hist.values()) != 4320:
        raise ValueError("Incomplete projection coverage.")
    return dict(
        dimension=d,
        configuration_size=len(X) + rational_rows,
        coordinate_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        native_projections_excluded=4320,
        proof_cases=dict(hist),
        maximum_subdivision_depth=maxdepth,
        required_squared_norm=64,
        max_squared_norm_bound_among_dim_1_and_2_cases=str(low_dim_max),
        omitted_rational_constraints=rational_rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", type=Path)
    print(json.dumps(certify(parser.parse_args().coordinates), indent=2))
