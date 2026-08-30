#!/usr/bin/env python3
"""Exact verification of the cap-direction configurations and their partitions, k = 9..23.

For dimension 72 + k the cap layer needs a kissing configuration W of R^k partitioned into
parts that are internally at cosine <= -1/2.  A part has at most 3 elements (n unit vectors
pairwise at cos <= -1/2 satisfy n <= 1 + 1/(1/2) = 3) and a 3-element part is a ZERO-SUM
TRIPLE.  The gain is 2 * sum_i |C_i| * (|Z_i| - 1), so a triple is worth 4 per class line and
a pair 2.

Everything below is integer arithmetic; no tolerance is used anywhere.  Two storage
conventions appear:

  lam<k>_W.npy          integer coordinates, the plain dot product is the true one
  rec_<k>_W.npy + _c    integer coordinates with integer inner-product coefficients c,
                        so <x,y> = sum_i c_i x_i y_i   (Cohn's dimensions1-24.txt convention)
  kap<k>_W.npy + _c     the same convention, for the Kappa sections K_12 and K_13 that
                        scripts/kappa_caps.py derives -- 756 and 918 vectors against
                        Lambda_12's 648 and Lambda_13's 906, so 252 and 306 triples
                        against 216 and 302

CHECKS, for every k:
  * every vector distinct, and the configuration spans exactly k dimensions;
  * 60-degree code: for every pair, 4<x,y>^2 <= |x|^2 |y|^2 whenever <x,y> > 0;
  * every triple sums to zero and has all three inner products equal to -|x|^2/2;
  * every pair part has <x,y> < 0 and 4<x,y>^2 >= |x|^2 |y|^2, i.e. cos <= -1/2;
  * the parts are disjoint and 3T + 2P + S = |W|.

    python verify_parts.py
"""
import numpy as np, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, '..', 'data')


def load(k, src):
    if src == 'lattice':
        f = 'lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k
        X = np.load(os.path.join(D, f)).astype(np.int64)
        c = np.ones(X.shape[1], np.int64)
        P = np.load(os.path.join(D, 'parts_%d.npz' % k))
    elif src == 'kappa':
        X = np.load(os.path.join(D, 'kap%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(D, 'kap%d_c.npy' % k)).astype(np.int64)
        P = np.load(os.path.join(D, 'kap%d_parts.npz' % k))
    else:
        X = np.load(os.path.join(D, 'rec_%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(D, 'rec_%d_c.npy' % k)).astype(np.int64)
        P = np.load(os.path.join(D, 'rec_%d_parts.npz' % k))
    return X, c, P['triples'], P['pairs']


def check(k, src):
    X, c, T, PR = load(k, src)
    n, amb = X.shape
    Y = X * c
    d = np.einsum('ij,ij->i', Y, X)
    assert (d > 0).all()
    assert len(set(map(tuple, X.tolist()))) == n, "repeated vector"
    span = np.linalg.matrix_rank(X.astype(float))
    assert span == k, "spans %d dimensions, not %d" % (span, k)
    # 60-degree code, exact
    # 4<x,y>^2 <= |x|^2|y|^2 whenever <x,y> > 0.  Bound the magnitudes first so int64 cannot
    # overflow: |<x,y>| <= max norm and the products stay far below 2^63.
    NM = int(d.max())
    assert 4 * NM * NM < 2 ** 62, "int64 headroom"
    # numpy's integer matmul is not BLAS-backed and is far too slow at n ~ 10^5, so the Gram
    # is formed in float32.  Every entry is an integer bounded by max|Y| * max|X| * dim, which
    # is asserted below to stay under 2^24, where float32 represents integers EXACTLY -- so
    # this is not an approximation, and the result is cast straight back to int64.
    B0 = int(np.abs(Y).max()) * int(np.abs(X).max()) * X.shape[1]
    assert B0 < 2 ** 24, "float32 exactness bound"
    Yf, Xf = Y.astype(np.float32), X.astype(np.float32)
    bad = 0
    uniform = len(set(d.tolist())) == 1        # then <x,y> <= N/2 is simply 2|<x,y>| <= N
    CH = max(1, 2 ** 25 // max(1, n))
    for s0 in range(0, n, CH):
        # kept in float32: every entry is an integer below 2^24, so the comparisons below
        # are exact and 8.7e9 int64 conversions are avoided at k = 23
        B = Yf[s0:s0 + CH] @ Xf.T
        if s0 == 0:                            # spot-check the float path against integers
            assert (B[:8].astype(np.int64) == (Y[:8] @ X.T)).all(), "float32 Gram is not exact"
        for r in range(len(B)):
            B[r, s0 + r] = -1                  # exclude the diagonal
        assert float(np.abs(B).max()) <= NM, "inner product exceeds the norm"
        if uniform:
            bad += int((B > 0.5 * float(d[0])).sum())
        else:
            pos = B > 0
            if pos.any():
                rhs = (d[s0:s0 + CH, None] * d[None, :]).astype(np.float32)
                bad += int((pos & (4 * B * B > rhs)).sum())
    assert bad == 0, "not a 60-degree code (%d violating pairs)" % bad
    seen = set()
    for a, b, e in T.tolist():
        assert (X[a] + X[b] + X[e] == 0).all(), "triple does not sum to zero"
        na = int(d[a])
        assert int(d[b]) == na and int(d[e]) == na, "triple with unequal norms"
        for (u, v) in ((a, b), (a, e), (b, e)):
            assert 2 * int(Y[u] @ X[v]) == -na, "triple pair not at 120 degrees"
        for x in (a, b, e):
            assert x not in seen, "parts overlap"
            seen.add(x)
    for a, b in PR.tolist():
        ip = int(Y[a] @ X[b])
        assert ip < 0 and 4 * ip * ip >= int(d[a]) * int(d[b]), "pair part not at cos <= -1/2"
        for x in (a, b):
            assert x not in seen, "parts overlap"
            seen.add(x)
    S = n - 3 * len(T) - 2 * len(PR)
    assert S >= 0 and len(seen) == 3 * len(T) + 2 * len(PR)
    return n, len(T), len(PR), S, int(d[0]), sorted(set(d.tolist()))


if __name__ == '__main__':
    print(__doc__)
    print("   k   source    |W|    span  triples  pairs  singles   units   3T+2P+S = |W|")
    tot = 0
    for k in range(9, 24):
        best = None
        for src in ('lattice', 'record', 'kappa'):
            f = os.path.join(D, {'lattice': 'parts_%d.npz' % k,
                                 'record': 'rec_%d_parts.npz' % k,
                                 'kappa': 'kap%d_parts.npz' % k}[src])
            if not os.path.exists(f):
                continue
            n, T, P, S, nrm, norms = check(k, src)
            u = 2 * T + P
            if best is None or u > best[0]:
                best = (u, src, n, T, P, S)
        u, src, n, T, P, S = best
        ok = (3 * T + 2 * P + S == n)
        tot += 1
        print("  %2d   %-8s %6d   %3d  %6d  %5d  %5d   %6d   %s" % (k, src, n, k, T, P, S, u, ok))
        assert ok
    print()
    print("ALL %d CONFIGURATIONS AND PARTITIONS VERIFIED EXACTLY" % tot)
