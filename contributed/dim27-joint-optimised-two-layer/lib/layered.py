# -*- coding: utf-8 -*-
"""Exact verification of a layered configuration in R^(24+k), k = 2, 3.

Coordinates: Cohn units, Leech minimal vectors integer of norm 32.  A cap point is
(y/3, (2/sqrt3) u) with y an INTEGER vector of norm 192 (a head, |x|^2 = 8/3 in norm-4
units) and u a unit direction in R^k; the axis points are (0, a) with |a| = 2.  Every
point has norm 32 (Cohn) and the kissing condition is <P,Q> <= 16 in Cohn units, i.e.

    head-equator :  y.z <= 48                                (integer)
    head-head    :  y.y'/9 + (32/3) u.u' <= 16,  i.e.  y.y' <= 144 - 96 u.u'   (integer once
                    u.u' is known exactly; here u.u' in {1, 1/2, 0, -1/2, -1})
    head-axis    :  (2/sqrt3) u.a <= 2 in norm-4 units (the cap's y-part has length 2/sqrt3,
                    the axis point length 2), decided exactly in Q(sqrt2, sqrt3) with sympy
    axis-axis    :  a.a' <= 2 (norm-4 units), exactly with sympy
    equator      :  Leech minus the owners; equator-equator is the Leech lattice.

Directions and axis points are given as sympy vectors; sides are indices into the list of
direction triangles.  The function returns the point count and raises AssertionError on
any failure; the caller prints.
"""
import numpy as np
import sympy as sp


def exact_max_cos(T1, T2):
    return max(sp.nsimplify(sp.simplify(u.dot(v))) for u in T1 for v in T2)


def verify(M, Y, side, triangles, axis, tau, dim, log):
    """M: Leech minimal vectors (int, Cohn).  Y: heads (int, Cohn, norm 192).  side: index of
    the direction triangle of each head.  triangles: list of lists of sympy unit vectors in
    R^k.  axis: list of sympy vectors of norm 2 in R^k.  tau: expected number of axis points."""
    M = np.asarray(M, dtype=np.int64); Y = np.asarray(Y, dtype=np.int64); side = np.asarray(side)
    K = len(triangles)
    assert len(axis) == tau
    assert (np.einsum('ij,ij->i', Y, Y) == 192).all(), "every head has |y|^2 = 192"
    assert side.min() >= 0 and side.max() < K
    # directions: unit, zero-sum triangles
    for T in triangles:
        assert len(T) == 3 and all(sp.simplify(u.dot(u) - 1) == 0 for u in T)
        assert all(sp.simplify(c) == 0 for c in (T[0] + T[1] + T[2]))
    # 1. head vs equator: deletions, exactly one per class head, none per free head, owners distinct
    D = (Y @ M.T) > 48
    cnt = D.sum(1)
    assert cnt.max() <= 1, "a head removes at most one minimal vector"
    owners = [int(np.nonzero(D[i])[0][0]) for i in range(len(Y)) if cnt[i] == 1]
    assert len(set(owners)) == len(owners), "owners distinct"
    nclass, nfree = len(owners), int((cnt == 0).sum())
    log("heads %d = %d class + %d free; owners distinct" % (len(Y), nclass, nfree))
    # 2. head vs head, by side pair, integer thresholds from exact direction cosines
    thr = np.zeros((K, K), dtype=np.int64)
    for a in range(K):
        for b in range(K):
            c = exact_max_cos(triangles[a], triangles[b])
            t = sp.nsimplify(144 - 96 * c)
            assert t.is_integer, "threshold must be an integer for these directions"
            thr[a, b] = int(t)
    log("side-pair thresholds (Cohn units): %s" % thr.tolist())
    G = Y @ Y.T
    np.fill_diagonal(G, -10 ** 6)
    T = thr[side[:, None], side[None, :]]
    assert (G <= T).all(), "head-head"
    tight = int((G == T).sum() // 2)
    log("all %d head pairs admissible; %d at the bound" % (len(Y) * (len(Y) - 1) // 2, tight))
    # 3. head vs axis and axis vs axis, exactly (norm-4 units: cap y-part has length 2/sqrt3)
    h = 2 / sp.sqrt(3)
    worst = max((sp.simplify(h * u.dot(a)) for T in triangles for u in T for a in axis), key=lambda e: sp.N(e, 30))
    gap = sp.simplify(worst - 2)
    assert gap == 0 or sp.N(gap, 50) < -sp.Rational(1, 10**20), "head-axis"
    log("head-axis maximum inner product %s = %.6f (bound 2)" % (worst, float(worst)))
    for i, a in enumerate(axis):
        assert sp.simplify(a.dot(a) - 4) == 0
        for b in axis[i + 1:]:
            assert sp.simplify(b.dot(a) - 2).is_nonpositive, "axis-axis"
    log("%d axis points of norm 4, pairwise inner product <= 2" % tau)
    # 4. count
    count = 196560 - nclass + 3 * len(Y) + tau
    log("count = 196560 - %d + 3 * %d + %d = %d" % (nclass, len(Y), tau, count))
    return count


def hexagon():
    tri = [[sp.Matrix([sp.cos(sp.pi * (30 + 60 * j) / 180), sp.sin(sp.pi * (30 + 60 * j) / 180)]) for j in js]
           for js in ((0, 2, 4), (1, 3, 5))]
    axis = [sp.Matrix([2 * sp.cos(sp.pi * k / 3), 2 * sp.sin(sp.pi * k / 3)]) for k in range(6)]
    return tri, axis, 6


def cuboctahedron():
    r2 = sp.sqrt(2)
    tris = [[(1, 1, 0), (-1, 0, -1), (0, -1, 1)], [(1, -1, 0), (-1, 0, 1), (0, 1, -1)],
            [(-1, -1, 0), (1, 0, -1), (0, 1, 1)], [(-1, 1, 0), (1, 0, 1), (0, -1, -1)]]
    tri = [[sp.Matrix(v) / r2 for v in T] for T in tris]          # the 12 vertices, in four zero-sum triangles
    assert len({tuple(v) for T in tris for v in T}) == 12
    # the cuboctahedron rotated by 45 degrees about the third axis, scaled to norm 2
    axis = [sp.Matrix(a) for a in [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0)]] + \
           [sp.Matrix([s1, s2, s3 * r2]) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
    return tri, axis, 12
