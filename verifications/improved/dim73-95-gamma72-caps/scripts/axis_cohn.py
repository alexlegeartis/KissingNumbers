#!/usr/bin/env python3
"""The axis layers of dimensions 87, 89 and 90, from the record configuration of R^k.

    python axis_cohn.py                  re-derive and report against the shipped layers
    python axis_cohn.py --write          re-derive and write data/poles_rec_<k>.npz
    python axis_cohn.py --write 17       just that k

WHY THIS FILE EXISTS.  scripts/axis_rotate.py makes the poles a rotated copy of the CAP
DIRECTIONS, so the layer can never exceed |Z|; scripts/axis_record.py drops that at k = 19,
20 and 21, where the record configuration of R^k can be rebuilt from the Golay code in the
frame the cap file is already written in.  At k = 15, 17 and 18 it cannot: the cap file is
Lambda_k written in more coordinates than it spans, and the record configurations are not of
that shape at all.  They are, however, PUBLISHED as exact integer coordinates -- Cohn's
dimensions1-24.txt gives all three with an integer metric diag(c) -- so the source is on
disk and the only question is whether the axis space can hold it.

    k     |Z| = Lambda_k     shipped layer     tau(k)      metric of the record
    15         2340              2340           2564       diag(1^14, 2)
    17         5346              5338           5730       diag(1^16, 2)
    18         7398              7374           7654       diag(1^16, 2, 6)

THE METRIC IS WHAT MAKES IT FIT.  A configuration is carried into the axis space by a
rational similarity M with M^T M = mu diag(c), and whether one exists is decided by the
discriminant: mu^k det(diag c) has to agree with the discriminant of the axis space modulo
rational squares.  In the STANDARD metric that fails at k = 18 -- the axis space has
discriminant 3, the rank is even, so mu^18 is a square and no frame exists, which is what
verify_record.py says of dimension 90.  The record configuration of R^18 is not written in
the standard metric.  Its own is diag(1^16, 2, 6), of discriminant 12 = 3 modulo squares,
and that is exactly the discriminant of the axis space.  The obstruction was never to the
CONFIGURATION; it was to the standard metric, and this configuration does not use it.

THE FRAME, IN CLOSED FORM.  Take v_1 ... v_r mutually orthogonal cap directions, all of norm
m, as many as there are (fourteen of fifteen at k = 15, seventeen of eighteen at k = 18 --
the discriminant again).  Then

    k = 15    v_1 ... v_14, and the complement line scaled to norm 16       mu = 8
    k = 17    v_1 +- v_2, ..., v_15 +- v_16 of norm 64, and 2 v_17          mu = 64
    k = 18    those sixteen, 2 v_17, and the complement line at norm 384    mu = 64

The complement line -- the orthogonal complement of the chosen caps inside the span -- is one
dimensional, so its norm class is forced, not chosen.  It comes out at 1 modulo squares at
k = 15 where the frame needs 1, and at 6 where the frame needs 6.  That is the discriminant
being right about something it could have been wrong about.

THE ROTATION IS A PRODUCT OF RATIONAL PLANE ROTATIONS.  A full-rank Cayley transform on a
rank-17 space has denominator about d^17, so a usable D allows d = 4 at best -- twelve per
cent per entry, against the one per cent the margins here need.  One plane rotation is
different.  In the plane of two frame axes with metric weights g_i and g_j,

    R = [[a, -(g_j/g_i) b], [b, a]] / D        with   a^2 + (g_j/g_i) b^2 = D^2

is rational, is an isometry of diag(c), and has denominator D alone.  Its angles are the
rational points of a conic, which are dense in the circle.  A greedy pass over planes and
angles, each round scored on a subsample and then counted in full, reaches the numbers above
in a handful of rotations: at k = 15 and 17 the FIRST one, with D = 13, is worth more than
the entire shipped layer's shortfall.

WHAT IS HEURISTIC AND WHAT IS NOT.  The frame is a search, the plane rotations are a search,
and the source is read from a published data set.  Nothing rests on the first two:
scripts/verify_record.py re-derives the frame's properties, checks the rotation is an
isometry of the source's own metric, and tests every kept pole against every cap in exact
integer arithmetic.  The source is shipped as data/pole_src_<k>.npy so that the package
verifies without the 24 MB data set, and is checked against the data set itself whenever it
can be found.
"""
import os
import sys
import time
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from axis_rotate import directions, rows                              # noqa: E402
from axis_block import frame as ortho_caps, span_basis, preserves_span, is_isometry  # noqa: E402

PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
sys.path.insert(0, os.path.join(PKG, '..', '..', '..', 'common'))
from published import TAU                                            # noqa: E402

KS = (15, 17, 18, 23)

# The frame's shape and metric per k, and the class the complement line has to carry.  All
# three are consequences of the discriminant and are ASSERTED below, not assumed.
#   'caps'   the frame vectors are cap directions as they stand
#   'pairs'  they are v_1 +- v_2, ..., v_15 +- v_16, 2 v_17, of norms 64 and 128
#   'block'  they are A V for an integer A with A A^T = diag(c), which is what a metric
#            with entries other than 1 and 2 needs
# k = 13 IS NO LONGER HERE.  Dimension 85 used to take the record configuration of R^13 as
# its pole source, mapped into the span of Lambda_13 by a frame, and got 1106 of tau(13) =
# 1154.  On 2026-08-30 its CAP set changed from Lambda_13 to the Kappa section K_13
# (scripts/kappa_caps.py), which is larger -- 918 against 906 -- and splits perfectly into
# zero-sum triples: 306 against 302.  That is worth 32 350 spheres, against the 224 poles the
# change costs, because a triple is worth four per class line and a pole is worth one.
#
# The frame does not survive the change and could not be made to.  A frame is a rational
# similarity, so it needs the source's quadratic space to be Q-similar to the axis space; the
# record configuration of R^13 is written in the standard metric, K_13's span carries
# diag(3,1,3,1,3,1,3,1,3,1,3,1,1), and their Hasse invariants differ at 2 and at 3.  Rank 13
# is odd, so the discriminant can always be matched by choosing mu -- and the Hasse invariant
# of a rank-13 form is unchanged by scaling, so no mu fixes it.  Dimension 85's poles are a
# rotated copy of K_13 instead (scripts/kappa_poles.py).
PLAN = {
        15: dict(mu=8, shape='caps', tail=1),
        17: dict(mu=64, shape='pairs', tail=None),
        18: dict(mu=64, shape='pairs', tail=6),
        23: dict(mu=32, shape='block', tail=None)}


# --------------------------------------------------------------------------- the source
def data_set(argv=()):
    """locate Cohn's dimensions1-24.txt, or None -- the same search extract.py uses"""
    for p in list(argv) + [os.environ.get('D124', '')]:
        if p and os.path.exists(p):
            return p
    up = HERE
    for _ in range(8):
        c = os.path.join(up, 'dimensions1-24.txt')
        if os.path.exists(c):
            return c
        up = os.path.dirname(up)
    return None


def record(dim, path):
    """(X, c, dropped): the integral points of one dimension of the data set

    Two of the blocks used here are not wholly integral -- dimension 13 writes 48 of its
    1154 points with a 2*sqrt(3) in them, and dimension 12 is given in floating point
    throughout.  A SUBSET of a 60-degree code is a 60-degree code, so the integral rows are
    a legal pole source on their own and the rest are dropped; how many were dropped is
    returned rather than hidden, because it is the difference between the layer and tau(k).
    """
    cur, want, coef, pts, drop = None, False, None, [], 0
    for line in open(path):
        s = line.strip()
        if s.startswith('Dimension:'):
            if want and (pts or drop):
                break
            cur = int(s.split(':')[1])
            want, pts, coef, drop = (cur == dim), [], None, 0
        elif s.startswith('Inner product'):
            if want:
                coef = [int(x) for x in s.split(':')[1].split(',')]
        elif s.startswith('Number') or s == 'Points:' or not s:
            continue
        elif want:
            if any(ch in s for ch in '.eEsqrt'):
                drop += 1                    # a float or a radical: not a rational point
            else:
                pts.append([int(x) for x in s.split(',')])
    if not pts:
        return None, None, drop
    return np.array(pts, dtype=np.int64), np.array(coef, dtype=np.int64), drop


def source(k, path=None):
    """(X, c): the record configuration of R^k, from the package or from the data set"""
    f = os.path.join(DATA, 'pole_src_%d.npy' % k)
    g = os.path.join(DATA, 'pole_src_%d_c.npy' % k)
    if os.path.exists(f) and os.path.exists(g):
        return np.load(f).astype(np.int64), np.load(g).astype(np.int64)
    path = path or data_set()
    assert path, ("k = %d: neither data/pole_src_%d.npy nor Cohn's dimensions1-24.txt was "
                  "found; pass the data set's path or set D124" % (k, k))
    X, c, _drop = record(k, path)
    assert X is not None, 'k = %d: the data set block has no integral point' % k
    return X, c


def is_code(X, c):
    """exact, one-sided: no two distinct points of X are closer than 60 degrees"""
    Y = X * c
    d = np.einsum('ij,ij->i', Y, X)
    n = len(X)
    B = max(1, int(8e6 // max(1, n)))
    bad = 0
    for s in range(0, n, B):
        G = Y[s:s + B] @ X.T
        for r in range(len(G)):
            G[r, s + r] = -10 ** 15
        idx = np.nonzero(G > 0)
        if len(idx[0]):
            lhs = 4 * G[idx].astype(object) ** 2
            rhs = d[idx[0] + s].astype(object) * d[idx[1]].astype(object)
            bad += int((lhs > rhs).sum())
    return bad == 0, sorted(set(int(v) for v in d))


# ---------------------------------------------------------------------------- the frame
def squarefree(n):
    n, s, d = int(n), 1, 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            s *= d
        d += 1
    return s * n


def complement_line(W, V):
    """a primitive integer vector of span(W) orthogonal to every row of V

    The complement of r orthogonal cap directions inside a span of rank r + 1 is a LINE, so
    its norm is not a choice: the frame needs a particular class modulo squares and either
    the line has it or no frame of that shape exists.
    """
    B = W[span_basis(W)].astype(object)
    A = [[F(int(v)) for v in row] for row in (V.astype(np.int64) @ B.T.astype(np.int64))]
    rows_, cols = len(A), len(A[0])
    piv, r = [], 0
    for cdx in range(cols):
        p = next((i for i in range(r, rows_) if A[i][cdx]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][cdx]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows_):
            if i != r and A[i][cdx]:
                f = A[i][cdx]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(cdx)
        r += 1
    free = [j for j in range(cols) if j not in piv]
    assert len(free) == 1, 'the complement of the frame is %d-dimensional, not 1' % len(free)
    j = free[0]
    y = [F(0)] * cols
    y[j] = F(1)
    for i, cdx in enumerate(piv):
        y[cdx] = -A[i][j]
    den = 1
    for v in y:
        den = den * v.denominator // gcd(den, v.denominator)
    u = np.array([int(v * den) for v in y], dtype=object) @ B
    g = 0
    for v in u:
        g = gcd(g, abs(int(v)))
    return np.array([int(v) // g for v in u], dtype=np.int64) if g else u.astype(np.int64)


def ortho_rows(cvec, dim, cap=3):
    """integer rows of Z^dim, mutually orthogonal, row i of norm cvec[i], or None

    Used by the 'block' frame: r mutually orthogonal cap directions of norm m span the axis
    space, and A rows of A A^T = diag(c) turn them into a frame of norms m c_i.  Rows of
    norm 1 are unit vectors and can be assigned outright; the rest are found by backtracking
    in the coordinates the unit rows leave free, which at k = 23 is a six-dimensional search.
    A determinant condition decides it in advance: det(A)^2 = prod(c), so a c whose product
    is not a perfect square admits no such A at all -- which is the case at k = 22, where
    3*6*2*2*6 = 432 is not a square.
    """
    prod = 1
    for v in cvec:
        prod *= int(v)
    if isqrt(prod) ** 2 != prod:
        return None
    ones = [i for i, v in enumerate(cvec) if int(v) == 1]
    rest = [i for i in range(len(cvec)) if int(cvec[i]) != 1]
    if len(ones) > dim:
        return None
    free = list(range(len(ones), dim))
    pool = {}
    for i in rest:
        n = int(cvec[i])
        if n in pool:
            continue
        vs = []
        for combo in np.ndindex(*([2 * cap + 1] * len(free))):
            v = np.array(combo, dtype=np.int64) - cap
            if int(v @ v) == n:
                vs.append(v)
        pool[n] = vs
        if not vs:
            return None
    A = np.zeros((len(cvec), dim), dtype=np.int64)
    for t, i in enumerate(ones):
        A[i, t] = 1
    chosen = []

    def rec(t):
        if t == len(rest):
            return True
        for v in pool[int(cvec[rest[t]])]:
            if all(int(v @ w) == 0 for w in chosen):
                chosen.append(v)
                A[rest[t], len(ones):] = v
                if rec(t + 1):
                    return True
                chosen.pop()
        return False

    if not rec(0):
        return None
    assert (A @ A.T == np.diag(np.asarray(cvec, dtype=np.int64))).all()
    return A


def frame(k, W, cX, seed=0, tries=800):
    """(M, mu): M^T M = mu diag(cX), the columns spanning the axis space

    Built in closed form from mutually orthogonal cap directions; see the module docstring.
    """
    rank = len(span_basis(W))
    assert rank == len(cX), 'k = %d: rank %d against a source in %d coordinates' % (k, rank, len(cX))
    m = int(W[0] @ W[0])
    one = np.ones(W.shape[1], dtype=np.int64)
    rng = np.random.default_rng(seed)
    want = rank - (1 if PLAN[k]['tail'] is not None else 0)
    V = None
    for _ in range(tries):
        got = ortho_caps(W, one, rng, want)
        if len(got) >= want:
            V = got[:want]
            break
    assert V is not None, 'k = %d: no %d mutually orthogonal cap directions' % (k, want)
    mu, tail = PLAN[k]['mu'], PLAN[k]['tail']
    cols = []
    if PLAN[k]['shape'] == 'caps':
        cols = [r.astype(object) for r in V]
    elif PLAN[k]['shape'] == 'block':
        A = ortho_rows(cX, rank)
        assert A is not None, \
            'k = %d: no orthogonal integer block for the metric %s' % (k, list(cX))
        cols = [np.asarray(row, dtype=object) for row in (A @ V.astype(np.int64))]
    else:
        for t in range(8):
            a, b = V[2 * t].astype(object), V[2 * t + 1].astype(object)
            cols += [a + b, a - b]
        cols.append(2 * V[16].astype(object))
    if tail is not None:
        u = complement_line(W, V)
        q = int(u @ u)
        assert squarefree(q) == tail, \
            ('k = %d: the complement line has norm %d, of class %d; the frame needs %d'
             % (k, q, squarefree(q), tail))
        need = mu * int(cX[-1])
        lam = F(isqrt(need // tail), isqrt(q // tail))
        assert lam * lam * q == need, \
            'k = %d: the complement line cannot be scaled to norm %d' % (k, need)
        col = np.array([F(int(v)) * lam for v in u], dtype=object)
        cols.append(col)
    den = 1
    for col in cols:
        for x in col:
            if isinstance(x, F):
                den = den * x.denominator // gcd(den, x.denominator)
    M = np.array([[int(F(x) * den) for x in col] for col in cols], dtype=np.int64).T
    mu = mu * den * den
    G = M.T @ M
    assert (G == mu * np.diag(cX)).all(), 'k = %d: M^T M is not mu diag(c)' % k
    return M, mu


# ----------------------------------------------------------------------- the exact test
def limits(mu, m, norms, D):
    """isqrt(3 mu m |t|^2 D^2) per source norm, with the perfect-square case ruled out"""
    lim = {}
    for na in sorted(set(int(v) for v in norms)):
        L = 3 * mu * m * na * D * D
        r = isqrt(L)
        assert r * r < L < (r + 1) ** 2, '3 mu m |t|^2 D^2 is a perfect square'
        lim[na] = r
    return lim


def cleared(X, cX, N, D, Zk, mu, m):
    """boolean per source point, in exact integer arithmetic; the verifier's own test

    A pole is M N t / D and a cap is w, so the inner product is (N t) . (M^T w) over D and
    the condition 4 <a,z>^2 <= 3 |a|^2 |z|^2 is  4 X^2 <= 3 mu m |t|_c^2 D^2.
    """
    A = X @ N.T
    norms = np.einsum('ij,ij->i', X * cX, X)
    lim = limits(mu, m, norms, D)
    lo = np.array([lim[int(v)] for v in norms], dtype=np.int64)
    SH = 31
    Ahi = A >> SH
    Alo = A - (Ahi << SH)
    _bd = (int(np.abs(Ahi).max()) + 2 ** SH) * int(np.abs(Zk).max()) * N.shape[0]
    assert _bd < 2 ** 50, 'D too large: verify_record.py could not split this in float64'
    Hf, Lf, Zf = Ahi.astype(np.float64), Alo.astype(np.float64), Zk.astype(np.float64)
    good = np.empty(len(X), bool)
    worst = 0
    B = max(1, int(8e6 // max(1, len(Zk))))
    for s in range(0, len(X), B):
        XH = Hf[s:s + B] @ Zf.T
        XL = Lf[s:s + B] @ Zf.T
        assert (XH == np.rint(XH)).all() and (XL == np.rint(XL)).all(), 'not integral'
        assert int(np.abs(XH).max()) < 2 ** 31, 'D too large: the high limb would overflow'
        Y = (XH.astype(np.int64) << SH) + XL.astype(np.int64)
        w = np.abs(Y).max(1)
        worst = max(worst, int(w.max()))
        good[s:s + B] = 2 * w <= lo[s:s + B]
    assert 2 * worst < 2 ** 63 - 1, 'the doubled inner product would overflow'
    return good, worst


# ------------------------------------------------------------------- the plane rotations
def conic(t, dmax):
    """(a, b, D) with a^2 + t b^2 = D^2 and gcd(a, b, D) = 1: one plane's rational angles"""
    out = set()
    for D in range(1, dmax + 1):
        for b in range(1, D + 1):
            r = D * D - t * b * b
            if r < 0:
                break
            a = isqrt(r)
            if a * a == r and gcd(gcd(abs(a), b), D) == 1:
                out.add((a, b, D))
                if a:
                    out.add((-a, b, D))
    return sorted(out, key=lambda x: (x[2], abs(x[0])))


def givens(kk, cX, i, j, a, b, D):
    """the integer numerator of the plane rotation of diag(cX) in the (i, j) plane"""
    t = int(cX[j]) // int(cX[i])
    G = D * np.eye(kk, dtype=object)
    G[i, i] = G[j, j] = a
    G[i, j] = -t * b
    G[j, i] = b
    return G


def search(k, X, cX, Zk, mu, m, dmax=25, rounds=60, budget=10 ** 10,
           s1=260, s2=900, keep1=24, keep2=5, seed=0, verbose=True,
           slack=2, patience=6):
    """(count, N, D, keep): a greedy product of rational plane rotations

    Every round scores each (plane, angle) on the rows that are failing plus a random sample,
    re-scores the best two dozen on a larger sample, and counts the best five in full.  The
    two samples decide only WHICH candidates are counted; the count that is kept is over the
    whole configuration, and the winner is re-tested exactly by cleared() before it is
    written.
    """
    n, kk = X.shape
    rng = np.random.default_rng(seed)
    tab = {}
    for i in range(kk):
        for j in range(i + 1, kk):
            r = F(int(cX[j]), int(cX[i]))
            if r.denominator == 1 and r not in tab:
                tab[r] = conic(int(r), dmax)
    Xf = X.astype(np.float64)
    Zf = Zk.astype(np.float64)
    lim = (np.sqrt(3.0) / 2.0
           * np.sqrt(mu * m * np.einsum('ij,ij->i', X * cX, X).astype(np.float64)))[:, None]
    B = Xf.copy()                                    # B = X N^T / D
    A = B @ Zf.T
    ratio = (np.abs(A) / lim).max(1)
    N, D, cnt = np.eye(kk, dtype=object), 1, int((ratio <= 1.0).sum())
    if verbose:
        print('   k=%2d: unrotated %d of %d' % (k, cnt, n), flush=True)
    hist, seen, stall = [], [(cnt, 1)], 0
    for _rd in range(rounds):
        bad = np.nonzero(ratio > 1.0)[0]
        r1 = np.unique(np.concatenate([bad, rng.choice(n, size=min(s1, n), replace=False)]))
        r2 = np.unique(np.concatenate([bad, rng.choice(n, size=min(s2, n), replace=False)]))
        pool = []
        for i in range(kk):
            for j in range(i + 1, kk):
                r = F(int(cX[j]), int(cX[i]))
                if r not in tab:
                    continue
                t = int(r)
                ab = [x for x in tab[r] if D * x[2] <= budget]
                if not ab:
                    continue
                Ar, Lr = A[r1], lim[r1]
                bi, bj = B[r1, i], B[r1, j]
                zi, zj = Zf[:, i], Zf[:, j]
                base = Ar - np.outer(bi, zi) - np.outer(bj, zj)
                for (a, b, d) in ab:
                    ni = (a * bi - t * b * bj) / d
                    nj = (b * bi + a * bj) / d
                    Y = np.abs(base + np.outer(ni, zi) + np.outer(nj, zj)) / Lr
                    pool.append((int((Y.max(1) <= 1.0).sum()), i, j, a, b, d))
        pool.sort(key=lambda c: (-c[0], c[5]))
        mid = []
        for _sc, i, j, a, b, d in pool[:keep1]:
            t = int(F(int(cX[j]), int(cX[i])))
            Ar, Lr = A[r2], lim[r2]
            bi, bj = B[r2, i], B[r2, j]
            zi, zj = Zf[:, i], Zf[:, j]
            ni = (a * bi - t * b * bj) / d
            nj = (b * bi + a * bj) / d
            Y = np.abs(Ar - np.outer(bi, zi) - np.outer(bj, zj)
                       + np.outer(ni, zi) + np.outer(nj, zj)) / Lr
            mid.append((int((Y.max(1) <= 1.0).sum()), i, j, a, b, d))
        mid.sort(key=lambda c: (-c[0], c[5]))
        best, alt = (cnt, None), []
        for _sc, i, j, a, b, d in mid[:keep2]:
            t = int(F(int(cX[j]), int(cX[i])))
            ni = (a * B[:, i] - t * b * B[:, j]) / d
            nj = (b * B[:, i] + a * B[:, j]) / d
            Ad = A + np.outer(ni - B[:, i], Zf[:, i]) + np.outer(nj - B[:, j], Zf[:, j])
            rt = (np.abs(Ad) / lim).max(1)
            full = int((rt <= 1.0).sum())
            if full > best[0]:
                best = (full, (i, j, a, b, d, ni, nj, Ad, rt))
            elif full >= cnt - slack:
                alt.append((full, (i, j, a, b, d, ni, nj, Ad, rt)))
        # THE PLATEAU.  Stopping at the first round where no single rotation improves is
        # what a greedy does, and it is not what the problem asks: two rotations can gain
        # where neither gains alone.  So when nothing improves, take a sideways move -- the
        # best of the candidates within `slack` of the current count -- and go on, keeping
        # the best PREFIX of the accepted sequence rather than its end.  The denominator is
        # what limits this, not the round count: every accepted rotation multiplies it.
        if cnt >= n:
            break                       # the whole source clears: there is nothing above it
        if best[1] is None:
            if not alt or stall >= patience:
                break
            stall += 1
            best = max(alt, key=lambda c: (c[0], -c[1][4]))
        else:
            stall = 0
        i, j, a, b, d, ni, nj, Ad, rt = best[1]
        B[:, i], B[:, j], A, ratio, cnt = ni, nj, Ad, rt, best[0]
        N = givens(kk, cX, i, j, a, b, d) @ N
        D *= d
        g = D
        for row in N:
            for v in row:
                g = gcd(g, abs(int(v)))
                if g == 1:
                    break
            if g == 1:
                break
        if g > 1:
            N, D = N // g, D // g
        hist.append((i, j, a, b, d))
        seen.append((cnt, D))
        if verbose:
            print('      (%2d,%2d) a=%4d b=%3d D=%-3d -> %5d of %5d   denominator %-12d%s'
                  % (i, j, a, b, d, cnt, n, D, '  (sideways)' if stall else ''), flush=True)
    # the best PREFIX, not the last state: a sideways move can be followed by a worse one
    _bi = max(range(len(seen)), key=lambda t: (seen[t][0], -seen[t][1]))
    if _bi < len(hist):
        N, D = np.eye(kk, dtype=object), 1
        for (i, j, a, b, d) in hist[:_bi]:
            N = givens(kk, cX, i, j, a, b, d) @ N
            D *= d
        g = D
        for row in N:
            for v in row:
                g = gcd(g, abs(int(v)))
        if g > 1:
            N, D = N // g, D // g
        cnt, hist = seen[_bi][0], hist[:_bi]
        if verbose:
            print('      keeping the best prefix: %d rotation(s), %d of %d, denominator %d'
                  % (_bi, cnt, n, D), flush=True)
    Ni = np.array(N.tolist(), dtype=np.int64)
    assert (Ni.astype(object) == N).all(), 'k = %d: N does not fit in int64' % k
    assert ((Ni.astype(object).T @ np.diag(cX.astype(object)) @ Ni.astype(object))
            == (D * D) * np.diag(cX.astype(object))).all(), 'N is not an isometry of diag(c)'
    good, _w = cleared(X, cX, Ni, D, Zk, mu, m)
    ex = int(good.sum())
    assert ex == cnt, ('k = %d: the float screen and the exact test disagree (%d vs %d)'
                       % (k, cnt, ex))
    return ex, Ni, D, np.nonzero(good)[0].astype(np.int32), hist


# ------------------------------------------------------------------------------ the pass
def shipped(k):
    """the largest layer already on disk for this k, over the file kinds that can hold one"""
    best = 0
    for fn in ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k, 'poles_rot_%d.npz' % k):
        p = os.path.join(DATA, fn)
        if os.path.exists(p):
            best = max(best, len(np.load(p, allow_pickle=True)['keep']))
    p = os.path.join(DATA, 'poles_%d.npy' % k)
    if os.path.exists(p):
        best = max(best, len(np.load(p)))
    return best


def main(write=False, ks=None, dmax=25, rounds=60, seed=0, src=None):
    print(__doc__)
    ROWS = rows()
    ds = data_set([src] if src else [])
    print('   Cohn\'s dimensions1-24.txt: %s' % (ds or 'not found -- using data/pole_src_*.npy'))
    print('=' * 92)
    print('  k  dim    |Z|   source  shipped   found    gain   mu    D          at tau(k)?')
    print('=' * 92)
    tot = 0
    for k in [k for k in KS if ks is None or k in ks]:
        t0 = time.time()
        W, c, m = directions(k, ROWS)
        assert (c == 1).all(), 'k = %d: the ambient metric is not the identity' % k
        X, cX = source(k, ds)
        assert len(X) <= TAU[k], \
            'k = %d: the source has %d points, more than tau(k) = %d' % (k, len(X), TAU[k])
        if len(X) < TAU[k]:
            print('   k=%2d: %d of the %d points of the record configuration are rational;'
                  ' the other %d are not usable' % (k, len(X), TAU[k], TAU[k] - len(X)))
        ok, norms = is_code(X, cX)
        assert ok, 'k = %d: the source is not a 60-degree code' % k
        M, mu = frame(k, W, cX, seed=seed)
        Zk = W @ M
        cnt, N, D, keep, hist = search(k, X, cX, Zk, mu, m, dmax=dmax, rounds=rounds, seed=seed)
        was = shipped(k)
        tot += max(0, cnt - was)
        print('  %2d  %3d %6d %8d %8d %7d %+7d %5d  %-10d %s   (%.0fs)'
              % (k, 72 + k, len(W), len(X), was, cnt, cnt - was, mu, D,
                 'YES' if cnt == TAU[k] else '%d short' % (TAU[k] - cnt), time.time() - t0))
        if write and cnt > was:
            np.save(os.path.join(DATA, 'pole_src_%d.npy' % k), X)
            np.save(os.path.join(DATA, 'pole_src_%d_c.npy' % k), cX)
            np.savez(os.path.join(DATA, 'poles_rec_%d.npz' % k),
                     N=N, D=np.array([D]), keep=keep, M=M, cs=cX,
                     code=np.zeros(0, dtype=np.int64))
            print('      wrote data/poles_rec_%d.npz and data/pole_src_%d.npy'
                  ' -- run scripts/seal_record.py' % (k, k))
    print('=' * 92)
    print('  total gain over dimensions %s: %+d spheres'
          % (', '.join(str(72 + k) for k in KS), tot))
    if not write:
        print('  nothing written: pass --write')
    return 0


if __name__ == '__main__':
    _ks = tuple(int(x) for x in sys.argv[1:] if x.isdigit()) or None
    sys.exit(main(write='--write' in sys.argv, ks=_ks))
