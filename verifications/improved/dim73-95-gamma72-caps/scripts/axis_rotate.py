#!/usr/bin/env python3
"""The axis layers of dimensions 81-95, as ROTATED copies of the cap directions.

    python axis_rotate.py                re-derive and report against the shipped layers
    python axis_rotate.py --write        re-derive and write data/poles_rot_<k>.npz
    python axis_rotate.py --write 7 10   just those k, cross-sections and root systems alike

WHAT AN AXIS LAYER HAS TO BE.  At cap level t = 2/3 a pole (0,a) meets a cap at
sqrt(1/3)<z,a>, so it needs 4<a,z>^2 <= 3|a|^2|z|^2 -- at least 30 degrees -- against every
cap direction z, and the poles must be a 60-degree code among themselves.

WHY A ROTATION.  Take the poles to be a ROTATED COPY of the cap direction set itself.  The
60-degree condition is then free, because a rotation preserves inner products, and the whole
question becomes how many of the |W| rotated directions clear the 30-degree caps.  That is
the same question dimension 38 asks, and the answer there was "all of them".

Here it is easier, because in these dimensions the caps are small.  A 30-degree cap covers
2.3e-3 of S^22, so a random rotation of the 93 150 directions of Lambda_23 is expected to
lose about two hundred of them -- against an axis layer of 248 points before this file
existed.  The layers built by the earlier shell rule ran at 0.27% of tau(k) at k = 23 and 40%
at k = 12; the rotated copies run at 99% and above, and reach tau(k) exactly wherever
Lambda_k does.

THE ROTATION IS RATIONAL, so the verification stays exact integer arithmetic.  For the metric
C = diag(c) -- c = 1 everywhere except k = 10, where the record configuration carries
diag(1,1,1,1,1,3,3,3,3,3) -- take K an INTEGER skew matrix with entries in {-1,0,1} and put

    S = C^{-1} K,        Q = (I - S)(I + S)^{-1} .

CS = K is skew, so Q^T C Q = C exactly: Q is an isometry of the form, it is rational, and
Hadamard bounds its denominator by k^{k/2}, which is 1.4e15 at k = 23.  Everything the check
touches then stays inside int64.

THAT CHART HAS A HOLE IN IT, and it is where the good rotations live.  I + S = 2(I + Q)^{-1},
so S runs to infinity exactly as Q acquires an eigenvalue of -1 -- and a rotation that lifts
the whole direction set clear of the caps is one that moves it a long way.  See charts():
rationalising Q E for E an integral isometry of the form puts the same rotation somewhere
finite, and E^T C E = C keeps N = N_0 E integral with N^T C N = D^2 C.  Without it k = 10
stalls at 464 of 510 for every denominator between 8 and 65536; with it, at 510.

THE TEST, exactly.  With Q = N/D and X = (N w)^T (C z), an integer,

    4<a,z>^2 <= 3 m^2   <=>   4X^2 <= 3 m^2 D^2   <=>   2|X| <= isqrt(3 m^2 D^2),

the last step because 3 m^2 D^2 is never a perfect square.  One big-integer square root per
k, computed once; the 1932^2 ... 93150^2 comparisons are int64.

WHAT IS HEURISTIC AND WHAT IS NOT.  Choosing K is a search over random integer skew matrices
and can come out differently on a different machine.  Nothing rests on it: verify_poles.py
takes the SHIPPED N, D and index set and checks every pair in exact integers.

k = 11 is not treated: dimension 83 uses the 604-point record realised over Z[sqrt2], whose
Gram is irrational, and the construction above is set up for rational metrics only.  Its axis
layer stays empty, as it was.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from published import TAU                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import capsha                                                      # noqa: E402
PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
KS = [k for k in range(9, 24) if k != 11]
SAMPLES_ROOT = 20000     # tries per root system; a try is one Fraction elimination


# SUPERSEDED FOR k = 13, 15, 17, 18, 21, 22 AND 23 (2026-08-27).  Those direction sets are
# written in more coordinates than they span -- Lambda_21's cross-section is 27 720 vectors in
# Lambda_24's 24 coordinates spanning 21 -- and the Cayley transform below is an isometry of
# the coordinates, not of the axis space, so it carries the poles out of R^k.  preserves_span
# now rejects such a candidate, so this file produces nothing at those k and they are derived
# by scripts/axis_block.py, whose rotation is built from the direction set itself and fixes
# the C-orthogonal complement pointwise.


def rows():
    return {r['k']: r for r in json.load(open(os.path.join(PKG, 'rows81-95.json')))}


def directions(k, ROWS):
    """(W, c, m): the cap directions dimension 72+k uses, their metric, and their norm"""
    if ROWS[k]['W'].startswith('record'):
        W = np.load(os.path.join(DATA, 'rec_%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(DATA, 'rec_%d_c.npy' % k)).astype(np.int64)
    elif ROWS[k]['W'] == 'kappa':
        W = np.load(os.path.join(DATA, 'kap%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(DATA, 'kap%d_c.npy' % k)).astype(np.int64)
    else:
        f = 'lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k
        W = np.load(os.path.join(DATA, f)).astype(np.int64)
        c = np.ones(W.shape[1], np.int64)
    keep = (np.abs(W).sum(0) > 0) | (c != 1)          # the arrays are zero-padded to 16 or 24
    W, c = W[:, keep], c[keep]
    m = int((W[0] * c) @ W[0])
    assert (np.einsum('ij,ij->i', W * c, W) == m).all(), 'the directions are not all norm m'
    return W, c, m


def _solve(G, K):
    """S = G^{-1} K in exact rationals, G an integer Gram (diagonal or not)"""
    n = len(G)
    aug = [[F(int(G[i][j])) for j in range(n)] + [F(int(K[i][j])) for j in range(n)]
           for i in range(n)]
    for col in range(n):
        p = next((r for r in range(col, n) if aug[r][col]), None)
        if p is None:
            return None
        aug[col], aug[p] = aug[p], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col]:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def cayley(K, G):
    """(N, D) with N/D = (I - S)(I + S)^{-1}, S = G^{-1} K, in exact rationals.

    CS = K is skew, which is exactly what makes Q^T G Q = G: the poles are then a copy of the
    directions under an isometry OF THE FORM, whatever the form is.  G is diag(c) for the
    cross-sections and the Cartan matrix for the root systems."""
    n = len(G)
    # Python ints, not numpy ones: Fraction accepts a numpy int64 and then carries it into
    # .denominator, so D comes back as an int64 and silently overflows at k = 10, the one
    # metric that is not the identity.
    G = [[int(v) for v in row] for row in G]
    K = [[int(v) for v in row] for row in K]
    S = _solve(G, K)
    if S is None:
        return None, None
    A = [[-S[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
    aug = [[S[i][j] + (1 if i == j else 0) for j in range(n)]
           + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    for col in range(n):
        p = next((r for r in range(col, n) if aug[r][col]), None)
        if p is None:
            return None, None                          # I + S singular; try another K
        aug[col], aug[p] = aug[p], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col]:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    Binv = [row[n:] for row in aug]
    X = [[sum((A[i][t] * Binv[t][j] for t in range(n)), F(0)) for j in range(n)]
         for i in range(n)]
    D = 1
    for row in X:
        for x in row:
            D = D * x.denominator // gcd(D, x.denominator)
    return np.array([[int(x * D) for x in row] for row in X], dtype=object), int(D)


def screen(Wf, CZf, m, N, D):
    """how many clear the caps, in floating point -- used to RANK candidates only.

    numpy has no BLAS for int64, so the exact pass costs a hundred times what the float one
    does at k = 23.  The search is a search either way; what makes the result a certificate
    is the exact pass below, run once on the winner, and verify_poles.py run on the file."""
    A = Wf @ np.array(N.tolist(), dtype=np.float64).T
    thr = np.sqrt(3.0) / 2.0 * m * D
    n = len(Wf)
    good = np.empty(n, bool)
    B = max(1, int(8e6 // max(1, n)))
    for s in range(0, n, B):
        good[s:s + B] = np.abs(A[s:s + B] @ CZf).max(1) <= thr
    return int(good.sum())


def cleared(W, c, m, N, D):
    """boolean per rotated direction, in EXACT integer arithmetic"""
    m, D = int(m), int(D)
    Ni = np.array(N.tolist(), dtype=np.int64)
    assert (Ni.astype(object) == N).all(), 'N does not fit in int64'
    A = W @ Ni.T                                    # D * (rotated directions)
    CZ = (W * c).T
    # The two bounds scripts/verify_poles.py asserts while splitting A at bit 31.  A layer
    # that fails them is one the certificate script cannot read, so it is rejected here
    # rather than written and discovered later -- rationalise() treats the AssertionError as
    # "this candidate does not qualify" and moves on.  2mD < 2^62 alone does NOT imply them,
    # and neither implies the other: at k = 3 the verifier stops at D of about 3e10 against
    # the search's 1e18, while at k = 10 and 14 the search's guard is the tighter of the two.
    _Ahi = A >> 31
    _bd = (int(np.abs(_Ahi).max()) + 2 ** 31) * int(np.abs(CZ).max()) * W.shape[1]
    assert _bd < 2 ** 50, 'D too large: verify_poles.py could not split this in float64'
    assert int(np.abs(_Ahi.astype(np.float64) @ CZ.astype(np.float64)).max()) < 2 ** 31, \
        'D too large: verify_poles.py could not shift the high limb'
    lim = isqrt(3 * m * m * D * D)
    assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2      # 3 m^2 D^2 is not a square
    good = np.empty(len(W), bool)
    worst = 0
    B = max(1, int(8e6 // max(1, len(W))))
    for s in range(0, len(W), B):
        X = A[s:s + B] @ CZ
        mx = np.abs(X).max(1)
        worst = max(worst, int(mx.max()))
        good[s:s + B] = 2 * mx <= lim
    assert 2 * worst < 2 ** 63 - 1, 'the int64 product could overflow'
    return good, lim, worst


def layer_file(k, ROWS):
    """the file this k's rotated layer lives in

    A layer is indices into W and a rotation OF W, so it belongs to a cap set and not to a
    dimension.  When the cap set at k = 12 and 13 became the Kappa section, the old
    poles_rot_<k>.npz stayed correct for Lambda_k and became meaningless for K_k -- it
    indexes vectors that are no longer there -- so the two are kept apart by name rather
    than overwritten."""
    return ('poles_kap_%d.npz' if ROWS[k]['W'] == 'kappa' else 'poles_rot_%d.npz') % k


def shipped_rotation(k, fn=None):
    """the rotation already in data/, as a starting point -- see borrowed() for the same idea

    Starting the descent from the best rotation known makes this file monotone: it cannot
    return less than what is shipped.  At k = 9 that is the difference between 286 poles and
    all 306."""
    for f in (fn or 'poles_rot_%d.npz' % k, 'poles_rot_k%d.npz' % k):
        p = os.path.join(DATA, f)
        if os.path.exists(p):
            z = np.load(p, allow_pickle=True)
            return np.array(z['N'].tolist(), dtype=object), int(z['D'][0])
    return None, None


def borrowed(k, W):
    """rotations already derived elsewhere in the repository for the SAME direction set.

    k = 14's cap directions are exactly dimension 38's 1932-point code, vector for vector, and
    dimension 38's rotation was pushed to all 1932 by DESCENT on SO(14) where sampling stalls
    at 1908.  So it is offered here first; the identity of the two sets is checked, not
    assumed."""
    if k != 14:
        return []
    d38 = os.path.join(HERE, '..', '..', 'dim38-leech-large-codimension', 'data')
    nf = os.path.join(d38, 'axis_rotation.txt')
    df = os.path.join(d38, 'axis_rotation_D.txt')
    wf = os.path.join(d38, 'dim14_directions.txt')
    if not all(os.path.exists(p) for p in (nf, df, wf)):
        return []
    A = np.loadtxt(wf, delimiter=',', dtype=np.int64)
    if {tuple(r) for r in A.tolist()} != {tuple(r) for r in W.tolist()}:
        return []
    return [(np.loadtxt(nf, dtype=np.int64).astype(object), int(open(df).read().strip()))]


def span_basis(W):
    """indices of rows of W spanning its row space, by fraction-free elimination

    The gcd step is not cosmetic: without it the entries double at every round and the
    24-coordinate direction sets produce numbers thousands of digits long."""
    rows_, piv, idx = [], [], []
    A = W.shape[1]
    for i in range(len(W)):
        r = [int(x) for x in W[i]]
        for p, b in zip(piv, rows_):
            if r[p]:
                r = [b[p] * x - r[p] * y for x, y in zip(r, b)]
                g = 0
                for x in r:
                    g = gcd(g, x)
                if g > 1:
                    r = [x // g for x in r]
        p = next((j for j in range(A) if r[j]), None)
        if p is None:
            continue
        piv.append(p)
        rows_.append(r)
        idx.append(i)
        if len(idx) == A:
            break
    return idx


def preserves_span(W, N):
    """Does N map the span of the directions into itself?

    THE ROTATION MUST BE AN ISOMETRY OF THE AXIS SPACE.  This file builds one on the
    coordinates W is written in, which is the axis space only when W has full rank there.
    At k = 13, 15, 17, 18, 21, 22 and 23 it does not -- Lambda_21's cross-section is 27 720
    vectors in Lambda_24's 24 coordinates spanning 21 -- and the rotation carries the poles
    into dimensions the configuration does not have.  Seven layers were written that way on
    2026-08-27 with every inner product correct.  Q is linear, so a basis settles it."""
    idx = span_basis(W)
    B = [[int(v) for v in W[i]] for i in idx]
    xs = [[int(v) for v in (W[i].astype(object) @ N.T)] for i in idx]
    k, A = len(B), len(B[0])
    M = [[F(B[r][col]) for r in range(k)] for col in range(A)]
    rhs = [[F(x[col]) for x in xs] for col in range(A)]
    r = 0
    for col in range(k):
        p = next((i for i in range(r, A) if M[i][col]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        rhs[r], rhs[p] = rhs[p], rhs[r]
        pv = M[r][col]
        M[r] = [v / pv for v in M[r]]
        rhs[r] = [v / pv for v in rhs[r]]
        for i in range(A):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
                rhs[i] = [a - f * b for a, b in zip(rhs[i], rhs[r])]
        r += 1
    if r != k:
        return False
    return all(rhs[i][t] == 0 for i in range(r, A) for t in range(len(xs)))


def check_orthogonal(N, D, c):
    n = len(c)
    C = np.zeros((n, n), dtype=object)
    for t in range(n):
        C[t, t] = int(c[t])
    # N^T C N = D^2 C, NOT N C N^T: the two coincide only when C is the identity, which is
    # why this passed at k = 9 and failed at k = 10, the one metric that is not.
    return ((N.T @ C @ N) == (D * D) * C).all()


# ---- the root systems of dimensions 75-80, through their Gram --------------------------
# A root is its coordinate vector in the SIMPLE-ROOT basis and the Gram is the Cartan matrix
# with 2 on the diagonal, so E_6 and E_7 need no embedding in R^8 and every quantity stays
# integral.  A rotation is then a G-isometry, which is what cayley() already builds.
CARTAN = {
    'A_3': [(0, 1), (1, 2)],
    'D_4': [(0, 1), (1, 2), (1, 3)],
    'D_5': [(0, 1), (1, 2), (2, 3), (2, 4)],
    'E_6': [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5)],
    'E_7': [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6)],
    'E_8': [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)],
}
ROOTNAME = {3: 'A_3', 4: 'D_4', 5: 'D_5', 6: 'E_6', 7: 'E_7', 8: 'E_8'}
ROOTTAU = {3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240}


def root_system(k):
    """(R, G): every root as a coordinate vector, and the Cartan matrix, both integral.

    Built by closing the simple roots under simple reflections and then SORTED, so the order
    is reproducible -- a shipped index set means nothing otherwise."""
    n = k
    G = 2 * np.eye(n, dtype=np.int64)
    for i, j in CARTAN[ROOTNAME[k]]:
        G[i, j] = G[j, i] = -1
    seen = set()
    for r in np.eye(n, dtype=np.int64):
        seen.add(tuple(int(x) for x in r))
        seen.add(tuple(-int(x) for x in r))
    frontier = list(seen)
    while frontier:
        new = []
        for r in frontier:
            rv = np.array(r, dtype=np.int64)
            for i in range(n):
                e = np.zeros(n, dtype=np.int64)
                e[i] = 1
                t = tuple(int(x) for x in rv - (rv @ G[i]) * e)
                if t not in seen:
                    seen.add(t)
                    new.append(t)
        frontier = new
    R = np.array(sorted(seen), dtype=np.int64)
    assert len(R) == ROOTTAU[k], (k, len(R))
    assert (np.einsum('ij,jk,ik->i', R, G, R) == 2).all(), k
    off = R @ G @ R.T
    np.fill_diagonal(off, -2)
    assert int(off.max()) <= 1, 'k=%d: the roots are not a 60-degree code' % k
    return R, G


# ------------------------------------------------------------------ descent, after sampling
# Sampling stalls: 70 of 72 at k = 6, 280 of 306 at k = 9, 208 of 240 at k = 8.  Descent
# reaches all of them.  This is dimension 38's search (dim38-.../scripts/axis_rotation.py)
# with the metric put back in: C = L L^T and v = L^T w turn <.,.>_C into the Euclidean form,
# an orthogonal O acts there, and Q = L^{-T} O L^T is the C-isometry it corresponds to.


def _euclid(C):
    """(V-map, Q-from-O, O-from-Q) for the metric C"""
    L = np.linalg.cholesky(C.astype(np.float64))
    LT = L.T
    LTi = np.linalg.inv(LT)
    return L, (lambda O: LTi @ O @ LT), (lambda Q: LT @ Q @ LTi)


def _band(V, O, thr, delta, block):
    """the (i, j) pairs within delta of the threshold -- the only ones with a gradient"""
    n = len(V)
    OV = V @ O.T
    ri, ci = [], []
    for s in range(0, n, block):
        r, c = np.nonzero(OV[s:s + block] @ V.T > thr - delta)
        if len(r):
            ri.append(r + s)
            ci.append(c)
    if not ri:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(ri), np.concatenate(ci)


def _cleared_float(V, O, thr, block):
    n, bad = len(V), 0
    OV = V @ O.T
    for s in range(0, n, block):
        bad += int(((OV[s:s + block] @ V.T).max(1) > thr).sum())
    return n - bad


MINIMAX_MAX = 2500      # |W| above which the dense |W|^2 objective is not affordable


def minimax(V, O, thr, steps=6000, eta=0.03, beta=200.0):
    """push the WORST cosine down, instead of the count up; (worst cosine, best O)

    descend() below maximises HOW MANY directions clear the caps.  That objective goes flat
    once nearly all of them do -- a softmax over the caps sees one bad row against hundreds
    of clear ones -- and it stalled two short of tau(7) and forty-six short of tau(10).
    Minimising the largest |<O v_i, v_j>| has a gradient everywhere, and with it BOTH
    ceilings turn out to be reachable in floating point with room to spare: the worst
    cosine comes down to 0.953 of the threshold at k = 7 and 0.977 at k = 10.  So the whole
    of that shortfall was in the RATIONAL APPROXIMATION and none of it in the search --
    which is why it went unfixed, the count objective reporting a plateau that was not there.

    The cost is |V|^2 per step with no band, so this is for the small direction sets.  The
    large ones keep descend(), where clearing all of W is not on offer in any case."""
    best, bO = float(np.abs((V @ O.T) @ V.T).max()), O.copy()
    for t in range(steps):
        X = (V @ O.T) @ V.T
        A = np.abs(X)
        P = np.exp(beta * (A - A.max()) / thr)          # a soft argmax over ALL pairs
        P /= P.sum()
        M = (V.T @ (P * np.sign(X)).T @ V) @ O.T        # d/dO of sum P_ij <O v_i, v_j>
        M = 0.5 * (M - M.T)
        w, U = np.linalg.eigh(1j * (-eta * M))          # expm of a real skew matrix
        O = np.real(U @ np.diag(np.exp(-1j * w)) @ U.conj().T) @ O
        u, _s, vt = np.linalg.svd(O)
        O = u @ vt
        g = float(np.abs((V @ O.T) @ V.T).max())
        if g < best:
            best, bO = g, O.copy()
        if t and t % 2000 == 0:                         # anneal: finer steps, sharper max
            eta *= 0.5
            beta *= 2.0
    return best, bO


def descend(V, O, thr, m, margin, steps, block, rescan=25):
    """maximise the number of rows of O V^T clearing thr; returns (best count, best O)"""
    n = len(V)
    # The band is what makes this affordable at k = 23, where the score matrix is 8.7e9
    # entries.  It is also what makes it WORSE than a dense descent at small k, where the
    # caps are large and most pairs carry gradient: 286 poles against 306 at k = 9.  So the
    # band is only narrow where it has to be.
    beta, eta, tgt = 40.0 / m, 0.35 / m, thr - margin
    delta = 1.0 * m if n <= 3000 else (0.6 * m if n <= 20000 else 0.30 * m)
    if n <= 3000:
        rescan = 5
    best, bestO = _cleared_float(V, O, thr, block), O.copy()
    rows = rmap = ri = ci = None
    for t in range(steps):
        if t % rescan == 0:
            ri, ci = _band(V, O, thr, delta, block)
            if len(ri) == 0:
                break
            rows = np.unique(ri)
            idx = dict((int(r), i) for i, r in enumerate(rows))
            rmap = np.array([idx[int(r)] for r in ri])
        s = (V[ri] * (V[ci] @ O)).sum(1)
        nr = len(rows)
        mx = np.full(nr, -1e18)
        np.maximum.at(mx, rmap, s)
        e = np.exp(beta * (s - mx[rmap]))
        den = np.zeros(nr)
        np.add.at(den, rmap, e)
        p = e / den[rmap]
        sm = np.zeros(nr)
        np.add.at(sm, rmap, p * s)
        g = 1.0 / (1.0 + np.exp(-(sm - tgt)))                      # softplus'
        a = p * (1.0 + beta * (s - sm[rmap])) * g[rmap]
        M = ((V[ci] * a[:, None]).T @ V[ri]) @ O.T
        M = 0.5 * (M - M.T)
        w, U = np.linalg.eigh(1j * (-eta * M))                     # expm of a real skew
        O = np.real(U @ np.diag(np.exp(-1j * w)) @ U.conj().T) @ O
        u, _sv, vt = np.linalg.svd(O)
        O = u @ vt
        if t % rescan == rescan - 1:
            c = _cleared_float(V, O, thr, block)
            if c > best:
                best, bestO = c, O.copy()
                if best == n:
                    break
    return best, bestO


DENOMS = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 28, 32, 40, 48, 64, 80, 96,
          128, 160, 192, 256)


def charts(C, rng, count=24):
    """integral E with E^T C E = C -- rational isometries to move the Cayley chart by

    THE CAYLEY TRANSFORM DOES NOT COVER O(k).  Q = (I - S)(I + S)^{-1} has
    I + S = 2 (I + Q)^{-1}, so S runs to infinity exactly as Q acquires an eigenvalue of -1
    -- and a rotation that carries the whole direction set clear of the caps is precisely
    one that moves it a long way.  At k = 10 the minimax optimum sits essentially at that
    singularity: |S| = 7.4e15, and every denominator from 8 to 65536 rounded to the SAME
    matrix, 434 of 510.  Reading that as a denominator ceiling was wrong; it was the chart.

    Rationalising Q E instead puts the same rotation at a different point of the chart.
    E is integral with E^T C E = C, so N = N_0 E is integral and N^T C N = D^2 C survives.
    With this, k = 7 reaches tau(7) = 126 and k = 10 reaches tau(10) = 510.

    For a diagonal form the sign matrices are isometries.  For a Cartan matrix the Weyl
    group is, generated by the simple reflections s_i = I - e_i (G row i): with roots as
    coordinate vectors in the simple-root basis, <r, alpha_i^dual> = (G r)_i.  Both are
    checked here rather than assumed."""
    k = len(C)
    out = [np.eye(k, dtype=np.int64)]
    if np.count_nonzero(C - np.diag(np.diag(C))) == 0:
        seen = {tuple([1] * k)}
        for _ in range(count * 8):
            e = tuple(int(x) for x in rng.choice([-1, 1], k))
            if e not in seen:
                seen.add(e)
                out.append(np.diag(np.array(e, dtype=np.int64)))
            if len(out) > count:
                break
    else:
        eye = np.eye(k, dtype=np.int64)
        refl = [eye - np.outer(eye[i], C[i]) for i in range(k)]
        for _ in range(count):
            E = eye.copy()
            for _t in range(int(rng.integers(1, 2 * k))):
                E = E @ refl[int(rng.integers(0, k))]
            out.append(E)
    for E in out:
        assert (E.T @ C @ E == C).all(), 'the chart matrix is not an isometry of the form'
    return out


def rationalise(Q, C, W, CZ, m, exact, rng=None, count=24, top=None, keep=4):
    """the best EXACT layer over (chart, denominator); (count, N, D)

    W and CZ are the direction set and (C w_j) as floats: the (chart, denominator) pairs are
    RANKED in floating point and only the best few are counted exactly, because there are
    now a hundred of them rather than twenty-six.  That is the same division of labour the
    rest of this file uses -- the float pass decides nothing, and verify_poles.py re-checks
    the file that is written."""
    k = len(C)
    I = np.eye(k)
    Ci = np.rint(C).astype(np.int64)
    assert (Ci == C).all(), 'the form must be integral'
    thr = np.sqrt(3.0) / 2.0 * m
    if rng is None:
        rng = np.random.default_rng(0)
    if top is None:
        # One Fraction elimination per (chart, denominator), and they are k x k, so this is
        # the file's whole runtime above k = 16.  The chart sweep earns its keep where the
        # MINIMAX objective runs, since that is what pushes a rotation to the eigenvalue-(-1)
        # singularity the Cayley chart cannot represent; the count objective does not, and
        # above MINIMAX_MAX the minimax does not run.  So the large k keep the identity chart
        # alone and the cost they had before this was added.
        top = 4 if k <= 12 else (2 if len(W) <= MINIMAX_MAX else 1)
    cand = []
    for E in charts(Ci, rng, count):
        QE = Q @ E.astype(np.float64)
        cand.append((float(np.linalg.svd(I + QE, compute_uv=False)[-1]), E, QE))
    cand.sort(key=lambda t: -t[0])                 # the best-conditioned charts first
    # The IDENTITY chart goes in whatever its conditioning, so that this search contains
    # the one it replaces and cannot return less than it.  charts() returns it first.
    ident = next(t for t in cand if (t[1] == np.eye(k, dtype=np.int64)).all())
    use = [ident] + [t for t in cand if t is not ident][:max(0, top - 1)]
    pool = []
    for smin, E, QE in use:
        if smin < 1e-9:
            continue
        try:
            S = np.linalg.solve((I + QE).T, (I - QE).T).T
        except np.linalg.LinAlgError:
            continue
        KR = C @ S
        KR = 0.5 * (KR - KR.T)
        if not np.isfinite(KR).all():              # an exactly singular chart
            continue
        nW = len(W)
        blk = max(1, int(8e6 // max(1, nW)))       # the n x n score matrix is 2.4 GB at k=20
        for d in DENOMS:
            if d * np.abs(KR).max() > 9e15:        # outside int64 before rounding
                continue
            Ki = np.rint(d * KR).astype(np.int64)
            Ki = np.triu(Ki, 1)
            N0, D = cayley((Ki - Ki.T).tolist(), (d * Ci).astype(np.int64))
            if N0 is None or 2 * m * D >= 2 ** 62:
                continue
            N = N0 @ E.astype(object)
            Qf = np.array([[float(int(v)) for v in row] for row in N]) / D
            QW = W @ Qf.T
            sc = sum(int((np.abs(QW[t:t + blk] @ CZ.T).max(1) <= thr).sum())
                     for t in range(0, nW, blk))
            pool.append((sc, D, N))
    pool.sort(key=lambda t: (-t[0], t[1]))         # most cleared, then smallest denominator
    best = (0, None, None)
    for _sc, D, N in pool[:keep]:
        try:
            cnt = exact(N, D)
        except (AssertionError, OverflowError):
            # OverflowError is numpy refusing to put N in int64.  2mD < 2^62 above bounds
            # |N| well inside it, so this is a belt-and-braces branch, not a live one.
            continue
        if cnt >= 0 and (cnt > best[0] or (cnt == best[0] and cnt > 0 and D < best[2])):
            best = (cnt, N, D)
    return best


def improve(W, C, m, N, D, cnt, exact, rng, starts=None, steps=900, margins=None):
    """descend from the best sampled rotation, and from a few random ones; keep the best

    The effort is scaled to the size of the direction set: the small ones are cheap and the
    margin that survives rationalisation varies with k, so they get a wider sweep."""
    n = len(W)
    if starts is None:
        starts = 24 if n < 300 else (8 if n < 1500 else (4 if n < 12000 else 2))
    if margins is None:
        margins = ((0.02, 0.05, 0.09, 0.14, 0.20) if n < 1500 else (0.02, 0.06, 0.12))
    L, q_of_o, o_of_q = _euclid(C)
    V = W.astype(np.float64) @ L
    thr = np.sqrt(3.0) / 2.0 * m
    block = max(1, int(4e7 // n))
    Cf = C.astype(np.float64)
    # (w_i) and (C w_j) as floats, for rationalise() to rank candidates with.  C w_j is the
    # diagonal product for the cross-sections and G r_j for the root systems, and W @ C is
    # both.  The Gram entries are integers well under 2^53, so the ranking is not even
    # approximate -- but nothing rests on it either way: exact() decides.
    Wf = W.astype(np.float64)
    CZf = Wf @ Cf
    O0 = o_of_q(np.array(N.tolist(), dtype=np.float64) / D)
    seeds = [O0]
    for _ in range(starts):
        A = rng.standard_normal((len(C), len(C)))
        Qr, _r = np.linalg.qr(A)
        seeds.append(Qr * np.sign(np.linalg.det(Qr)))
    best = (cnt, N, D)

    def _better(c2, N2, D2):
        """a strictly larger layer, or the same layer on a SMALLER denominator

        D is not part of the claim -- verify_poles.py re-checks whatever is shipped -- but a
        smaller one keeps the integer products further from the int64 ceiling, and sampling
        turns up denominators as small as 19 -- at k = 7 -- where rationalising a float
        optimum returns eight to twelve."""
        return N2 is not None and (c2 > best[0] or (c2 == best[0] and D2 < best[2]))

    # (1) The minimax objective, where the dense pass is affordable.  It is the only stage
    # that ever clears the WHOLE direction set, and it leaves margin, which is exactly what
    # the rational approximation then has to spend.
    # best[0] < n as well: where the layer is ALREADY the whole direction set there is
    # nothing above it, and at k = 14 this stage is 1932^2 per gradient step -- half an hour
    # to re-confirm 1932.  A layer short of n still gets the full sweep, and a run that
    # reaches n keeps going through the remaining seeds, since a later one may reach it on a
    # smaller denominator.
    if n <= MINIMAX_MAX and best[0] < n:
        for O in seeds:
            _g, Ob = minimax(V, O, thr)
            c2, N2, D2 = rationalise(q_of_o(Ob), Cf, Wf, CZf, m, exact, rng)
            if _better(c2, N2, D2):
                best = (c2, N2, D2)
    # (2) The count objective.  It is not redundant: where clearing everything is out of
    # reach -- every k from 17 up -- maximising the count is the right question, and there
    # it beats the minimax rotation.  Where it IS reached there is nothing above it, since
    # the layer is a copy of Z, so the descent is skipped rather than run to confirm it.
    for margin in (margins if best[0] < n else ()):
        for O in seeds:
            _b, Ob = descend(V, O, thr, m, margin * m, steps, block)
            c2, N2, D2 = rationalise(q_of_o(Ob), Cf, Wf, CZf, m, exact, rng)
            if _better(c2, N2, D2):
                best = (c2, N2, D2)
            if best[0] >= n:
                break
        if best[0] >= n:
            break
    return best


def roots_pass(write, rng, ks=tuple(range(3, 9))):
    """dimensions 75-80: the same trick on the root systems, where it helps"""
    print()
    print('=' * 92)
    print('  DIMENSIONS 75-80: the root systems')
    print('  k  dim  roots  shipped  found  gain     D      at tau(k)?')
    print('=' * 92)
    tot = 0
    for k in ks:
        R, G = root_system(k)
        m = 2
        # What is on disk: the rotated layer where one has been written, the shell layer
        # otherwise.  This used to be the SHELL count alone, so a re-run that came out below
        # the shipped ROTATED layer would have overwritten it -- 76 was the number being
        # compared against while 124 was the number shipped.  Nothing was ever lost, because
        # the search is seeded from the shipped rotation and so cannot return less than it,
        # but the guard was one seeding failure away from a silent regression.
        shipped = 2 * len(np.load(os.path.join(DATA, 'poles_k%d.npy' % k)))
        _rotf = os.path.join(DATA, 'poles_rot_k%d.npz' % k)
        shipped_D = None
        if os.path.exists(_rotf):
            _z = np.load(_rotf)
            shipped = max(shipped, len(_z['keep']))
            shipped_D = int(_z['D'][0])
        Rf = R.astype(np.float64)
        RGf = (R @ G).astype(np.float64)
        best = (-1, None, None)
        _Ns, _Ds = shipped_rotation(k)
        if _Ns is not None and 2 * m * _Ds < 2 ** 62 and preserves_span(R, _Ns):
            _A = Rf @ np.array(_Ns.tolist(), dtype=np.float64).T
            best = (int((np.abs(_A @ RGf.T).max(1) <= np.sqrt(3) / 2 * m * _Ds).sum()),
                    _Ds, _Ns)
        # The root systems are small enough that this is worth doing properly: a try is one
        # 8x8 Fraction elimination and one 240x8 product.  It is also the only stage that
        # produces a genuinely SMALL denominator -- D = 21 clears all 126 roots of E_7,
        # against ten digits for anything the rationalisation of a float optimum returns --
        # so the comparison takes the smaller D at an equal count rather than the first seen.
        for _t in range(SAMPLES_ROOT):
            _rad = 1 + (_t % 3)
            K = rng.integers(-_rad, _rad + 1, (k, k))
            K = np.triu(K, 1)
            K = (K - K.T).tolist()
            N, D = cayley(K, G)
            if N is None or 2 * m * D >= 2 ** 62:
                continue
            A = Rf @ np.array(N.tolist(), dtype=np.float64).T
            sc = int((np.abs(A @ RGf.T).max(1) <= np.sqrt(3) / 2 * m * D).sum())
            if (sc > best[0] or (sc == best[0] and best[1] is not None and D < best[1])) \
                    and preserves_span(R, N):
                best = (sc, D, N)
        _sc, D, N = best
        # Descent, from the best sampled rotation.  Sampling stalls two roots short at k = 6
        # and thirty-two short at k = 8; descent reaches all of them.
        def _exact_r(Nx, Dx, _R=R, _G=G, _m=m):
            _Ni = np.array(Nx.tolist(), dtype=np.int64)
            assert (_Ni.astype(object) == Nx).all(), 'N does not fit in int64'
            _lim = isqrt(3 * _m * _m * Dx * Dx)
            assert _lim * _lim < 3 * _m * _m * Dx * Dx < (_lim + 1) ** 2
            if not preserves_span(_R, Nx):
                return -1
            _X = (_R @ _Ni.T) @ (_R @ _G).T
            return int((2 * np.abs(_X).max(1) <= _lim).sum())
        _sc, N, D = improve(R, G, m, N, D, _sc, _exact_r, rng)
        Ni = np.array(N.tolist(), dtype=np.int64)
        lim = isqrt(3 * m * m * D * D)
        assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2
        X = (R @ Ni.T) @ (R @ G).T
        good = 2 * np.abs(X).max(1) <= lim
        cnt = int(good.sum())
        assert cnt == _sc, ('k=%d: float screen and exact test disagree' % k)
        C = G.astype(object)
        assert ((Ni.astype(object).T @ C @ Ni.astype(object)) == (D * D) * C).all(), k
        gain = max(0, cnt - shipped)
        tot += gain
        # `shipped` is the larger of the shell layer and the rotated one, so the note
        # cannot say "the shell layer is kept" any more: at k = 3 that layer is the rotated
        # 12 against a shell of 8.  Say which of the three cases this is.
        if cnt > shipped:
            _note = '   -- an improvement; --write would take it'
        elif cnt == shipped and shipped_D is not None and D < shipped_D:
            _note = '   -- the same layer on a smaller denominator (%d)' % shipped_D
        elif cnt == shipped:
            _note = '   -- the shipped layer, re-derived'
        else:
            _note = '   -- below the shipped layer, which is kept'
        print('  %2d  %3d  %5d  %7d  %5d  %+6d  %6d   %s%s'
              % (k, 72 + k, len(R), shipped, cnt, cnt - shipped, D,
                 'YES' if cnt == ROOTTAU[k] else '%d short' % (ROOTTAU[k] - cnt), _note))
        if write and (cnt > shipped
                      or (cnt == shipped and shipped_D is not None and D < shipped_D)):
            np.savez_compressed(os.path.join(DATA, 'poles_rot_k%d.npz' % k),
                                N=Ni, D=np.array([D]), keep=np.nonzero(good)[0].astype(np.int32))
    print('=' * 92)
    print('  gain over dimensions 75-80: %+d spheres' % tot)
    return tot


def main(write=False, seed=20260827, skip=False, ks=None):
    print(__doc__)
    ROWS = rows()
    rng = np.random.default_rng(seed)
    out, tot = [], 0
    # A named list of k runs just those, cross-sections and root systems alike; the default
    # is everything.  Writes are gated on an improvement either way, so this is a budget
    # knob and not a claim about which k were checked.
    _cross = [k for k in KS if ks is None or k in ks]
    _roots = [k for k in range(3, 9) if ks is None or k in ks]
    print('=' * 92)
    print('  k  dim   |W|     shipped  rotated  gain     D (digits)   at tau(k)?')
    print('=' * 92)
    for k in _cross:
        if skip and os.path.exists(os.path.join(DATA, layer_file(k, ROWS))):
            _z = np.load(os.path.join(DATA, layer_file(k, ROWS)))
            print('  %2d  %3d   already written: %d poles' % (k, 72 + k, len(_z['keep'])))
            continue
        W, c, m = directions(k, ROWS)
        n = W.shape[1]
        # THE SPAN DECIDES WHETHER THIS FAMILY CAN WORK AT ALL.  Everything below builds a
        # Cayley transform on the n coordinates W is written in; that is an isometry of the
        # AXIS SPACE only when W has full rank in them.  Where it does not, preserves_span
        # rejects every candidate, and the right thing is to say so once rather than to
        # discover it a candidate at a time -- which is what happened between 2026-08-27,
        # when the shipped rotation at those k became a valid one from axis_block.py, and
        # today: the "no isometry" branch below stopped firing, and the file ran a full
        # descent against a family that cannot beat what is already there.
        if len(span_basis(W)) < n:
            print('  %2d  %3d   |Z| spans %d of %d coordinates: no isometry of R^%d in this '
                  'family,\n              see scripts/axis_block.py'
                  % (k, 72 + k, len(span_basis(W)), n, k))
            continue
        # blocked: the full Gram is 18.5 GiB at k = 22
        CZi = (W * c).T
        _mx = -m
        for _s in range(0, len(W), max(1, int(8e6 // max(1, len(W))))):
            _B = W[_s:_s + max(1, int(8e6 // max(1, len(W))))] @ CZi
            for _t in range(len(_B)):
                _B[_t, _s + _t] = -m
            _mx = max(_mx, int(_B.max()))
        assert 2 * _mx <= m, 'W is not a 60-degree code'
        t0 = time.time()
        # more tries buy only a handful of points; the cost is |W|^2 per try, which is
        # 8.7e9 integer products at k = 23, so the budget shrinks as |W| grows
        tries = 60 if len(W) < 6000 else (20 if len(W) < 20000 else 3)
        Wf = W.astype(np.float64)
        CZf = (W * c).T.astype(np.float64)
        best = (-1, None, None)
        _Ns, _Ds = shipped_rotation(k, layer_file(k, ROWS))
        if _Ns is not None and 2 * m * _Ds < 2 ** 62 and preserves_span(W, _Ns):
            best = (screen(Wf, CZf, m, _Ns, _Ds), _Ds, _Ns)
        for Nx, Dx in borrowed(k, W):
            best = max(best, (screen(Wf, CZf, m, Nx, Dx), Dx, Nx), key=lambda t: t[0])
        for _ in range(tries):
            K = rng.integers(-1, 2, (n, n))
            K = np.triu(K, 1)
            K = (K - K.T).tolist()
            N, D = cayley(K, np.diag(c))
            if N is None or 2 * m * D >= 2 ** 62:
                continue
            sc = screen(Wf, CZf, m, N, D)
            # A rotation that leaves the span is not a rotation of the axis space, however
            # many caps it clears.  Seven layers were written before this was asked.
            if (sc > best[0] or (sc == best[0] and best[1] is not None and D < best[1])) \
                    and preserves_span(W, N):
                best = (sc, D, N)
        _sc, D, N = best
        if N is None:
            # Every candidate left the axis space, which is what happens exactly when the
            # direction set is written in more coordinates than it spans.  Those k are
            # derived by scripts/axis_block.py, whose generator is built from the direction
            # set itself and fixes the orthogonal complement pointwise.
            print('  %2d  %3d   no isometry of R^%d in this family: see scripts/axis_block.py'
                  % (k, 72 + k, k))
            continue
        # Descent, from the best sampled rotation.  It cannot help where the ambient
        # coordinates are larger than the span -- the Cayley rationalisation is an ambient
        # isometry there and preserves_span rejects it -- and those k are derived by
        # scripts/axis_block.py instead.
        def _exact_x(Nx, Dx, _W=W, _c=c, _m=m):
            if not preserves_span(_W, Nx):
                return -1
            _g, _l, _w = cleared(_W, _c, _m, Nx, Dx)
            return int(_g.sum())
        _sc, N, D = improve(W, np.diag(c), m, N, D, _sc, _exact_x, rng)
        good, lim, worst = cleared(W, c, m, N, D)      # the exact pass, once
        cnt = int(good.sum())
        assert cnt == _sc, ('k=%d: the float screen and the exact test disagree (%d vs %d)'
                            % (k, _sc, cnt))
        assert check_orthogonal(N, D, c), 'k=%d: N^T C N != D^2 C' % k
        # rows81-95.json is written by final81-95.py FROM these files, so reading the
        # layer back out of it makes this guard depend on which driver ran last.  Read the
        # file the layer actually lives in.
        was = ROWS[k]['poles']
        _rotf = os.path.join(DATA, layer_file(k, ROWS))
        # rows81-95.json's count belongs to the layer that row was built from, which after a
        # cap-set change was a DIFFERENT W: comparing against it would ask a 756-vector set
        # to beat a 1106-pole layer built on another one, and refuse to write anything.
        if ROWS[k]['W'] == 'kappa' and not os.path.exists(_rotf):
            was = 0
        was_D = None
        if os.path.exists(_rotf):
            _z = np.load(_rotf)
            was = max(was, len(_z['keep']))
            was_D = int(_z['D'][0])
        tot += cnt - was
        out.append((k, cnt, N, D, good))
        # only on an improvement: scripts/axis_block.py writes the same files from a
        # different family, and neither dominates the other
        if write and (cnt > was or (cnt == was and was_D is not None and D < was_D)):
            # The digest of W travels with the layer: indices into a cap set and a
            # rotation of it are both meaningless against a different one, and since
            # 2026-08-30 two cap sets exist at k = 12 and 13.  scripts/capsha.py.
            np.savez_compressed(os.path.join(DATA, layer_file(k, ROWS)),
                                N=np.array(N.tolist(), dtype=np.int64), D=np.array([D]),
                                keep=np.nonzero(good)[0].astype(np.int32),
                                wsha=capsha.digest(W, c))
        print('  %2d  %3d  %6d   %6d   %6d   %+7d   %5d        %s   (%.0fs)'
              % (k, 72 + k, len(W), was, cnt, cnt - was, len(str(D)),
                 'YES' if cnt == TAU[k] else '%d short' % (TAU[k] - cnt), time.time() - t0))
    print('=' * 92)
    print('  total gain over dimensions 81-95: %+d spheres' % tot)
    print('  (dimension 83 is not here: its k = 11 configuration lives over Z[sqrt2])')
    if _roots:
        roots_pass(write, rng, _roots)
    if write:
        print('  each layer was written to data/poles_rot_<k>.npz as it was found')
    return 0


if __name__ == '__main__':
    _ks = tuple(int(x) for x in sys.argv[1:] if x.isdigit()) or None
    sys.exit(main(write='--write' in sys.argv, skip='--skip-existing' in sys.argv, ks=_ks))
