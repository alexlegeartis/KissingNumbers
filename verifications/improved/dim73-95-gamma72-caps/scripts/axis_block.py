#!/usr/bin/env python3
"""Axis layers whose rotation is an isometry of the AXIS SPACE, not of the ambient coordinates.

    python axis_block.py            derive and compare with what is shipped
    python axis_block.py --write    derive and overwrite data/poles_rot_<k>.npz

WHY THIS FILE EXISTS.  The first version of the rotated axis layers (axis_rotate.py,
2026-08-27) built the rotation as a Cayley transform of an integer skew matrix on the
coordinates the direction set is written in.  For most k that is the axis space and the
construction is sound.  For k = 13, 15, 17, 18, 21, 22 and 23 it is not: those direction sets
are written in more coordinates than they span -- Lambda_21's cross-section is 27 720 vectors
written in Lambda_24's 24 coordinates and spanning 21 of them -- and an isometry of R^24 is
not an isometry of R^21.  It carries the poles out of the axis space.  Every inner product
those layers were checked on was correct; the configurations needed R^96 instead of R^93.

Exact ranks, before the fix:

    k                 13  15  17  18  21  22  23
    Z spans           13  15  17  18  21  22  23
    Z and QZ span     14  16  24  24  24  24  24

WHY NOT JUST WORK IN A BASIS OF THE SPAN.  Because of the denominator.  Write W = U P with P
a lattice basis of the span and carry the metric across as G = P C P^T; then every isometry
of (Z^k, G) is an isometry of the axis space, which is right.  But S = G^{-1} K has
denominators dividing det(G), and det(G) is 2^59 at k = 17 and 2^71 at k = 23.  det(G) is a
lattice invariant, so no choice of basis avoids it, and the exact test overflows int64.

WHAT IS DONE INSTEAD.  Build the rotation out of the direction set.  For v and w rows of W,

    S = v (Cw)^T - w (Cv)^T

is an integer matrix; CS = (Cv)(Cw)^T - (Cw)(Cv)^T is skew, which is exactly the condition
for (I - S)(I + S)^{-1} to preserve the form.  S sends v to -m w and w to m v and annihilates
everything C-orthogonal to both, so on the plane of v and w it is m times the rotation
generator and it is zero elsewhere.  Taking r pairs, all 2r vectors mutually C-orthogonal and
each of norm m, and scaling by 1/lambda, the Cayley transform is a direct sum of plane
rotations through 2 arctan(m/lambda) and can be written down:

    N = D I + sum_t [ -2m (v_t (Cv_t)^T + w_t (Cw_t)^T)
                      + 2 lambda (v_t (Cw_t)^T - w_t (Cv_t)^T) ],      D = lambda^2 + m^2

with N^T C N = D^2 C.  No inversion, no rationals, and D between 80 and 5120 in practice
against the 2^71 the basis route forces.  Every v_t and w_t lies in the span, so Q fixes its
C-orthogonal complement pointwise and the poles cannot leave R^k -- and that is checked here
rather than assumed, on a basis of the span, which settles it because Q is linear.

WHAT IS HEURISTIC AND WHAT IS NOT.  The choice of frame and of lambda is a search and can be
beaten.  Nothing rests on it: scripts/verify_poles.py takes the shipped N, D and index set
and re-checks every cap pair in exact integer arithmetic, the isometry, and the span.
"""
import os
import sys
import time
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
sys.path.insert(0, HERE)
from axis_rotate import directions, rows                              # noqa: E402

# the direction sets written in more coordinates than they span
AMBIENT_ABOVE_SPAN = (13, 15, 17, 18, 21, 22, 23)
LAMS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 16, 19, 21, 24, 27, 31, 32, 37, 43, 53, 64)


def frame(W, c, rng, want):
    """as many mutually C-orthogonal rows of W as a greedy pass finds"""
    Cz = W * c
    chosen = []
    for i in rng.permutation(len(W)):
        if not chosen or not (W[chosen] @ Cz[i]).any():
            chosen.append(int(i))
            if len(chosen) >= want:
                break
    return W[chosen]


def block_rotation(Fm, c, m, lam):
    """(N, D): plane rotations through 2 arctan(m/lambda), in closed form"""
    amb = Fm.shape[1]
    D = lam * lam + m * m
    co = np.array(c, dtype=object)
    N = D * np.eye(amb, dtype=object)
    for t in range(len(Fm) // 2):
        v, w = Fm[2 * t].astype(object), Fm[2 * t + 1].astype(object)
        cv, cw = v * co, w * co
        N = (N - 2 * m * (np.outer(v, cv) + np.outer(w, cw))
             + 2 * lam * (np.outer(v, cw) - np.outer(w, cv)))
    return N, D


def is_isometry(N, D, c):
    C = np.diag(np.array(c, dtype=object))
    return bool(((N.T @ C @ N) == (D * D) * C).all())


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
    """Does N map the span of the directions into itself?  Q is linear, so a basis settles it."""
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


def exact_keep(W, c, m, N, D):
    """boolean per rotated direction, in exact integer arithmetic"""
    Ni = np.array(N.tolist(), dtype=np.int64)
    assert (Ni.astype(object) == N).all(), 'N does not fit in int64'
    lim = isqrt(3 * m * m * D * D)
    assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2, '3 m^2 D^2 is a perfect square'
    A = W @ Ni.T
    CZ = (W * c).T
    n = len(W)
    B = max(1, int(8e6 // max(1, n)))
    good = np.empty(n, bool)
    worst = 0
    for s in range(0, n, B):
        X = A[s:s + B] @ CZ
        mx = np.abs(X).max(1)
        worst = max(worst, int(mx.max()))
        good[s:s + B] = 2 * mx <= lim
    assert 2 * worst < 2 ** 63 - 1, 'the int64 product could overflow'
    return int(good.sum()), good


def screen(Wf, CZf, m, N, D, rows_=None):
    """float ranking only: nothing is decided here, the winner is re-tested in integers"""
    Nf = np.array(N.tolist(), dtype=np.float64)
    X = (Wf if rows_ is None else Wf[rows_])
    n = len(X)
    B = max(1, int(8e6 // max(1, len(Wf))))
    thr = np.sqrt(3.0) / 2.0 * m * D
    A = X @ Nf.T
    return sum(int((np.abs(A[s:s + B] @ CZf).max(1) <= thr).sum()) for s in range(0, n, B))


def derive(k, rng, trials=None, sub=3000, top=3):
    """(count, N, D, keep) for the best block rotation found at this k"""
    W, c, m = directions(k, rows())
    n, amb = len(W), W.shape[1]
    trials = trials or (30 if n < 3000 else (10 if n < 20000 else 5))
    Wf = W.astype(np.float64)
    CZf = (W * c).T.astype(np.float64)
    sel = None if n <= sub else rng.choice(n, size=sub, replace=False)
    cands = []
    for _ in range(trials):
        Fm = frame(W, c, rng, 2 * (amb // 2))
        if len(Fm) < 2:
            continue
        for lam in LAMS:
            N, D = block_rotation(Fm, c, m, lam)
            if 2 * m * D >= 2 ** 62:
                continue
            cands.append((screen(Wf, CZf, m, N, D, sel), N, D))
    if not cands:
        return None
    cands.sort(key=lambda t: -t[0])
    best = (-1, None, None)
    for _s, N, D in (cands[:top] if sel is not None else cands[:1]):
        sc = screen(Wf, CZf, m, N, D)
        if sc > best[0]:
            best = (sc, N, D)
    sc, N, D = best
    cnt, good = exact_keep(W, c, m, N, D)
    assert cnt == sc, 'k=%d: the float screen and the exact test disagree (%d vs %d)' % (k, sc, cnt)
    assert is_isometry(N, D, c), 'k=%d: N^T C N is not D^2 C' % k
    assert preserves_span(W, N), 'k=%d: the rotation does not preserve the axis space' % k
    return cnt, N, D, np.nonzero(good)[0].astype(np.int32)


def main(write=False, seed=20260827, ks=AMBIENT_ABOVE_SPAN):
    print(__doc__)
    rng = np.random.default_rng(seed)
    print('=' * 92)
    print('  k  dim   |Z|   ambient  span   shipped   block   D       in R^k?  kept?')
    print('=' * 92)
    for k in ks:
        W, c, m = directions(k, rows())
        f = os.path.join(DATA, 'poles_rot_%d.npz' % k)
        shipped = len(np.load(f, allow_pickle=True)['keep']) if os.path.exists(f) else 0
        t0 = time.time()
        got = derive(k, rng)
        if got is None:
            print('  %2d  %3d  NO USABLE ROTATION' % (k, 72 + k))
            continue
        cnt, N, D, keep = got
        print('  %2d  %3d %6d %8d %5d %9d %7d  %-7d YES     %-7s %.0fs'
              % (k, 72 + k, len(W), W.shape[1], len(span_basis(W)), shipped, cnt, D,
                 'yes' if cnt > shipped else 'no', time.time() - t0))
        sys.stdout.flush()
        # Write only if it beats what is there.  Three families produce these layers -- the
        # lattice shell, the ambient Cayley of axis_rotate.py, and the blocks here -- and
        # none dominates: the blocks win at k = 12 and lose at k = 10.  Writing only on an
        # improvement lets them run in any order.
        if write and cnt > shipped:
            # int64, not object: verify_poles.py loads these WITHOUT allow_pickle, and an
            # object array cannot be read that way.  N fits by construction -- exact_keep
            # has already asserted it -- and the assert here says so rather than assuming.
            Ni = np.array(N.tolist(), dtype=np.int64)
            assert (Ni.astype(object) == N).all(), 'k=%d: N does not fit in int64' % k
            np.savez_compressed(f, N=Ni, D=np.array([D]), keep=keep)
    print('=' * 92)
    if write:
        print('  written to data/poles_rot_<k>.npz; scripts/verify_poles.py checks them')
    else:
        print('  nothing written: pass --write to overwrite the shipped layers')


if __name__ == '__main__':
    # --seed lets the search be re-run for a different frame sequence.  Because a layer is
    # written only when it beats the one on disk, repeated runs improve the data and can
    # never damage it, so this is the intended way to spend more effort on it.
    _a = [x for x in sys.argv[1:] if x.isdigit()]
    _seed = 20260827
    if '--seed' in sys.argv:
        _seed = int(sys.argv[sys.argv.index('--seed') + 1])
        _a = [x for x in _a if int(x) != _seed]
    _ks = tuple(int(x) for x in _a) or AMBIENT_ABOVE_SPAN
    main('--write' in sys.argv, seed=_seed, ks=_ks)
