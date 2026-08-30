#!/usr/bin/env python3
"""The axis layer of dimension 83, over Z[sqrt 2].

    python axis_sqrt2.py             re-derive and report
    python axis_sqrt2.py --write     re-derive and write data/poles_sqrt2_11.npz

DIMENSION 83 HAS NO AXIS LAYER, and the reason on record is arithmetic rather than
geometric: its 604 cap directions are the record configuration of R^11, whose coordinates
are (A + B sqrt2)/3, so the Gram lands in Z[sqrt 2] and the rational machinery of
scripts/axis_rotate.py has nothing to compare with.  The geometry is no different from any
other k -- a pole needs 30 degrees from the caps and 60 degrees from the other poles -- and
the caps are just as small here as everywhere else.

THE ROTATION STAYS RATIONAL.  Only the CONFIGURATION lives in Z[sqrt 2]; the rotation need
not.  Take Q = N/D with N integral and N^T C N = D^2 C for the same diagonal metric C the
cap file carries.  Then Q x has coordinates in the same ring, the pole-pole condition is
inherited from the 604 points being a 60-degree code (scripts/verify11.py checks that), and
the only thing left is the cap condition -- one comparison in Z[sqrt 2] per pole and cap.

THE COMPARISON.  With 9 D <x_j, N x_i> = P + Q sqrt2 for integers P, Q, and both norms equal
to 4, the condition 4 <z,a>^2 <= 3 |z|^2 |a|^2 is

    (P + Q sqrt2)^2 <= 972 D^2,   i.e.   (972 D^2 - P^2 - 2 Q^2) - 2 P Q sqrt2 >= 0,

and a + b sqrt2 >= 0 is decided by the signs of a and b and, when they disagree, by a^2
against 2 b^2.  Integers only, and Python integers at that, because D^2 has fifteen digits.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from axis_rotate import _euclid, minimax, rationalise               # noqa: E402

PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
K = 11
TAU11 = 604


def load():
    """(A, B, c): the 604 cap directions as (A + B sqrt2)/3, and the diagonal metric"""
    z = np.load(os.path.join(DATA, 'rec11_exact.npz'))
    A, B, c = z['A'].astype(np.int64), z['B'].astype(np.int64), z['c'].astype(np.int64)
    P = (A * c) @ A.T + 2 * ((B * c) @ B.T)
    Q = (A * c) @ B.T + (B * c) @ A.T
    N = int(P[0, 0])
    assert N == 36 and (np.diag(P) == N).all() and (np.diag(Q) == 0).all(), \
        'the cap directions are not of one rational norm'
    return A, B, c


def nonneg(a, b):
    """a + b sqrt2 >= 0, elementwise, for object arrays of Python integers"""
    a2, b2 = a * a, b * b
    both_up = (a >= 0) & (b >= 0)
    both_dn = (a < 0) & (b <= 0)
    mixed_a = (a >= 0) & (b < 0)                    # a >= -b sqrt2  <=>  a^2 >= 2 b^2
    mixed_b = (a < 0) & (b > 0)                     # b sqrt2 >= -a  <=>  2 b^2 >= a^2
    out = np.zeros(a.shape, dtype=bool)
    out |= both_up
    out |= mixed_a & (a2 >= 2 * b2)
    out |= mixed_b & (2 * b2 >= a2)
    assert not (both_dn & out).any()
    return out


def cleared(A, B, c, N, D):
    """boolean per pole: is (N/D) x_i at least 30 degrees from every cap?"""
    No = np.array(N.tolist(), dtype=object)
    Ai, Bi = A.astype(object) @ No.T, B.astype(object) @ No.T
    Ac, Bc = (A * c).astype(object), (B * c).astype(object)
    P = Ac @ Ai.T + 2 * (Bc @ Bi.T)                 # P[j, i], the cap j against the pole i
    Q = Ac @ Bi.T + Bc @ Ai.T
    D = int(D)
    lim = 972 * D * D
    good = nonneg(lim - P * P - 2 * Q * Q, -2 * P * Q)
    return good.all(0)                              # a pole must clear EVERY cap


def main():
    write = '--write' in sys.argv
    A, B, c = load()
    n = len(A)
    C = np.diag(c).astype(np.int64)
    X = (A.astype(np.float64) + np.sqrt(2.0) * B.astype(np.float64)) / 3.0
    m = 4
    assert abs(float((X[0] * c) @ X[0]) - m) < 1e-9
    print(__doc__)
    print('   %d cap directions in R^11, metric diag%s, norm %d'
          % (n, tuple(int(v) for v in sorted(set(c.tolist())))[:4], m))

    def exact(Nn, Dd):
        return int(cleared(A, B, c, Nn, Dd).sum())

    L, q_of_o, o_of_q = _euclid(C)
    V = X @ L
    thr = np.sqrt(3.0) / 2.0 * m
    CZf = X @ C.astype(np.float64)
    rng = np.random.default_rng(0)
    best = (0, None, None)
    for t in range(12):
        O = np.linalg.qr(rng.standard_normal((K, K)))[0]
        O = O * np.sign(np.linalg.det(O))
        t0 = time.time()
        w, O = minimax(V, O, thr, steps=4000, eta=0.03, beta=200.0)
        cnt, Nn, Dd = rationalise(q_of_o(O), C, X, CZf, m, exact, rng=rng, top=4, keep=4)
        print('   start %2d: worst cosine %.6f (%.4f of the threshold), layer %4d  (%.0f s)'
              % (t, w / m, w / thr, cnt, time.time() - t0))
        if cnt > best[0] or (cnt == best[0] and cnt and Dd < best[2]):
            best = (cnt, Nn, Dd)
        if best[0] == n:
            break
    cnt, Nn, Dd = best
    print()
    print('   dimension 83: %d poles of a possible %d   (D = %s)%s'
          % (cnt, TAU11, Dd, '   AT THE CEILING' if cnt == TAU11 else ''))
    if write and cnt:
        keep = np.nonzero(cleared(A, B, c, Nn, Dd))[0].astype(np.int32)
        assert len(keep) == cnt
        np.savez(os.path.join(DATA, 'poles_sqrt2_11.npz'),
                 N=np.array(Nn.tolist(), dtype=np.int64),
                 D=np.array([Dd], dtype=np.int64), keep=keep)
        print('   wrote data/poles_sqrt2_11.npz')
    return 0


if __name__ == '__main__':
    sys.exit(main())
