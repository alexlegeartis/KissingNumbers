#!/usr/bin/env python3
"""The axis layers of dimensions 84 and 85, as rotated copies of the Kappa cap sets.

    python kappa_poles.py             re-derive and report against the shipped layers
    python kappa_poles.py --write     re-derive and write data/poles_kap_<k>.npz
    python kappa_poles.py --write 12  just that k

WHY NOT axis_rotate.py.  That file builds its rotation as a Cayley transform of a random
integer skew matrix and then descends on SO(k) with a minimax objective, which is the right
tool when nothing is known about the direction set.  Here something is: K_12 is an EISENSTEIN
lattice, its coordinates come in six blocks with metric diag(1, 3), and each block is a copy
of Z[omega] written as a + b sqrt(-3).  A rotation that respects that structure is a product
of one plane rotation per block and has denominator D alone -- against the k^(k/2) of a full
Cayley transform, which is 1.4e6 at k = 12 before the descent even starts.

THE PLANE ROTATIONS.  In a block with metric weights (g_i, g_j),

    R = [[a, -(g_j/g_i) b], [b, a]] / D        with   a^2 + (g_j/g_i) b^2 = D^2

is rational and an isometry of that block.  For (1, 3) that is a^2 + 3 b^2 = D^2 outright; for
(3, 1) the same conic reappears after b = 3c, as R = [[a, -c], [3c, a]] / D.  So both cap sets
are driven by ONE conic, whose rational points are dense in the circle.

WHAT THE IDENTITY BLOCK MEANS.  (a, b) = (D, 0) leaves a block alone, and taking it in every
block is the identity rotation, under which every pole coincides with a cap and nothing
clears.  Taking a = 1, b = 1, D = 2 in every block is multiplication by the unit
(1 + sqrt(-3))/2 of Z[omega], which is an automorphism of K_12 -- the poles are the caps
again, in a different order, and again nothing clears.  The layer comes from choosing the
blocks DIFFERENTLY from each other, which is exactly what breaks the Eisenstein structure the
cap set has.

THE SEARCH is a greedy over blocks with a plateau escape: each round tries every (block,
angle) against the exact count, takes the best strict improvement, and when there is none
takes the best sideways move instead, keeping the best PREFIX of the accepted sequence rather
than its end.  A round that improves nothing and has no sideways move left stops it.

WHAT IS HEURISTIC AND WHAT IS NOT.  The search is a heuristic and may come out differently on
a different machine.  Nothing rests on it: the layer is written with the digest of the cap set
it was built from (scripts/capsha.py), and scripts/verify_poles.py re-checks the shipped N, D
and index set against every cap in exact integer arithmetic.
"""
import os
import sys
import time
from math import gcd, isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
sys.path.insert(0, HERE)
import capsha                                                          # noqa: E402

KS = (12, 13)


def capset(k):
    """(W, c, m): the Kappa cap directions, their metric, and their common norm"""
    W = np.load(os.path.join(DATA, 'kap%d_W.npy' % k)).astype(np.int64)
    c = np.load(os.path.join(DATA, 'kap%d_c.npy' % k)).astype(np.int64)
    m = int((W[0] * c) @ W[0])
    assert (np.einsum('ij,ij->i', W * c, W) == m).all(), 'the directions are not of one norm'
    return W, c, m


def conic(dmax):
    """{D: [(a, b)]} with a^2 + 3 b^2 = D^2, b > 0, gcd(a, b, D) = 1, plus the identity"""
    out = {}
    for D in range(2, dmax + 1):
        sols = []
        for b in range(1, D + 1):
            r = D * D - 3 * b * b
            if r < 0:
                break
            a = isqrt(r)
            if a * a == r and gcd(gcd(a, b), D) == 1:
                sols.append((a, b))
        if sols:
            out[D] = sols
    return out


def blocks(c):
    """the coordinate pairs a plane rotation can act on, as (i, j, weights)"""
    out = []
    for t in range(0, len(c) - 1, 2):
        gi, gj = int(c[t]), int(c[t + 1])
        if {gi, gj} == {1, 3}:
            out.append((t, t + 1, gi, gj))
    return out


def rotation(c, choice, D):
    """the block-diagonal integer N with N^T C N = D^2 C"""
    n = len(c)
    N = np.eye(n, dtype=object) * D
    for (i, j, gi, gj), (a, b) in zip(blocks(c), choice):
        if b == 0:
            continue
        if (gi, gj) == (1, 3):
            N[i, i], N[i, j], N[j, i], N[j, j] = a, -3 * b, b, a
        else:                                    # (3, 1): the same conic after b = 3c
            N[i, i], N[i, j], N[j, i], N[j, j] = a, -b, 3 * b, a
    return N


def check_isometry(N, c, D):
    C = np.zeros((len(c), len(c)), dtype=object)
    for t in range(len(c)):
        C[t, t] = int(c[t])
    return bool(((N.T @ C @ N) == (D * D) * C).all())


def cleared(W, c, m, N, D):
    """boolean per direction, exactly: 4 <Qw, z>^2 <= 3 m^2, with Q = N/D"""
    lim = isqrt(3 * m * m * D * D)
    assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2, '3 m^2 D^2 is a perfect square'
    A = W @ np.array(N.tolist(), dtype=np.int64).T
    CZ = (W * c).T
    assert 2 * int(np.abs(A).max()) * int(np.abs(CZ).max()) * W.shape[1] < 2 ** 62, \
        'the integer product would overflow'
    good = np.empty(len(W), bool)
    B = max(1, int(8e6 // max(1, len(W))))
    for s in range(0, len(W), B):
        X = A[s:s + B] @ CZ
        good[s:s + B] = 2 * np.abs(X).max(1) <= lim
    return good


def search(k, W, c, m, dmax=30, rounds=40, budget=600.0):
    """(count, N, D, good): a block rotation clearing as many caps as possible"""
    bl = blocks(c)
    cs = conic(dmax)
    t0 = time.time()
    best = (-1, None, None, None)
    for D in sorted(cs):
        angles = [(D, 0)] + cs[D]
        choice = [(D, 0)] * len(bl)
        cur = -1
        seen = []
        stall = 0
        for _ in range(rounds):
            if time.time() - t0 > budget:
                break
            gain, alt = None, []
            for bi in range(len(bl)):
                for ang in angles:
                    if choice[bi] == ang:
                        continue
                    trial = list(choice)
                    trial[bi] = ang
                    N = rotation(c, trial, D)
                    if not check_isometry(N, c, D):
                        continue
                    n = int(cleared(W, c, m, N, D).sum())
                    if n > cur:
                        if gain is None or n > gain[0]:
                            gain = (n, bi, ang)
                    else:
                        alt.append((n, bi, ang))
            if gain is not None:
                cur, bi, ang = gain
                choice[bi] = ang
                stall = 0
            elif alt and stall < 3:
                cur, bi, ang = max(alt, key=lambda t: t[0])
                choice[bi] = ang
                stall += 1
            else:
                break
            seen.append((cur, list(choice)))
            if cur >= len(W):
                break
        if seen:
            b = max(seen, key=lambda t: t[0])
            if b[0] > best[0]:
                N = rotation(c, b[1], D)
                assert check_isometry(N, c, D), 'the winning rotation is not an isometry'
                good = cleared(W, c, m, N, D)
                assert int(good.sum()) == b[0], 'the score and the exact test disagree'
                best = (b[0], N, D, good)
                print('      D = %-4d %6d of %d' % (D, b[0], len(W)))
        if best[0] >= len(W):
            break
    return best


def main(write, ks):
    print(__doc__)
    print('   k  dim    |W|   shipped   rotated   D    at |W|?')
    print('   ' + '-' * 62)
    for k in ks:
        W, c, m = capset(k)
        f = os.path.join(DATA, 'poles_kap_%d.npz' % k)
        was = 0
        if os.path.exists(f):
            z = np.load(f)
            if capsha.matches(z, W, c):
                was = len(z['keep'])
        cnt, N, D, good = search(k, W, c, m)
        if N is None:
            print('   %2d  %3d %6d   %7d   no rotation found' % (k, 72 + k, len(W), was))
            continue
        print('   %2d  %3d %6d   %7d   %7d %4d   %s'
              % (k, 72 + k, len(W), was, cnt, D, 'YES' if cnt == len(W) else
                 '%d short' % (len(W) - cnt)))
        if write and cnt > was:
            np.savez_compressed(f, N=np.array(N.tolist(), dtype=np.int64),
                                D=np.array([D]),
                                keep=np.nonzero(good)[0].astype(np.int32),
                                wsha=capsha.digest(W, c))
            print('        wrote data/poles_kap_%d.npz (%d poles, D = %d)' % (k, cnt, D))
    return 0


if __name__ == '__main__':
    _ks = [int(a) for a in sys.argv[1:] if a.isdigit()] or list(KS)
    sys.exit(main('--write' in sys.argv, _ks))
