#!/usr/bin/env python3
"""Derive the axis-layer rotation, instead of shipping it as a matrix that fell from the sky.

    python axis_rotation.py            re-derive and compare with the shipped rotation
    python axis_rotation.py --write    re-derive and overwrite data/axis_rotation*.txt

WHAT THE ROTATION IS FOR.  The axis layer of the dimension-38 construction is a set of unit
directions in R^14 that must (a) be a 60-degree code among themselves and (b) sit at least 30
degrees away from every one of the 1932 cap directions.  Taking a ROTATED copy Q.P of the cap
direction set itself makes (a) automatic, because Q is orthogonal, and leaves (b) as the only
question: how many of the 1932 rotated directions land outside every 30-degree cap.

WHY THIS FILE EXISTS.  The first version of this package took Q by sampling and kept 1908,
losing 24 points, and the README argued that the last 24 were out of reach:

    "a random rotation loses about 24 points and best-of-N would recover only a handful"

The first half is right and the second is wrong.  Sampling really is hopeless -- 400 random
rotations here lose 56 on average with a standard deviation of 11, and the caps cover about
1.24% of S^13, so an independent-points model gives exp(-24) for a clean rotation.  But SO(14)
is 91-dimensional and the objective is differentiable once the maximum over cap directions is
softened, and DESCENT walks from 1908 to 1932 in a couple of hundred steps.  The distinction
matters: "sampling cannot find it" is a measurement, "it is not there" is not, and the second
does not follow from the first.

    S[i,j] = <Q w_i, w_j>,   m_i = softmax_beta_j S[i,j],   f(Q) = sum_i softplus(m_i - thr)

with dS[i,j]/dQ = w_j w_i^T, so the gradient is the 14x14 matrix P^T C^T P and the step
Q <- expm(-eta skew(G Q^T)) Q stays exactly on the group.

FROM REAL TO RATIONAL.  The shipped check is integer arithmetic, so Q has to be N/D with N
integer and N N^T = D^2 I.  Cayley gives that for free: for any RATIONAL skew-symmetric S,
(I - S)(I + S)^{-1} is rational and orthogonal.  So the real solution is turned into
S = (I - Q)(I + Q)^{-1}, rounded to (1/d)Z, and rebuilt exactly in Fractions.  Rounding
perturbs every inner product, which is why the descent is run against a MARGIN of 0.15
rather than against the threshold itself -- a solution sitting on the boundary would not
survive being rounded.  d is kept small because D grows like d^14 and the verification wants
|<N w_i, w_j>| <= 8D to stay inside int64.

WHAT IS HEURISTIC AND WHAT IS NOT.  Everything in this file is a search: it can fail, and a
different BLAS or a different numpy may walk a different path.  Nothing rests on it.  What
the claim rests on is verify.py, which takes the shipped N and D and checks all 1932 x 1932
inner products in exact integer arithmetic.  This file explains where they came from.
"""
import os
import sys
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
n = 14
THR = 4.0 * np.sqrt(3.0)          # 8 cos(30 degrees): the cap-axis threshold, times |w||a|/8

P = np.loadtxt(os.path.join(DATA, "dim14_directions.txt"), delimiter=",", dtype=np.int64)
Pf = P.astype(np.float64)
assert P.shape == (1932, n) and ((P * P).sum(1) == 8).all()


def kept(Q, margin=0.0):
    """(how many rotated directions clear every cap, the smallest clearance)"""
    M = np.abs((Pf @ Q.T) @ Pf.T).max(1)
    return int((M <= THR - margin).sum()), float((THR - M).min())


def descend(Q, iters, margin, eta=2e-4, beta=6.0, alpha=40.0, verbose=True):
    from scipy.linalg import expm
    best, bq = kept(Q, margin)[0], Q.copy()
    for t in range(iters):
        S = (Pf @ Q.T) @ Pf.T
        mx = S.max(1, keepdims=True)
        E = np.exp(beta * (S - mx))
        Z = E.sum(1, keepdims=True)
        m = mx[:, 0] + np.log(Z[:, 0]) / beta
        u = alpha * (m - THR + margin)
        sig = 1.0 / (1.0 + np.exp(-np.clip(u, -60, 60)))
        G = Pf.T @ ((E / Z) * sig[:, None]).T @ Pf
        Om = G @ Q.T
        Q = expm(-eta * 0.5 * (Om - Om.T)) @ Q
        U, _, Vt = np.linalg.svd(Q)                      # keep it exactly orthogonal
        Q = U @ Vt
        k = kept(Q, margin)[0]
        if k > best:
            best, bq = k, Q.copy()
        if verbose and t % 100 == 0:
            print("      iteration %4d   clearing the margin: %4d of 1932" % (t, k))
    return best, bq


def cayley_exact(M, d):
    """(I - M/d)(I + M/d)^{-1} in exact rationals, M an integer skew-symmetric matrix"""
    A = [[F(-M[i][j], d) + (1 if i == j else 0) for j in range(n)] for i in range(n)]
    aug = [[F(M[i][j], d) + (1 if i == j else 0) for j in range(n)]
           + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if aug[r][c])
        aug[c], aug[p] = aug[p], aug[c]
        pv = aug[c][c]
        aug[c] = [x / pv for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c]:
                f = aug[r][c]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    Binv = [row[n:] for row in aug]
    return [[sum((A[i][k] * Binv[k][j] for k in range(n)), F(0)) for j in range(n)]
            for i in range(n)]


def clear_denominators(X):
    D = 1
    for row in X:
        for x in row:
            D = D * x.denominator // gcd(D, x.denominator)
    return [[int(x * D) for x in row] for row in X], D


def exact_report(N, D):
    """the same test verify.py runs, so this file can say whether its own output passes"""
    Ni = np.array(N, dtype=np.int64)
    T = isqrt(48 * D * D)
    IP = (P @ Ni.T) @ P.T
    return int((np.abs(IP) > T).sum()), int(np.abs(IP).max()), T


def main(write=False, seed=20260827, iters=600, margin=0.15, dmax=21):
    print(__doc__)
    print("=" * 78)
    Nsh = np.loadtxt(os.path.join(DATA, "axis_rotation.txt"), dtype=np.int64)
    Dsh = int(open(os.path.join(DATA, "axis_rotation_D.txt")).read().strip())
    k, sl = kept(Nsh.astype(np.float64) / Dsh)
    print("the SHIPPED rotation clears %d of 1932 cap directions (slack %.6f)" % (k, sl))
    print()

    rng = np.random.default_rng(seed)
    print("[1] sampling, for the record")
    ss = []
    for _ in range(200):
        A = rng.standard_normal((n, n))
        Q, R = np.linalg.qr(A)
        ss.append(kept(Q * np.sign(np.diag(R)))[0])
    ss = np.array(ss)
    print("      200 random rotations: mean %.1f, sd %.1f, best %d -- so sampling loses about"
          % (ss.mean(), ss.std(), ss.max()))
    print("      %d and the best of 200 loses %d.  This is why the first version kept 1908."
          % (1932 - int(round(ss.mean())), 1932 - int(ss.max())))
    print()

    print("[2] descent from the shipped rotation, against a margin of %.2f" % margin)
    best, Q = descend(Nsh.astype(np.float64) / Dsh, iters, margin)
    k, sl = kept(Q)
    print("      %d of 1932 clear the margin; all %d clear the threshold, slack %.6f"
          % (best, k, sl))
    if k < 1932:
        print("      the descent did not reach 1932 on this run; nothing is written")
        return 1
    print()

    print("[3] Cayley, smallest denominator whose D stays inside int64")
    S = (np.eye(n) - Q) @ np.linalg.inv(np.eye(n) + Q)
    S = 0.5 * (S - S.T)
    found = None
    for d in range(4, dmax + 1):
        M = np.rint(S * d).astype(np.int64)
        M = np.triu(M, 1)
        M = [[int(x) for x in row] for row in (M - M.T)]
        X = cayley_exact(M, d)
        N, D = clear_denominators(X)
        if 8 * D >= 2 ** 63 - 1:
            continue
        kk, sle = kept(np.array([[float(x) for x in row] for row in X]))
        if kk != 1932:
            continue
        viol, mx, T = exact_report(N, D)
        print("      d = %-3d D = %-20d exact violations %d   (float slack %.6f)"
              % (d, D, viol, sle))
        if viol == 0:
            found = (N, D, d)
            break
    if found is None:
        print("      no denominator up to %d worked on this run" % dmax)
        return 1
    N, D, d = found
    print()
    print("[4] the result: Q = N/D with d = %d, D = %d" % (d, D))
    viol, mx, T = exact_report(N, D)
    print("      max |<N w_i, w_j>| = %d   threshold isqrt(48 D^2) = %d   violations %d"
          % (mx, T, viol))
    if write:
        np.savetxt(os.path.join(DATA, "axis_rotation.txt"), np.array(N, dtype=np.int64),
                   fmt="%d")
        open(os.path.join(DATA, "axis_rotation_D.txt"), "w").write("%d\n" % D)
        np.savetxt(os.path.join(DATA, "axis_keep.txt"), np.arange(1932), fmt="%d")
        print("      written to data/")
    else:
        same = (np.array(N, dtype=np.int64) == Nsh).all() and D == Dsh
        print("      identical to the shipped rotation: %s" % same)
        print("      (a different one is expected -- the descent is a search, and what the")
        print("       claim rests on is verify.py checking the SHIPPED matrix exactly)")
    return 0


if __name__ == "__main__":
    sys.exit(main(write="--write" in sys.argv))
