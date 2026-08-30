#!/usr/bin/env python3
"""EXACT verification of the upper bound on the class problem: at most 425 lines.

    python certificate.py            # ~1 min, needs numpy; all decisions in exact arithmetic

THE QUESTION.  How many lines of Leech minimal vectors can be pairwise at |cos| <= 1/4?
This is the number the whole cap family turns on: raising it by one LINE gives +2, +8, +16,
+32, +52, +96, +168 points in dimensions 25...31 respectively -- half of that, +1, +4, +8,
+16, +26, +48, +84, per extra VECTOR, a line being an antipodal pair.  The record is **248** lines
(Ma et al., PackingStar, arXiv:2511.13391); the best upper bound is **425**.

Two independent routes to 425 are checked here, both exactly.

  ROUTE 1 -- the Hoffman ratio bound on the conflict graph.
    The 98280 Leech lines carry a 4-class association scheme, by |<u,u'>| in {32,16,8,0} at
    squared norm 32.  The conflict relation is |<u,u'>| = 16 (cos = +-1/2), of valency 4600.
    Its eigenvalues are the eigenvalues of the 4x4 INTERSECTION MATRIX B_1, an integer
    matrix, and they come out 4600, 1000, 76, -20 -- verified here by checking that the
    characteristic polynomial of B_1 factors exactly as (x-4600)(x-1000)(x-76)(x+20) over
    the integers.  Hoffman then gives

        alpha <= N * (-lambda_min) / (k - lambda_min) = 98280 * 20 / 4620 = 4680/11

    in exact rational arithmetic, so alpha <= 425.

  ROUTE 2 -- an explicit degree-4 polynomial certificate, with no Leech input at all.
    Put f(t) = (4992/11) * t^2 * (t^2 - 1/16).  Expanded in the Gegenbauer basis for
    S^23 (normalised P_k(1) = 1) it is

        f = 1 * P_0 + (5083/77) * P_2 + (27600/77) * P_4,

    all coefficients non-negative and f_0 = 1.  On the admissible inner products of a line
    system with |cos| <= 1/4 -- that is, t^2 in [0, 1/16] -- f is <= 0, with equality only
    at t = 0 and t^2 = 1/16.  Delsarte's inequality then gives

        N <= f(1) / f_0 = 4680/11 = 425.4545...

    Every step of this is checked below in Fraction arithmetic: the Gegenbauer recurrence,
    the expansion coefficients, their non-negativity, the sign of f on the admissible set,
    and the final quotient.  **This route uses no property of the Leech lattice**, so it
    bounds ANY line system in R^24 with |cos| <= 1/4.

WHY THIS FILE EXISTS.  `scheme_lp.py` in this directory builds the same scheme and runs the
Delsarte LP, the Lovasz theta and the ratio bound numerically.  Its intersection numbers are
correct -- it verifies them constant over random base pairs -- but its eigenvalue extraction
is numerically unstable at this size and returns garbage (it reports eigenvalues near
8801.87, -0.0468, 19489.5 where the true ones are 4600, 1000, 76, -20, and a bound of 6.22
where the true one is 425.45).  Do not quote its output.  This script replaces it for the
part that matters.

Reference for the published form of the bound: H. Cohn, Y. Jiao, A. Kumar, S. Torquato,
"Rigidity of spherical codes", Geom. Topol. 15 (2011), section 7, which states it as
|S| <= floor(9360/11) = 850 vectors.
"""
import os
import sys
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from leech import load as leech_load                      # noqa: E402

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


# ======================================================================= route 1
print(__doc__)
print("=" * 78)
print("ROUTE 1  the Hoffman ratio bound on the conflict graph")
print("=" * 78)

A = leech_load().astype(np.int64)
nz = A != 0
first = np.argmax(nz, axis=1)
sgn = np.sign(A[np.arange(len(A)), first])
L = np.unique(A * sgn[:, None], axis=0)        # one representative per line
check(L.shape == (98280, 24), "98280 Leech lines, built from the Golay code: %s" % (L.shape,))


def relrow(i):
    """relation of every line to line i: 0 self, 1 |ip|=16, 2 |ip|=8, 3 orthogonal"""
    d = np.abs(L @ L[i])
    r = np.zeros(len(L), dtype=np.int8)
    r[d == 8] = 2
    r[d == 16] = 1
    r[d == 32] = 0
    r[d == 0] = 3
    return r


r0 = relrow(0)
val = [int((r0 == i).sum()) for i in range(4)]
check(val == [1, 4600, 47104, 46575], "valencies %s" % val)

rng = np.random.default_rng(0)
B = np.zeros((4, 4), dtype=np.int64)
consistent = True
for h in range(4):
    idx = np.nonzero(r0 == h)[0]
    rows = []
    for _ in range(3):                          # three independent base pairs per relation
        v = int(rng.choice(idx))
        rv = relrow(v)
        rows.append([int(((r0 == 1) & (rv == j)).sum()) for j in range(4)])
    consistent = consistent and all(x == rows[0] for x in rows)
    B[h] = rows[0]
check(consistent, "intersection numbers constant over 3 independent base pairs per relation")
print("      B_1 = %s" % B.tolist())

# characteristic polynomial of B_1, exactly, by expanding prod (x - lambda) and comparing
EV = [4600, 1000, 76, -20]
poly = [F(1)]
for lam in EV:                                  # multiply by (x - lam)
    poly = [poly[i] - (lam * poly[i + 1] if i + 1 < len(poly) else 0)
            for i in range(len(poly))] + [F(0)]
    poly = poly[:-1] + [F(1)]
# build coefficients of prod(x - ev) directly instead, to avoid index gymnastics
coef = [F(1)]
for lam in EV:
    new = [F(0)] * (len(coef) + 1)
    for i, c in enumerate(coef):
        new[i + 1] += c
        new[i] += -lam * c
    coef = new                                  # coef[i] is the coefficient of x^i


def charpoly_int(M):
    """characteristic polynomial of a 4x4 integer matrix, exactly (Faddeev-LeVerrier)"""
    n = M.shape[0]
    Mf = [[F(int(x)) for x in row] for row in M]
    I = [[F(1 if i == j else 0) for j in range(n)] for i in range(n)]

    def mul(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    Mk = [[F(0)] * n for _ in range(n)]
    cs = [F(1)]
    for k in range(1, n + 1):
        Mk = [[Mf[i][j] if k == 1 else sum(Mf[i][t] * Mk[t][j] for t in range(n))
               for j in range(n)] for i in range(n)] if k == 1 else mul(Mf, Mk)
        tr = sum(Mk[i][i] for i in range(n))
        c = -tr / k
        cs.append(c)
        for i in range(n):
            Mk[i][i] += c
    # cs = [1, c1, ..., cn] for x^n + c1 x^(n-1) + ... + cn
    return list(reversed(cs))                   # ascending powers


got = charpoly_int(B)
check(got == coef,
      "char poly of B_1 factors exactly as (x-4600)(x-1000)(x-76)(x+20) over Z")
print("      eigenvalues of the conflict relation: %s" % EV)

N, kdeg, lmin = 98280, 4600, -20
hoff = F(N * (-lmin), kdeg - lmin)
check(hoff == F(4680, 11),
      "Hoffman: %d * %d / %d = %s = %.4f" % (N, -lmin, kdeg - lmin, hoff, float(hoff)))
check(int(hoff) == 425, "so at most floor(4680/11) = %d lines" % int(hoff))

# ======================================================================= route 2
print()
print("=" * 78)
print("ROUTE 2  the degree-4 polynomial certificate, with no Leech input")
print("=" * 78)

DIM = 24


def gegenbauer(kmax, n=DIM):
    """P_0..P_kmax for S^(n-1), normalised P_k(1) = 1, as exact coefficient lists."""
    P = [[F(1)], [F(0), F(1)]]
    for k in range(1, kmax):
        a, b, c = F(2 * k + n - 2), F(k), F(k + n - 2)
        shifted = [F(0)] + P[k]                              # t * P_k
        nxt = [(a * (shifted[i] if i < len(shifted) else 0)
                - b * (P[k - 1][i] if i < len(P[k - 1]) else 0)) / c
               for i in range(k + 2)]
        P.append(nxt)
    return P[:kmax + 1]


P = gegenbauer(4)
check(all(sum(p) == 1 for p in P), "Gegenbauer recurrence gives P_k(1) = 1 for k = 0..4")

# f(t) = (4992/11) t^2 (t^2 - 1/16) = (4992/11) t^4 - (312/11) t^2
f = [F(0), F(0), F(-4992, 11 * 16), F(0), F(4992, 11)]
check(f[2] == F(-312, 11), "f(t) = (4992/11) t^4 - (312/11) t^2")

# expand f in the Gegenbauer basis, highest degree first
rem = list(f)
coeffs = {}
for k in (4, 3, 2, 1, 0):
    lead = P[k][k]
    ck = (rem[k] if k < len(rem) else F(0)) / lead
    coeffs[k] = ck
    for i, c in enumerate(P[k]):
        rem[i] -= ck * c
check(all(x == 0 for x in rem), "f expands exactly in the Gegenbauer basis")
print("      f = %s" % "  +  ".join("(%s) P_%d" % (coeffs[k], k)
                                    for k in (0, 2, 4) if coeffs[k] != 0))
check(coeffs[0] == 1 and coeffs[2] == F(5083, 77) and coeffs[4] == F(27600, 77),
      "coefficients are 1, 5083/77, 27600/77 as published")
check(all(coeffs[k] >= 0 for k in coeffs), "every Gegenbauer coefficient is non-negative")
check(coeffs[1] == 0 and coeffs[3] == 0, "odd coefficients vanish, as an even f requires")


def fval(t):
    return sum(c * t ** i for i, c in enumerate(f))


# f <= 0 exactly on the admissible set: t^2 in [0, 1/16].  f = (4992/11) t^2 (t^2 - 1/16),
# and t^2 >= 0 while (t^2 - 1/16) <= 0 there, so the product is <= 0.  Sampled as a check.
worst = max(fval(F(j, 400)) for j in range(-100, 101))
check(worst <= 0, "f(t) <= 0 for every |t| <= 1/4  (max over 201 exact samples: %s)" % worst)
check(fval(F(0)) == 0 and fval(F(1, 4)) == 0 and fval(F(-1, 4)) == 0,
      "f vanishes exactly at t = 0 and t = +-1/4, the attained inner products")
check(fval(F(1)) == F(4680, 11), "f(1) = %s" % fval(F(1)))

bound = fval(F(1)) / coeffs[0]
check(bound == F(4680, 11), "Delsarte: N <= f(1)/f_0 = %s = %.4f" % (bound, float(bound)))
check(int(bound) == 425, "so at most %d lines, by a route that never mentions the Leech "
                         "lattice" % int(bound))

print()
print("=" * 78)
print("Both routes give 4680/11 = 425.4545..., so at most 425 lines.")
print("The record is 248.  Closing that gap is the open problem -- see README.md.")
print("=" * 78)
print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
raise SystemExit(0 if ok else 1)
