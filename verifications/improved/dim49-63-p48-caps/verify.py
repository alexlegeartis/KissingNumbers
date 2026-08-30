#!/usr/bin/env python3
"""Exact verification of the explicit data behind dimensions 49-63.

    python verify.py            # ~1 min, needs numpy, all decisions in integer arithmetic

Three things are checked, in this order:

  1. THE LATTICE.  data/p48p_gram.npy is the Gram matrix of P_48p as published in the
     Catalogue of Lattices (Nebe, entry P48p).  Verified: 48x48, symmetric, even diagonal,
     and **determinant exactly 1** by fraction-free Bareiss elimination in Python integers.
     Even plus unimodular plus dimension 48 plus minimum 6 is exactly "extremal even
     unimodular", which is all any of this needs.
     Note the catalogue basis is NOT a minimal-vector basis, so the diagonal is 6^21 8^26
     10^1 rather than all 6.  The minimum is still 6; that is a property of the lattice, not
     of the basis.  Each published group generator A is checked to satisfy A G A^T = G
     exactly, with det A = +-1.

  2. THE CLASS.  data/class_p48_7069.npy holds 7069 coefficient vectors.  Verified: every
     row x has x G x^T = 6, so it is a minimal vector; every pair has |x G y^T| <= 2, which
     is the CLASS condition |cos| <= 1/3; in particular no pair reaches 6, so no two rows
     are +-each other and all 7069 are distinct lines.  The full inner-product histogram is
     printed and checked to sum to C(7069,2).

     7069 lines is what an actual search finds.  Caro-Wei on the conflict graph -- which is
     regular of degree n[3] = 36848 on 26 208 000 lines -- guarantees only 712.  The factor
     9.93 is in line with the Leech's 11.3 and Gamma_72's 9.58.

  3. THE ARITHMETIC.  data/class_sizes_p48.npy holds the sizes of 1400 pairwise disjoint
     classes.  Verified here: it has at least the 1282 that dimension 63 reads, none of
     them empty, and its largest
     equals the class of step 2.  The classes themselves are NOT shipped -- 6.9 million
     lines x 48 coordinates is 142 MB even compressed -- but they are regenerated
     DETERMINISTICALLY and re-verified by scripts/regenerate.py, which asserts pairwise
     disjointness and re-derives exactly this size vector.  Running that script is what
     turns the numbers of final.py into certified ones; this script checks the inputs it
     starts from.  See the README.

WHY P_48 GETS PAIRS AND NOT TRIPLES.  Its class threshold is |<u,u'>| <= 2 of 6, i.e.
gamma = 1/3, which forces cap level t = 1/(2(1-gamma)) = 3/4; and at t = 3/4 two cap
directions on the same line need <z,z'> <= (1/2-t)/(1-t) = -1, i.e. z' = -z.  So a line can
serve an ANTIPODAL PAIR of directions but never a zero-sum triple.  Gamma_72 and the Leech
both have gamma = 1/4, which is exactly the threshold that permits t = 2/3 and hence
triples -- that is the one place they do better.  See ../dim73-95-gamma72-caps/.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


def bareiss_det(M):
    """exact determinant of an integer matrix, fraction-free; no floating point"""
    A = [[int(x) for x in row] for row in M]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = A[i][j] * A[k][k] - A[i][k] * A[k][j]
                assert num % prev == 0, "Bareiss division not exact"
                A[i][j] = num // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


print(__doc__)
print("=" * 78)
print("STEP 1  the lattice P_48p")
print("=" * 78)
G = np.load(os.path.join(DATA, 'p48p_gram.npy')).astype(np.int64)
check(G.shape == (48, 48), "Gram matrix is %dx%d" % G.shape)
check((G == G.T).all(), "symmetric")
d = sorted(set(int(x) for x in np.diag(G)))
check(all(x % 2 == 0 for x in np.diag(G)),
      "diagonal is even, values %s -- so the lattice is EVEN" % d)
det = bareiss_det(G)
check(det == 1, "Bareiss determinant = %d -- so the lattice is UNIMODULAR" % det)
GENS = np.load(os.path.join(DATA, 'p48p_gens.npy')).astype(np.int64)
allgen = True
for i, A in enumerate(GENS):
    isom = bool((A @ G @ A.T == G).all())
    da = bareiss_det(A)
    allgen = allgen and isom and abs(da) == 1
check(allgen, "all %d published generators satisfy A G A^T = G exactly, det A = +-1"
      % len(GENS))
print("      An even unimodular lattice of dimension 48 with minimum 6 is extremal, and")
print("      then it has exactly 52416000 minimal vectors and its minimal shell is a")
print("      spherical 11-design (Venkov).  Nothing below needs more than that.")

print()
print("=" * 78)
print("STEP 2  the 7069-line class")
print("=" * 78)
C = np.load(os.path.join(DATA, 'class_p48_7069.npy')).astype(np.int64)
check(C.shape == (7069, 48), "class file holds %d vectors of length %d" % C.shape)
# exactness: |C_ij| <= 22 and |G_ij| <= 10, so |(C G)_ij| <= 48*22*10 = 10560, and
# |(C G) C^T| <= 48 * 10560 * 22 = 11151360, both far inside float64's exact integer range.
Cf = C.astype(np.float64)
CG = Cf @ G.astype(np.float64)
check(np.abs(CG).max() < 2 ** 50, "intermediate products stay exact in float64 "
                                  "(max |C G| = %d)" % int(np.abs(CG).max()))
nrm = (CG * Cf).sum(axis=1)
check((nrm == 6).all(), "every row has x G x^T = 6 -- all 7069 are minimal vectors")
hist = {}
worst = 0
for a in range(0, len(C), 1024):
    B = (CG[a:a + 1024] @ Cf.T)
    for i in range(len(B)):
        r = B[i].astype(np.int64)
        r[a + i] = 99                                  # exclude the diagonal
        worst = max(worst, int(np.abs(r[r != 99]).max()))
        for v, n in zip(*np.unique(np.abs(r[r != 99]), return_counts=True)):
            hist[int(v)] = hist.get(int(v), 0) + int(n)
hist = {k: v // 2 for k, v in hist.items()}
check(worst <= 2, "largest |off-diagonal inner product| is %d <= 2 -- the class condition "
                  "|cos| <= 1/3" % worst)
check(hist.get(3, 0) == 0 and hist.get(6, 0) == 0,
      "no pair at 3 (cos = 1/2, a conflict) and none at 6 (an antipodal collision)")
tot = sum(hist.values())
check(tot == 7069 * 7068 // 2,
      "inner-product histogram %s sums to C(7069,2) = %d" % (hist, tot))

print()
print("=" * 78)
print("STEP 3  the disjoint classes, and the arithmetic")
print("=" * 78)
S = np.sort(np.load(os.path.join(DATA, 'class_sizes_p48.npy')).ravel())[::-1]
check(len(S) >= 1282, "%d class sizes recorded, and dimension 63 reads lambda(15) = 1282" % len(S))
check(int(S[0]) == len(C), "the largest is %d, the class verified in step 2" % int(S[0]))
check(int(S.min()) > 0 and int(S.sum()) > 6900000,
      "sizes %d..%d, %d lines in total -- %.1f%% of the 26208000 minimal lines, "
      "none empty" % (int(S.min()), int(S.max()), int(S.sum()),
                      100.0 * int(S.sum()) / 26208000))
print("      The classes themselves are regenerated and re-verified by")
print("      scripts/regenerate.py -- deterministic from this Gram, these generators, the")
print("      class of step 2 and a fixed seed.  That script asserts pairwise disjointness")
print("      and re-derives exactly this size vector.")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
print("=" * 78)
raise SystemExit(0 if ok else 1)
