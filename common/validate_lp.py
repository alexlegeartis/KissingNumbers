#!/usr/bin/env python3
"""ADVERSARIAL VALIDATION of the k-point moment LP.

    python validate_lp.py            # ~5 min, needs numpy + scipy

`kpoint_lp.py` produces LOWER bounds on cross-section sizes.  A lower bound is worthless if
it is too small, but it is WRONG if it is too large -- and it becomes too large the moment
any single constraint in the program is not actually valid.  That is the failure mode this
file exists to catch, and it catches it the only way that is convincing: by running the
program on every configuration whose true size is known independently, and checking

        LP minimum   <=   the true value   <=   LP maximum.

If any constraint were too strong, some truth would fall outside its own bracket.  Every
cross-section whose size is known independently is checked -- eighteen data points, from four
lattices and three separate sources of truth: direct enumeration from coordinates (5), the
classical laminated-lattice values (7), and other people's published computations (5).  The
eighteenth, Gamma_72 at k = 1, has no external source -- there the moment system has a
unique solution and pins the count outright, so it tests consistency rather than a bound.

WHAT WOULD FAIL, AND WHAT THAT WOULD MEAN

  * a truth BELOW its LP minimum  ->  some constraint is invalid; every bound this
    machinery has ever produced is suspect, including dimensions 46, 47, 70 and 71;
  * a truth ABOVE its LP maximum  ->  the same, in the other direction;
  * the two equivalent Gram matrices disagreeing  ->  the moment family is incomplete,
    which is exactly the bug that was present for a day and moved three results.

THE EQUIVALENT-GRAM TEST is the cheap decisive one and is run first.  The Gram matrices
(3,-3,0) and 3*A_3 satisfy U^T G U = 3*A_3 for an integer U of determinant -1, so they
describe the same configuration and MUST give the same answer.  With the moment family
restricted to all-even exponent vectors -- which is what an earlier revision imposed -- they
gave 6462480 and 5194969.  Any future edit to kpoint_lp.py should be run through this file.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpoint_lp import run3                                # noqa: E402
from theta import onepoint                                # noqa: E402

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


# one-point marginals, each recomputed from Venkov's theorem rather than quoted
def marginals(n, mu, N, t):
    sol, sok, _ = onepoint(n, mu, N, t)
    assert sok, (n, "surplus equations fail")
    return {c: int(sol[c]) for c in sol}


M1E = marginals(8, 2, 240, 7)
M1L = marginals(24, 4, 196560, 11)
M1P = marginals(48, 6, 52416000, 11)
M1G = marginals(72, 8, 6218175600, 11)

print(__doc__)
print("=" * 86)
print("TEST 1  equivalent Gram matrices must give the same answer")
print("=" * 86)
G1 = [[6, 3, -3], [3, 6, 0], [-3, 0, 6]]          # the (3,-3,0) presentation
G2 = [[6, 3, 3], [3, 6, 3], [3, 3, 6]]            # 3 * A_3
U = np.array([[-1, -1, 0], [0, 0, -1], [-1, 0, 0]], dtype=np.int64)
check(int(round(np.linalg.det(U))) in (1, -1), "U is unimodular, det = %d"
      % int(round(np.linalg.det(U))))
check((U.T @ np.array(G1) @ U == np.array(G2)).all(),
      "U^T (3,-3,0) U = 3*A_3 exactly, so the two Grams are the same configuration")
r1 = run3(48, 52416000, 6, [1, 2, 3], M1P, G1)
r2 = run3(48, 52416000, 6, [1, 2, 3], M1P, G2)
check(r1[2] == r2[2], "both give %d and %d" % (r1[2], r2[2]))
print("      (with an all-even moment family they gave 6462480 and 5194969)")

print()
print("=" * 86)
print("TEST 2  every independently known cross-section lies inside its own LP bracket")
print("=" * 86)
print("  %-34s %-13s %-13s %-13s %s" % ("configuration", "LP min", "TRUTH", "LP max", "source of truth"))

CASES = [
    # (label, n, N, mu, IPS, M1, Gram, truth, design order, source)
    ("E_8, k=1", 8, 240, 2, [1], M1E, None, 126, 7, "E_7, classical"),
    ("E_8, k=2 orthogonal", 8, 240, 2, [1], M1E, [[2, 0], [0, 2]], 60, 7, "D_6, classical"),
    ("E_8, k=2 60-degree", 8, 240, 2, [1], M1E, [[2, 1], [1, 2]], 72, 7, "E_6, classical"),
    ("Leech, k=1", 24, 196560, 4, [1, 2], M1L, None, 93150, 11, "Lambda_23, classical"),
    ("Leech, k=2 orthogonal", 24, 196560, 4, [1, 2], M1L, [[4, 0], [0, 4]], 43164, 11,
     "counted from coordinates"),
    ("Leech, k=2 cos 1/4", 24, 196560, 4, [1, 2], M1L, [[4, 1], [1, 4]], 44550, 11,
     "counted from coordinates"),
    ("Leech, k=2 60-degree", 24, 196560, 4, [1, 2], M1L, [[4, 2], [2, 4]], 49896, 11,
     "Lambda_22, classical"),
    # BOTH ends of this interval are attained, and the earlier note here was wrong to call
    # 19962 "the LP's upper end read back as a count".  It is a genuine count, for a
    # different orthogonal triple.  In coordinates where a minimal vector has raw norm 32,
    # 4e1+4e2, 4e3+4e4, 4e5+4e6 gives 19530 while 4e1+4e2, 4e1-4e2, 4e3+4e4 gives 19962;
    # both are pairwise orthogonal of norm mu, so they share the Gram matrix mu*I_3.  The
    # cross-section is therefore NOT a function of the Gram matrix, and this row records one
    # embedding.  leech_truth.py's find_tuple takes a single seeded random tuple, so which
    # value it reports is an accident of the seed.  Sampling 500 orthogonal triples gives
    # 19530 in 97.2% of them and 19962 in the rest; at k = 4 the counts 8616, 8760 and 9192
    # all occur.  The paper's Table 2 says so.
    ("Leech, k=3 orthogonal", 24, 196560, 4, [1, 2], M1L,
     [[4, 0, 0], [0, 4, 0], [0, 0, 4]], 19530, 11, "leech_truth.py, from coordinates"),
    ("Leech, k=3 60-degree triple", 24, 196560, 4, [1, 2], M1L,
     [[4, 2, 2], [2, 4, 2], [2, 2, 4]], 27720, 11, "Lambda_21, classical"),
    # k = 4 on the Leech is what licenses the dimension-68 bound, which is a k = 4
    # cross-section of Gamma_72.  All three counts are leech_truth.py's, from coordinates.
    ("Leech, k=4 orthogonal", 24, 196560, 4, [1, 2], M1L,
     [[4, 0, 0, 0], [0, 4, 0, 0], [0, 0, 4, 0], [0, 0, 0, 4]], 8616, 11,
     "leech_truth.py, from coordinates"),
    ("Leech, k=4 2*A_4", 24, 196560, 4, [1, 2], M1L,
     [[4, 2, 2, 2], [2, 4, 2, 2], [2, 2, 4, 2], [2, 2, 2, 4]], 15540, 11,
     "leech_truth.py, from coordinates"),
    ("Leech, k=4 2*D_4", 24, 196560, 4, [1, 2], M1L,
     [[4, -2, 0, 0], [-2, 4, -2, -2], [0, -2, 4, 0], [0, -2, 0, 4]], 17400, 11,
     "Lambda_20, classical; leech_truth.py agrees"),
    ("P_48, k=1", 48, 52416000, 6, [1, 2, 3], M1P, None, 23766960, 11,
     "Boyvalenkov-Cherkashin eq. (3); Sun-Wang Table 3"),
    ("P_48, k=2 60-degree", 48, 52416000, 6, [1, 2, 3], M1P, [[6, 3], [3, 6]], 12309600, 11,
     "Ozeki degree-3 Siegel; Sun-Wang 3A_2"),
    ("P_48, k=3 (3A_3 embedding)", 48, 52416000, 6, [1, 2, 3], M1P,
     [[6, 3, 3], [3, 6, 3], [3, 3, 6]], 6687648, 11, "Sun-Wang Table 3, complement 3A_3"),
    ("P_48, k=3 (3D_3 embedding)", 48, 52416000, 6, [1, 2, 3], M1P,
     [[6, 3, 3], [3, 6, 3], [3, 3, 6]], 6702080, 11, "Sun-Wang Table 3, complement 3D_3"),
    ("P_48, k=4 (3D_4)", 48, 52416000, 6, [1, 2, 3], M1P,
     [[6, -3, 0, 0], [-3, 6, -3, -3], [0, -3, 6, 0], [0, -3, 0, 6]], 4153600, 11,
     "Sun-Wang Table 3, complement 3D_4"),
    ("Gamma_72, k=1", 72, 6218175600, 8, [1, 2, 3, 4], M1G, None, 2603658750, 11,
     "Venkov, one-point, over-determined"),
]


# The docstring above says how many data points there are.  Counts in prose go stale:
# it is the sentence, not the list, that a reader trusts, so make the list check it.
_WORDS = {'sixteen': 16, 'seventeen': 17, 'eighteen': 18, 'nineteen': 19, 'twenty': 20}
_said = next((v for w, v in _WORDS.items() if w + ' data points' in __doc__), None)
assert _said == len(CASES), (
    'the docstring says %s data points and CASES has %d'
    % (_said, len(CASES)))
nbelow = nabove = 0
for lab, n, N, mu, IPS, M1, G, truth, tdes, src in CASES:
    t0 = time.time()
    if G is None:
        lo = hi = M1[0]                     # k = 1: the system has a unique solution
    else:
        rl = run3(n, N, mu, IPS, M1, G, tdes=tdes)
        rh = run3(n, N, mu, IPS, M1, G, tdes=tdes, sense=-1)
        if rl is None or rh is None or rl[0] != 'cert' or rh[0] != 'cert':
            check(False, "%s: no certificate produced" % lab)
            continue
        lo, hi = rl[2], rh[2]
    inside = lo <= truth <= hi
    if truth < lo:
        nbelow += 1
    if truth > hi:
        nabove += 1
    mark = "ok" if inside else "*** OUTSIDE ***"
    print("  %-34s %-13d %-13d %-13d %s" % (lab, lo, truth, hi, src))
    check(inside, "    %s: %d <= %d <= %d   %s  [%.0fs]"
          % (lab, lo, truth, hi, mark, time.time() - t0))

print()
check(nbelow == 0, "no truth falls BELOW its LP minimum "
                   "(that would make a claimed lower bound invalid)")
check(nabove == 0, "no truth falls ABOVE its LP maximum "
                   "(that would make an 'exact' claim invalid)")

print()
print("=" * 86)
print("TEST 3  the two claims that assert EXACTNESS really have min = max")
print("=" * 86)
for lab, n, N, mu, IPS, M1, G, want in (
        ("P_48, 60-degree pair -> dimension 46", 48, 52416000, 6, [1, 2, 3], M1P,
         [[6, 3], [3, 6]], 12309600),):
    rl = run3(n, N, mu, IPS, M1, G)
    rh = run3(n, N, mu, IPS, M1, G, sense=-1)
    check(rl[2] == rh[2] == want,
          "%s: min = max = %d" % (lab, rl[2]))
print("      Dimension 47 and dimension 71 are exact for a different reason: at k = 1 the")
print("      moment system has a unique solution, so there is nothing to minimise.")
print("      Dimension 70 is NOT exact -- its bracket is checked in its own package.")

print()
print("=" * 86)
print("ALL CHECKS PASSED" if ok else "*** FAILURES -- DO NOT SHIP ***")
print("=" * 86)
raise SystemExit(0 if ok else 1)
