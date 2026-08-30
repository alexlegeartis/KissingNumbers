#!/usr/bin/env python3
"""Dimensions 46 and 47: cross-sections of an extremal even unimodular lattice of
dimension 48.

    tau(47) >= 23766960        Cohn's table: 9741412     a factor 2.44
    tau(46) >= 12309600        Cohn's table: 5318060     a factor 2.31

Both are LOWER BOUNDS on the kissing numbers.  What is exact is the SIZE OF THE
CONFIGURATION: the 47-dimensional cross-section of P_48 by one minimal vector contains
exactly 23766960 minimal vectors, and the 46-dimensional cross-section by a 60-degree pair
contains exactly 12309600.  tau(46) and tau(47) themselves are not known, unlike
tau(8) = 240 and tau(24) = 196560.

Neither number is new to the literature -- see the README -- but Cohn's table carries
neither.

STRUCTURE OF THE ARGUMENT

  Step 0  The theta series of an extremal even unimodular lattice of dimension 48 is the
          unique weight-24 modular form for SL_2(Z) with constant term 1 and vanishing
          coefficient of q; it gives minimum 6 and N = 52416000 minimal vectors.  Recomputed
          here from E_4 and Delta, exactly.  Four such lattices are known -- P_48p, P_48q
          (Conway-Sloane), P_48m, P_48n (Nebe) -- and nothing below distinguishes them.

  Step 1  VENKOV: the minimal shell is a spherical 11-design, so every even moment up to
          degree 10 is exact.  Distinct minimal vectors have <v,w> in {0,+-1,+-2,+-3} (and
          +-6 for w = +-v), because |v-w|^2 = 12 - 2<v,w> is an even integer >= 6.  Four
          unknowns, six equations: over-determined, and it solves consistently in
          non-negative integers.  n[0] = 23766960 is the dimension-47 cross-section.

  Step 2  For dimension 46 the same theorem is applied to PRODUCTS of inner products with
          two fixed minimal vectors, which is genuine two-point information, together with
          the exact one-point marginals, antipodal symmetry, lattice translations and
          lattice-range constraints.  Minimising -- and separately maximising -- the
          all-zero cell by linear programming with an exact rational dual certificate gives
          the same number both ways, so the cross-section size is exact.

  Step 3  Realizability.  A bound for a prescribed Gram matrix is worth nothing unless two
          minimal vectors with that Gram exist.  Here they do, and Step 1 proves it: the
          60-degree pair needs n[3] > 0, and n[3] = 36848.

    python derive.py

Runs in a few seconds; needs numpy and scipy (for the floating-point LP that finds a dual
vector -- the bound itself is then certified in exact rational arithmetic, so the solver
cannot corrupt it).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common'))
from theta import minimum_and_kissing, onepoint          # noqa: E402
from kpoint_lp import run3                               # noqa: E402
from published import COHN                               # noqa: E402

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


print(__doc__)
print("=" * 78)
print("VALIDATION -- the identical code on cross-sections whose size is known")
print("=" * 78)
for n, mu, N, t, truth, name in ((8, 2, 240, 7, 126, "E_8   -> E_7"),
                                 (24, 4, 196560, 11, 93150, "Leech -> Lambda_23")):
    sol, sok, ns = onepoint(n, mu, N, t)
    check(sol[0] == truth and sok,
          "%-20s one-point n[0] = %-9s known %-9s (%d surplus equation(s) hold: %s)"
          % (name, sol[0], truth, ns, sok))
M1L = {0: 93150, 1: 47104, 2: 4600}
for G, name, truth in (([[4, 0], [0, 4]], "Leech, orthogonal pair", 43164),
                       ([[4, 1], [1, 4]], "Leech, cos 1/4 pair", 44550),
                       ([[4, 2], [2, 4]], "Leech, 60-degree pair", 49896)):
    r = run3(24, 196560, 4, [1, 2], M1L, G)
    check(r and r[0] == 'cert' and r[2] == truth,
          "%-22s two-point LP = %-9s true %-9s" % (name, r[2] if r else None, truth))
print("   (49896 is Lambda_22, kissing number of the 22-dimensional Leech cross-section;")
print("    43164 and 44550 are counted directly from explicit Leech coordinates.)")

print()
print("=" * 78)
print("STEP 0  the extremal theta series in dimension 48")
print("=" * 78)
mu, N = minimum_and_kissing(48)
check((mu, N) == (6, 52416000),
      "minimum %d, kissing number %d  (expected 6 and 52416000)" % (mu, N))

print()
print("=" * 78)
print("STEP 1  the one-point distribution, and DIMENSION 47")
print("=" * 78)
sol, sok, ns = onepoint(48, 6, 52416000, 11)
for c in sorted(sol):
    print("      n[%+d] = n[%+d] = %d" % (c, -c, sol[c]))
print("      n[+6] = n[-6] = 1   (v = +-u)")
total = sum((2 if c else 1) * sol[c] for c in sol) + 2
check(all(v.denominator == 1 and v >= 0 for v in sol.values()),
      "the solution is integral and non-negative")
check(total == N, "the counts total %d = N" % total)
check(sok, "%d surplus moment equation(s) hold -- a genuine consistency check" % ns)
D47 = int(sol[0])
check(D47 == 23766960, "n[0] = %d" % D47)
print()
print("      tau(47) >= %d       Cohn's table: %d, a factor %.2f"
      % (D47, COHN[47], D47 / COHN[47]))

print()
print("=" * 78)
print("STEP 2 and 3  the two-point LP, and DIMENSION 46")
print("=" * 78)
M1P = {0: int(sol[0]), 1: int(sol[1]), 2: int(sol[2]), 3: int(sol[3])}
check(M1P[3] > 0, "a 60-degree pair of minimal vectors EXISTS: n[3] = %d > 0 "
                  "(realizability, and it is Step 1 that proves it)" % M1P[3])
lo = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 3], [3, 6]])
hi = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 3], [3, 6]], sense=-1)
check(lo and lo[0] == 'cert', "minimum of the all-zero cell, exact certificate: %s" % lo[2])
check(hi and hi[0] == 'cert', "maximum of the all-zero cell, exact certificate: %s" % hi[2])
check(lo[2] == hi[2], "minimum = maximum = %d, so the cross-section size is EXACT" % lo[2])
D46 = lo[2]
check(D46 == 12309600, "the 46-dimensional cross-section has %d minimal vectors" % D46)
print()
print("      tau(46) >= %d       Cohn's table: %d, a factor %.2f"
      % (D46, COHN[46], D46 / COHN[46]))
print()
print("      For comparison, the ORTHOGONAL pair gives a smaller cross-section:")
orth = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 0], [0, 6]])
orthmax = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 0], [0, 6]], sense=-1)
print("        orthogonal pair: [%d, %d]  -- not pinned, and Ozeki's degree-3 table shows"
      % (orth[2], orthmax[2]))
print("        it genuinely varies from pair to pair, so the interval cannot be collapsed.")
check(orth[2] < D46, "the 60-degree pair really is the better rank-2 cross-section")

print()
print("=" * 78)
print("INDEPENDENT CONFIRMATION of 12309600, from modular forms")
print("=" * 78)
a3330 = 23775066324172800000       # Ozeki, Tsukuba J. Math. 40 (2016), Table 3
a333 = 1931424768000
check(a3330 % a333 == 0 and a3330 // a333 == D46,
      "Ozeki's degree-3 Siegel Fourier coefficients: %d / %d = %d, exactly and with no "
      "remainder" % (a3330, a333, a3330 // a333))
print("      That computation uses modular forms and no linear programming at all.")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
print("=" * 78)
raise SystemExit(0 if ok else 1)
