#!/usr/bin/env python3
"""Dimensions 70 and 71: cross-sections of Nebe's extremal 72-dimensional lattice.

    tau(71) >= 2 603 658 750       previously 331 737 984   -- a factor 7.85
    tau(70) >= 1 249 778 250       previously 331 737 984   -- a factor 3.77

Both are LOWER BOUNDS on the kissing numbers.  What is exact for dimension 71 is the SIZE
OF THIS PARTICULAR CONFIGURATION: the 71-dimensional cross-section of Gamma_72 by one
minimal vector contains exactly 2 603 658 750 minimal vectors.  tau(71) itself is not
known, unlike tau(8) = 240 and tau(24) = 196560.  For dimension 70 the LP brackets the
cross-section in [1 249 778 250, 1 250 506 441] and only the lower end is claimed.

WHAT "PREVIOUSLY" MEANS HERE, WHICH IS THE EASIEST THING TO GET WRONG.  The published
tables have NO ENTRY for any dimension between 65 and 71.  Cohn's table covers dimensions
1-48 plus 72; the Nebe-Sloane table it replaced lists 1-40, 42, 44, 48, 64, 72, 80, 128 and
nothing in between.  So the best previously published value in dimensions 65-71 is the one
inherited by monotonicity from dimension 64: 331 737 984, the Edel-Rains-Sloane
construction of 1998.  The ceiling is dimension 72's 6 218 175 600.

  Is Edel-Rains-Sloane really not better at n = 70, 71?  Their construction is generic in n
  but was evaluated by its authors only at 32, 36, 40, 44, 64, 80, 128.  It reverse-
  engineers exactly -- both of their large published values reproduce to the digit,
  2^28 + 30828*2048 + 10416*16 + 128 = 331737984 and
  2^30 + 143780*2048 + 20540*16 + 160 = 1368532064 -- and evaluating it at n = 70 and 71,
  generously handing it the dimension-80 constant-weight code size, gives at most
  563 125 646.  Independently, the construction is monotone in n, so it is at most its
  dimension-80 value 1 368 532 064 for every n <= 80; dimension 71's 2 603 658 750 clears
  even that, beating it in EVERY dimension up to 80.

STRUCTURE OF THE ARGUMENT

  Step 0  The theta series of an extremal even unimodular lattice of dimension 72 is the
          unique weight-36 modular form for SL_2(Z) with constant term 1 and vanishing
          coefficients of q^1, q^2, q^3.  Recomputed here from E_4 and Delta: minimum 8 and
          6 218 175 600 minimal vectors, matching Nebe (2012).

  Step 1  VENKOV: for n = 72 = 0 (mod 24) the minimal shell is a spherical 11-design, so
          every even moment up to degree 10 is exact.  For a fixed minimal vector u the
          inner product <u,v> is an integer c with |c| <= 4 (because |u-v|^2 = 16-2c must
          be at least 8), except v = +-u where c = +-8.  FIVE unknowns, SIX equations: the
          system is over-determined and solves consistently in non-negative integers.  The
          surplus equation is a genuine check, and it holds.

  Step 2  Gamma_72 n u^perp is a 71-dimensional lattice of minimum 8 whose minimal vectors
          are exactly the n[0] minimal vectors of Gamma_72 orthogonal to u.  So
          tau(71) >= n[0].

  Step 3  For dimension 70, two minimal vectors at 60 degrees (<u1,u2> = 4).  Such a pair
          EXISTS because n[4] = 127800 > 0, which Step 1 establishes exactly.  The two-point
          LP -- design moments including the mixed ones, exact one-point marginals,
          antipodal symmetry, lattice translations, lattice range -- is solved and an exact
          rational dual certificate extracted.  A sweep over all five possible inner
          products confirms that the 60-degree pair is the best of them.

    python derive.py

Runs in a few seconds; needs numpy and scipy.  Everything that decides anything is exact
rational arithmetic; see common/PROOF-kpoint.md sections 3 and 6.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common'))
from theta import minimum_and_kissing, onepoint          # noqa: E402
from kpoint_lp import run3                               # noqa: E402
from published import floor_for                          # noqa: E402

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
                                 (24, 4, 196560, 11, 93150, "Leech -> Lambda_23"),
                                 (48, 6, 52416000, 11, 23766960, "P_48  -> dim 47")):
    sol, sok, ns = onepoint(n, mu, N, t)
    check(sol[0] == truth and sok,
          "%-20s one-point n[0] = %-11s known %-11s (%d surplus eq(s) hold: %s)"
          % (name, sol[0], truth, ns, sok))
M1L = {0: 93150, 1: 47104, 2: 4600}
for G, name, truth in (([[4, 0], [0, 4]], "Leech, orthogonal pair", 43164),
                       ([[4, 2], [2, 4]], "Leech, 60-degree pair", 49896)):
    r = run3(24, 196560, 4, [1, 2], M1L, G)
    check(r and r[0] == 'cert' and r[2] == truth,
          "%-22s two-point LP = %-9s true %-9s" % (name, r[2] if r else None, truth))
for G, name, truth in (([[2, 0], [0, 2]], "E_8, orthogonal pair", 60),
                       ([[2, 1], [1, 2]], "E_8, 60-degree pair", 72)):
    r = run3(8, 240, 2, [1], {0: 126, 1: 56}, G, tdes=7)
    check(r and r[0] == 'cert' and r[2] == truth,
          "%-22s two-point LP = %-9s true %-9s  (D_6 and E_6)"
          % (name, r[2] if r else None, truth))

print()
print("=" * 78)
print("STEP 0  the extremal theta series in dimension 72")
print("=" * 78)
mu, N = minimum_and_kissing(72)
check((mu, N) == (8, 6218175600),
      "minimum %d, kissing number %d  (Nebe 2012: 8 and 6218175600)" % (mu, N))

print()
print("=" * 78)
print("STEP 1  the one-point distribution of Gamma_72's minimal shell")
print("=" * 78)
sol, sok, ns = onepoint(72, 8, 6218175600, 11)
for c in sorted(sol):
    print("      n[%+d] = n[%+d] = %d" % (c, -c, sol[c]))
print("      n[+8] = n[-8] = 1   (v = +-u)")
total = sum((2 if c else 1) * sol[c] for c in sol) + 2
check(all(v.denominator == 1 and v >= 0 for v in sol.values()),
      "the solution is integral and non-negative")
check(total == N, "the counts total %d = N" % total)
check(sok, "%d surplus moment equation(s) hold -- a genuine consistency check" % ns)

print()
print("=" * 78)
print("STEP 2  DIMENSION 71")
print("=" * 78)
D71 = int(sol[0])
prev71, why71 = floor_for(71)
check(D71 == 2603658750, "the 71-dimensional cross-section has n[0] = %d minimal vectors"
      % D71)
print("      previously: %d" % prev71)
print("        (%s)" % why71)
print("      tau(71) >= %d       a factor %.2f" % (D71, D71 / prev71))
check(D71 <= 6218175600, "at most tau(72) = 6218175600, as monotonicity requires")

print()
print("=" * 78)
print("STEP 3  DIMENSION 70, and the sweep over all five Gram matrices")
print("=" * 78)
M1G = {c: int(sol[c]) for c in sol}
check(M1G[4] > 0, "a 60-degree pair of minimal vectors EXISTS: n[4] = %d > 0 "
                  "(realizability, and it is Step 1 that proves it)" % M1G[4])
best = (None, -1)
for ip in (0, 1, 2, 3, 4):
    r = run3(72, 6218175600, 8, [1, 2, 3, 4], M1G, [[8, ip], [ip, 8]])
    v = r[2] if r and r[0] == 'cert' else None
    print("      <u1,u2> = %d  (cos = %s):  certified %s" % (ip, ip / 8, v))
    if v is not None and v > best[1]:
        best = (ip, v)
check(best[0] == 4, "the 60-degree pair is the best of the five, at %d" % best[1])
D70 = best[1]
check(D70 == 1249778250, "certified minimum %d" % D70)
hi = run3(72, 6218175600, 8, [1, 2, 3, 4], M1G, [[8, 4], [4, 8]], sense=-1)
print("      LP bracket for that cross-section: [%d, %d] -- unlike dimension 46 the two "
      "ends" % (D70, hi[2]))
print("      do not meet, so only the lower one is claimed.")
prev70, why70 = floor_for(70)
print()
print("      previously: %d" % prev70)
print("        (%s)" % why70)
print("      tau(70) >= %d       a factor %.2f" % (D70, D70 / prev70))
check(D70 <= D71, "at most the dimension-71 value, as monotonicity requires")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
print("=" * 78)
raise SystemExit(0 if ok else 1)
