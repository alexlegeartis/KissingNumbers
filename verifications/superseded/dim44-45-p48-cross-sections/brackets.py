#!/usr/bin/env python3
"""Dimensions 44 and 45: WITHDRAWN, and why the withdrawal is itself worth keeping.

This project certified

    tau(45) >= 6 687 320        (P_48, k = 3, Gram 3*A_3)
    tau(44) >= 4 147 967        (P_48, k = 4, Gram 3*D_4)

Note the dimension-44 Gram: 3*D_4, not 3*A_4.  They are different lattices and give
different cross-sections -- 3*A_4 certifies only 3 694 470 -- and Sun-Wang's dimension-44
complement is the D_4 one.  Comparing against the wrong Gram makes their value look like it
falls OUTSIDE the bracket, which is how the mistake announces itself.

Both are valid lower bounds and both are BEATEN by Sun-Wang, arXiv:2607.20359v3 (18 Aug
2026), Table 3, who compute the twelve Conway-Sloane (1982) cross-sections of P_48 directly
by sweeping all 52 416 000 minimal vectors:

    dim 45   3A_3 -> 6 687 648 ,  3D_3 -> 6 702 080          (and an antipode packing at 7 379 838)
    dim 44   3D_4 -> 4 153 600

So dimensions 44 and 45 are withdrawn.  What survives is a sharp INDEPENDENT CHECK on the
k-point moment machinery: this script recomputes the LP brackets and confirms that every
one of Sun-Wang's values lies inside them.  A bound that brackets an independently computed
truth, from a completely different method, is evidence the machinery is sound -- which
matters, because the same machinery carries dimensions 46, 47, 70 and 71.

TWO THINGS THE COMPARISON SETTLES

  * The count depends on the EMBEDDING, not just the abstract Gram matrix.  Sun-Wang's
    3A_3 and 3D_3 are the same abstract lattice yet give 6 687 648 and 6 702 080.  The LP
    sees only the Gram, so it can never do better than bracket them -- and the bracket here
    is [6 687 320, 6 702 736], which contains both.  The width of that bracket is not slack
    in the method; it is the real spread of the geometry.

  * Ozeki's Theorem 6.3 falls out.  With the mixed moments the three-point system for a
    60-degree triple has affine solution space of dimension exactly 1, along the line
        n[0,0,0] + 164 * n[3,3,3] = 6732912.
    Its two ends are the bracket above, and the corresponding n[3,3,3] range 184..278 is
    Ozeki's bound on tau = a(T4,3)/a(T3,1) -- obtained here by linear programming with an
    exact rational dual certificate rather than by modular forms.

    python brackets.py            # ~1 min, needs numpy + scipy
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common'))
from kpoint_lp import run3                               # noqa: E402
from published import COHN                               # noqa: E402

# P_48 one-point marginals, forced by Venkov's theorem; see common/theta.py
M1P = {0: 23766960, 1: 12608784, 2: 1678887, 3: 36848}
N, MU, DIM = 52416000, 6, 48

# 3*A_3, the Gram of three minimal vectors pairwise at 60 degrees.  A_3 = D_3 as abstract
# lattices, so this single Gram covers BOTH of Sun-Wang's dimension-45 complements, 3A_3 and
# 3D_3 -- which is exactly why the LP can only bracket them: it sees the Gram, and they
# differ only in how they sit inside P_48.
A3 = [[6, 3, 3], [3, 6, 3], [3, 3, 6]]
# 3*D_4.  NOT 3*A_4: Sun-Wang's dimension-44 complement is the D_4 root lattice scaled to
# minimum 6, whose Gram has off-diagonals in {0,-3}, and its cross-section is a different
# number from A_4's.  Getting this wrong is easy and was got wrong once here.
D4 = [[6, -3, 0, 0], [-3, 6, -3, -3], [0, -3, 6, 0], [0, -3, 0, 6]]
A4 = [[6, 3, 3, 3], [3, 6, 3, 3], [3, 3, 6, 3], [3, 3, 3, 6]]

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


print(__doc__)
print("=" * 78)
print("THE LP BRACKETS")
print("=" * 78)
res = {}
for G, k, dim, name in ((A3, 3, 45, "3*A_3 = 3*D_3"), (D4, 4, 44, "3*D_4")):
    lo = run3(DIM, N, MU, [1, 2, 3], M1P, G)
    hi = run3(DIM, N, MU, [1, 2, 3], M1P, G, sense=-1)
    res[dim] = (lo[2], hi[2])
    print("  dim %d  (k = %d, Gram %-13s):  [%d, %d]" % (dim, k, name, lo[2], hi[2]))
    print("           Cohn's table %d;  this project claimed the lower end" % COHN[dim])

print()
print("=" * 78)
print("SUN-WANG'S VALUES, AND WHETHER THEY LIE INSIDE")
print("=" * 78)
SW = [(45, "3A_3", 6687648), (45, "3D_3", 6702080), (44, "3D_4", 4153600)]
print("  (dimension 45's two entries share one Gram matrix -- A_3 and D_3 are the same")
print("   abstract lattice -- so the LP cannot separate them and can only bracket both.)")
for dim, name, v in SW:
    lo, hi = res[dim]
    check(lo <= v <= hi, "dim %d, complement %-5s: %d is inside [%d, %d]"
          % (dim, name, v, lo, hi))
    check(v > lo, "        and it beats this project's %d, so the claim is withdrawn" % lo)
print()
print("  Sun-Wang further report an ANTIPODE packing with kissing number 7379838 in")
print("  dimension 45, which is outside the bracket -- correctly, because that")
print("  configuration is not a cross-section of P_48 at all.")

print()
print("=" * 78)
print("OZEKI'S THEOREM 6.3, REDERIVED BY LINEAR PROGRAMMING")
print("=" * 78)
t3lo = run3(DIM, N, MU, [1, 2, 3], M1P, A3, target=(3, 3, 3))
t3hi = run3(DIM, N, MU, [1, 2, 3], M1P, A3, target=(3, 3, 3), sense=-1)
lo45, hi45 = res[45]
print("  n[3,3,3] certified in [%d, %d]" % (t3lo[2], t3hi[2]))
check(6732912 - 164 * 278 == lo45,
      "6732912 - 164*278 = %d = the low end of the dimension-45 bracket"
      % (6732912 - 164 * 278))
check(6732912 - 164 * 184 == hi45,
      "6732912 - 164*184 = %d = the high end of the dimension-45 bracket"
      % (6732912 - 164 * 184))
check(t3hi[2] == 278, "the certified maximum of n[3,3,3] is exactly 278, Ozeki's upper end")
check(t3lo[2] == 184, "the certified minimum of n[3,3,3] is exactly 184, Ozeki's lower end")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
print("=" * 78)
raise SystemExit(0 if ok else 1)
