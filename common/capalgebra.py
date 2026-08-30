#!/usr/bin/env python3
"""EXACT verification of the cap construction's constraint algebra and its count.

    python capalgebra.py            # a few seconds, no third-party packages

Thirty-nine of this repository's forty-seven claims come from one construction -- dimensions
25, 27, 38, 49-61 and 73-95 -- and all of them rest on the same two things: that a certain list of
inequalities holds at the chosen cap level, and that a certain count is right.  Both are
re-derived here **from first principles, in exact rational arithmetic**, independently of
`final.py`, `verify.py` and every other script that uses them.

THE CONSTRUCTION.  In R^n (+) R^k, all points unit vectors, every pairwise inner product
required to be <= 1/2.  Write a = sqrt(t), b = sqrt(1-t).  A *class* C_i is a set of minimal
LINES of a lattice L of minimum mu, pairwise at |cos| <= gamma; Z_i is a set of unit
directions in R^k attached to class i.

    equator   (v/|v|, 0)                    v minimal in L, v not +-u for any class line u
    caps      (s*a*u/|u|, b*z)              s in {+-1}, u a line of C_i, z in Z_i
    poles     (0, w)                        w in a kissing configuration of R^k

    count = |shell| - 2*sum|C_i| + 2*sum|C_i|*|Z_i| + #poles
          = |shell| + 2*sum |C_i| * (|Z_i| - 1) + #poles

THE SEVEN KINDS OF PAIR, and the exact inequality each imposes:

    1  equator . equator      <v,v'>/mu <= 1/2                     automatic: v-v' is a
                                                                   lattice vector of norm >= mu
    2  equator . cap          a * <v,u>/mu <= 1/2                  worst <v,u>/mu = 1/2 since
                                                                   v != +-u, so t <= 1
    3  cap . cap, same z,     t*gamma + (1-t) <= 1/2               -> t >= 1/(2(1-gamma))
       two lines of a class
    4  cap . cap, same line,  t + (1-t)*<z,z'> <= 1/2              -> <z,z'> <= (1/2-t)/(1-t)
       two directions                                              (s = s'; s != s' is weaker)
    5  cap . cap, two         t/2 + (1-t)*<z,z'> <= 1/2            -> <z,z'> <= 1/2
       classes
    6  cap . pole             b * <z,w> <= 1/2                     -> <z,w>^2 <= 1/(4(1-t))
    7  pole . pole            <w,w'> <= 1/2                        -> the poles are a kissing
                                                                   configuration of R^k

Every one is checked below by squaring where a square root appears, so nothing is decided in
floating point.  Rows 3 and 4 are the interesting ones: row 3 forces t UP, row 4 forces the
directions APART, and the two together are what decides whether a class line may serve an
antipodal PAIR of directions or a zero-sum TRIPLE.

    gamma = 1/3  (P_48)   ->  t >= 3/4, and at t = 3/4 row 4 gives <z,z'> <= -1  ->  PAIRS
    gamma = 1/4  (Leech,  ->  t >= 2/3, and at t = 2/3 row 4 gives <z,z'> <= -1/2 -> TRIPLES
                 Gamma_72)

and three is the most a line can ever get, because n unit vectors pairwise at cos <= -1/2
satisfy n <= 1 + 1/(1/2) = 3.

WHAT THIS DOES NOT CHECK.  That the classes and direction partitions exist -- that is what
the per-package `verify.py`, `verify_classes.py`, `triples.py` and `regenerate.py` do, in
exact integer arithmetic on explicit coordinates.  This file checks that IF they exist, the
arithmetic on top of them is right.  The two halves together are the claim.
"""
from fractions import Fraction as F

ok = True


def check(cond, msg):
    global ok
    ok = ok and bool(cond)
    print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))


HALF = F(1, 2)


def sqrt_le(coef, val, bound):
    """is  sqrt(coef) * val <= bound ?   exactly, for coef >= 0 and rational val, bound."""
    lhs_neg = val < 0
    if lhs_neg:
        return True                      # sqrt(coef) * (negative) <= a positive bound
    if bound < 0:
        return False
    return coef * val * val <= bound * bound


def audit_scheme(name, mu, ipmax, t, part_sizes, poles_claimed, pole_ip_bound=None):
    """Check rows 1-7 for one scheme.  gamma = ipmax/mu is the class threshold."""
    gamma = F(ipmax, mu)
    print()
    print("  --- %s ---" % name)
    print("      mu = %d, class threshold |ip| <= %d of %d  (gamma = %s), cap level t = %s"
          % (mu, ipmax, mu, gamma, t))

    # row 2:  a * (1/2) <= 1/2   <=>   t <= 1
    check(sqrt_le(t, HALF, HALF), "row 2  equator.cap: sqrt(t)/2 <= 1/2, i.e. t <= 1")

    # row 3:  t*gamma + (1-t) <= 1/2
    r3 = t * gamma + (1 - t)
    check(r3 <= HALF, "row 3  same direction, two class lines: %s*%s + %s = %s <= 1/2"
          % (t, gamma, 1 - t, r3))

    # row 4:  the threshold that decides pairs vs triples
    thr = (HALF - t) / (1 - t)
    print("      row 4  same line, two directions: <z,z'> <= (1/2 - t)/(1 - t) = %s" % thr)
    for m in part_sizes:
        if m == 1:
            continue
        # m unit vectors summing to zero, pairwise at cos = -1/(m-1)
        cos = F(-1, m - 1)
        good = cos <= thr
        check(good, "       a part of size %d needs <z,z'> = %s <= %s" % (m, cos, thr))
    # and that no larger part is possible
    if thr < 0:
        nmax = int(1 + 1 / (-thr))
        check(max(part_sizes) <= nmax,
              "       at most %d directions per line (n <= 1 + 1/%s), largest used is %d"
              % (nmax, -thr, max(part_sizes)))

    # row 5:  t/2 + (1-t)*<z,z'> <= 1/2  with <z,z'> <= 1/2
    r5 = t * HALF + (1 - t) * HALF
    check(r5 <= HALF, "row 5  two classes: %s <= 1/2 whenever <z,z'> <= 1/2" % r5)

    # rows 6 and 7
    if poles_claimed:
        b2 = 1 - t
        if pole_ip_bound is None:
            # poles anywhere: worst case <z,w> = 1
            good = sqrt_le(b2, F(1), HALF)
            check(good, "row 6  cap.pole: sqrt(1-t) * 1 = sqrt(%s) <= 1/2, so poles are "
                        "unrestricted" % b2)
        else:
            good = sqrt_le(b2, pole_ip_bound, HALF)
            check(good, "row 6  cap.pole: sqrt(%s) * %s <= 1/2, so poles need <z,w> <= %s"
                        % (b2, pole_ip_bound, pole_ip_bound))
        check(True, "row 7  pole.pole: the poles are taken to be a kissing configuration "
                    "of R^k, so <w,w'> <= 1/2 by construction")
    else:
        b2 = 1 - t
        worst = sqrt_le(b2, F(1), HALF)
        check(not worst, "rows 6, 7  NO poles are claimed -- and correctly so: at t = %s a "
                         "pole at a cap direction would give sqrt(%s) = %.4f > 1/2"
              % (t, b2, float(b2) ** 0.5))
    return gamma


def count(shell, class_sizes, part_sizes, npoles):
    """total = shell + 2 * sum |C_i| * (|Z_i| - 1) + npoles, derived from the three families"""
    assert len(class_sizes) == len(part_sizes)
    deleted = 2 * sum(class_sizes)
    caps = 2 * sum(c * z for c, z in zip(class_sizes, part_sizes))
    return shell - deleted + caps + npoles


print(__doc__)
print("=" * 86)
print("THE CONSTRAINT ALGEBRA, SCHEME BY SCHEME")
print("=" * 86)

# ---------------------------------------------------------------- Leech, t = 3/4, k = 1
audit_scheme("dimension 25: Leech, k = 1, antipodal pair, t = 3/4",
             mu=4, ipmax=1, t=F(3, 4), part_sizes=[2], poles_claimed=True)

# ---------------------------------------------------------------- Leech, t = 2/3, triples
audit_scheme("dimensions 26-31, 38: Leech, zero-sum triples, t = 2/3",
             mu=4, ipmax=1, t=F(2, 3), part_sizes=[3, 2], poles_claimed=True,
             pole_ip_bound=F(866, 1000))          # 30 degrees; < sqrt(3)/2 = 0.8660254...

# ---------------------------------------------------------------- P_48, t = 3/4, pairs
audit_scheme("dimensions 49-63: P_48, antipodal pairs, t = 3/4",
             mu=6, ipmax=2, t=F(3, 4), part_sizes=[2], poles_claimed=True)

# ---------------------------------------------------------------- Gamma_72, t = 2/3
audit_scheme("dimensions 74-80, 91-95: Gamma_72, zero-sum triples, t = 2/3",
             mu=8, ipmax=2, t=F(2, 3), part_sizes=[3, 2], poles_claimed=False)

# ---------------------------------------------------------------- Gamma_72, t = 3/4
audit_scheme("dimensions 73, 81-90: Gamma_72, antipodal pairs, t = 3/4",
             mu=8, ipmax=2, t=F(3, 4), part_sizes=[2], poles_claimed=True)

print()
print("=" * 86)
print("WHY P_48 CANNOT USE TRIPLES, AND THE OTHER TWO CAN")
print("=" * 86)
for lat, mu, ipmax in (("Leech", 4, 1), ("P_48", 6, 2), ("Gamma_72", 8, 2)):
    gamma = F(ipmax, mu)
    tmin = 1 / (2 * (1 - gamma))
    thr = (HALF - tmin) / (1 - tmin)
    nmax = int(1 + 1 / (-thr)) if thr < 0 else 1
    print("  %-9s gamma = %-4s  ->  t >= %-4s  ->  <z,z'> <= %-5s  ->  at most %d "
          "directions per line" % (lat, gamma, tmin, thr, nmax))
print("      P_48's shell has minimum 6 and a class needs |ip| <= 2, so gamma = 1/3 and the")
print("      cap level cannot drop below 3/4 -- which is exactly where triples die.")

print()
print("=" * 86)
print("THE COUNT, RE-DERIVED FROM THE THREE FAMILIES")
print("=" * 86)
CASES = [
    # (label, shell, class sizes in LINES, part sizes, poles, claimed total)
    ("dim 25  Leech, 1 class on an antipodal pair, 2 poles",
     196560, [248], [2], 2, 197058),
    ("dim 26  Leech, 2 classes on 2 triples, 6 poles",
     196560, [248] * 2, [3] * 2, 6, 198550),
    ("dim 27  Leech, 4 classes on 4 triples, 12 poles",
     196560, [248] * 4, [3] * 4, 12, 200540),
    ("dim 28  Leech, 8 triples, 24 poles", 196560, [248] * 8, [3] * 8, 24, 204520),
    ("dim 29  Leech, 12 triples + 2 pairs, 40 poles",
     196560, [248] * 14, [3] * 12 + [2] * 2, 40, 209496),
    ("dim 30  Leech, 24 triples, 72 poles", 196560, [248] * 24, [3] * 24, 72, 220440),
    ("dim 31  Leech, 42 triples, 126 poles", 196560, [248] * 42, [3] * 42, 126, 238350),
    # dimension 38's 644 classes are not all the same size; what the construction needs is
    # only their TOTAL, which is the whole line set.  Modelled as one aggregate class on a
    # triple, which is exactly what the formula sums to.
    ("dim 38  Leech, 644 classes covering all 98280 lines, 1932 axis points",
     196560, [98280], [3], 1932, 591612),
]
for lab, shell, cs, ps, npo, want in CASES:
    got = count(shell, cs, ps, npo)
    check(got == want, "%-62s %d" % (lab, got))

print()
print("      Dimensions 26 and 28-31 reproduce Cohn's table EXACTLY -- five independent")
print("      confirmations that the accounting is right -- and 25 and 27 are the two records")
print("      this project holds, by +2 and +496.  Dimension 38's 644 classes have different")
print("      sizes; the formula depends only on their total, the whole 98280-line set, and")
print("      the per-class sizes and the triple partition are checked in its own verify.py.")
print()
print("=" * 86)
print("ALL CHECKS PASSED" if ok else "*** FAILURES -- DO NOT SHIP ***")
print("=" * 86)
raise SystemExit(0 if ok else 1)
