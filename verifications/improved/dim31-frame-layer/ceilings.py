#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ceilings.py --- which factors of this configuration are at a ceiling, and on what evidence

    python ceilings.py 31            (a few seconds)

verify.py certifies the LOWER BOUND.  This script asks the other question -- how much of the
configuration is finished -- and it is written so that every line of the answer says whether it
is a PROOF, an EXHAUSTIVE search, or a measurement.  Nothing here is a sample.

    K(24 + k)  =  196560  +  2 w L  +  |axis|  +  96k          (all classes the same size)

with `L` the size of one owner class and `w = sum_i (|G_i| - 1)` the direction WEIGHT: a group
of three directions is worth 2 per owner line, a group of two 1, a singleton 0.

[1] THE DIRECTION WEIGHT is at `floor(2 tau(k) / 3)`, and that is a theorem.

    Two owner classes that share a direction `z` would meet at `<u,u'> <= 1` for all their
    vectors -- the cap constraint `(2/3)<u,u'> + (4/3)<z,z> <= 2` -- and so would be ONE class.
    So the groups are disjoint.  Inside a group the directions are pairwise at `<= -1/2`, and
    `|sum z|^2 >= 0` then caps a group at THREE.  Across groups they are pairwise at `<= 1/2`.
    So the union of the groups is a 60-degree code in R^k: with `t` groups of three and `q` of
    two, `3t + 2q <= tau(k)`, and maximising `w = 2t + q` over that gives

        w  <=  floor(2 tau(k) / 3)   exactly   (t = floor(tau/3) is optimal, LP value 2 tau/3).

    THE SCOPE.  tau(4) = 24 is PROVED, so dimension 28 is at an unconditional ceiling.
    tau(5), tau(6) and tau(7) are not: 40, 72 and 126 are the best KNOWN kissing numbers and
    44, 77 and 134 the best proven upper bounds, so in dimensions 29, 30 and 31 the ceiling is
    conditional on a 5-, 6- or 7-dimensional kissing number that nobody has proved.  One more
    direction in tau(k) is worth two more owner lines' worth of points here, and the report
    below prices the whole gap.

[2] THE AXIS IS MAXIMAL -- exhaustively, not by sampling.

    A new axis point is a UNIT vector `a` of R^k with `<a,a_i> <= 1/2` against the axis already
    there, `<a,z> <= sqrt3/2` against the cap directions and `|<a,w>| <= 1/sqrt2` against the
    layer frame.  Every right-hand side is positive, so `P = {a : C a <= b}` is a polytope with
    the origin in its INTERIOR, and `|a|` is convex, so

        max_{a in P} |a|  is attained at a VERTEX of P.

    Enumerate the vertices.  If the largest is below 1 then P contains no unit vector at all and
    no point can be added -- over the whole continuum, not over a sample.  The margin here is
    large (13 to 23 per cent), so the floating-point enumeration is not close to a decision; it
    is repeated with the right-hand sides perturbed, and checked against an independent scan
    that needs no enumeration at all (for a direction `d` the ray leaves P at
    `t(d) = 1/max_j (c_j.d / b_j)`, so `max_d t(d)` is the same maximum, from below).

    THE CONTROL.  Drop one axis point and the enumeration must find it again.  It does, and it
    finds NOTHING ELSE: the recovered unit vector is always the point removed.  So the axis is
    maximal, and still maximal after any one-for-one swap.

    THE SCOPE.  This is maximality for THIS axis code against THIS layer frame and THESE cap
    directions.  A different 60-degree code of R^k, or a different rotation of this one, is a
    different polytope; that search is what chose the rotation this package ships.

[3] WHAT IS NOT AT A CEILING is printed at the end: the class, the packing where it is not
    closed, and tau(k) itself where it is not proved.
"""
import numpy as np, sys, os, json, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
DIM = int(sys.argv[1]) if len(sys.argv) > 1 else 31
if DIM not in (29, 30, 31):
    sys.exit("usage: python ceilings.py 29|30|31")

# The best KNOWN kissing numbers and the best PROVEN upper bounds, dimensions 4-7.  Both
# columns are H. Cohn, Table of kissing number bounds; the copy audit.py checks this against
# is common/published.py (COHN and COHN_UPPER).  Dimension 4 is solved, the others are not.
TAU_KNOWN = {4: 24, 5: 40, 6: 72, 7: 126}
TAU_PROVED = {4: 24, 5: 44, 6: 77, 7: 134}

FAILS = []


def req(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg, flush=True)
    if not ok:
        FAILS.append(msg)


def qval(t, den):
    a, b, c, d = t
    return (a + b*np.sqrt(2) + c*np.sqrt(3) + d*np.sqrt(6)) / den


M = json.load(open(os.path.join(HERE, 'data', 'geom%d.json' % DIM)))
k, den = M['k'], M['denominator']
A = np.array([[qval(t, den) for t in a] for a in M['axis']])
Zall = np.array([[qval(t, den) for t in z] for z in M['directions']])
W = np.array([[qval(t, den) for t in w] for w in M['frame']])
GRP = M['groups']
used = sorted({i for g in GRP for i in g})
Z = Zall[used]
L = np.load(os.path.join(HERE, 'data', 'bounds%d.npy' % DIM))
SIZES = np.diff(L)

print("== dimension %d: k = %d, %d classes, %d axis points, %d cap directions of %d =="
      % (DIM, k, len(GRP), len(A), len(Z), len(Zall)))

# --------------------------------------------------------------------------- [1] the weight
print("\n[1] the direction weight")
flat = [i for g in GRP for i in g]
req(len(set(flat)) == len(flat), "the direction groups are pairwise disjoint (%d directions)"
    % len(flat))
req(max(len(g) for g in GRP) <= 3, "no group has more than three directions -- and cannot, since "
                                   "they are pairwise at cosine <= -1/2")
COS = Zall @ Zall.T
req(all(COS[a, b] <= -0.5 + 1e-12 for g in GRP for a, b in itertools.combinations(g, 2)),
    "inside a group the directions are at cosine <= -1/2")
req(all(COS[a, b] <= 0.5 + 1e-12 for i, g in enumerate(GRP) for h in GRP[i+1:]
        for a in g for b in h), "across groups they are at cosine <= 1/2, so the union of the "
                                "groups is a 60-degree code in R^%d" % k)
t3 = sum(1 for g in GRP if len(g) == 3)
q2 = sum(1 for g in GRP if len(g) == 2)
w = 2*t3 + q2
req(w == sum(len(g) - 1 for g in GRP), "weight w = 2t + q = %d  (t = %d triples, q = %d pairs)"
    % (w, t3, q2))
ck, cp = 2*TAU_KNOWN[k]//3, 2*TAU_PROVED[k]//3
req(3*t3 + 2*q2 <= TAU_KNOWN[k], "the groups use %d directions, within tau(%d) = %d"
    % (3*t3 + 2*q2, k, TAU_KNOWN[k]))
req(w == ck, "w = %d IS the ceiling floor(2 tau(%d)/3) = %d for the best KNOWN tau(%d) = %d"
    % (w, k, ck, k, TAU_KNOWN[k]))
if TAU_PROVED[k] == TAU_KNOWN[k]:
    print("   tau(%d) = %d is PROVED, so this ceiling is unconditional." % (k, TAU_KNOWN[k]))
else:
    print("   CONDITIONAL: tau(%d) is only known to lie in [%d, %d], so the UNCONDITIONAL "
          "ceiling is floor(2 x %d/3) = %d, which is %d above the %d realised -- worth %d points "
          "at %d lines a class.  It needs a better %d-dimensional kissing number."
          % (k, TAU_KNOWN[k], TAU_PROVED[k], TAU_PROVED[k], cp, cp - w, w,
             2*(cp - w)*int(SIZES.max()), int(SIZES.max()), k))

# --------------------------------------------------------------------------- [2] the axis
print("\n[2] the axis, exhaustively")
from scipy.spatial import HalfspaceIntersection

C = np.vstack([A, Z, W, -W])
b = np.concatenate([np.full(len(A), 0.5), np.full(len(Z), np.sqrt(3)/2),
                    np.full(2*len(W), 1/np.sqrt(2))])
req(bool((b > 0).all()), "every constraint has a positive right-hand side, so the origin is "
                         "interior and the polytope is the right object")


def maxnorm(Ax, eps=0.0):
    Cc = np.vstack([Ax, Z, W, -W])
    bb = np.concatenate([np.full(len(Ax), 0.5), np.full(len(Z), np.sqrt(3)/2),
                         np.full(2*len(W), 1/np.sqrt(2))]) * (1.0 + eps)
    V = HalfspaceIntersection(np.hstack([Cc, -bb[:, None]]), np.zeros(k)).intersections
    return float(np.linalg.norm(V, axis=1).max()), len(V)


m0, nv = maxnorm(A)
print("   %d vertices;  max |a| over the admissible polytope = %.9f" % (nv, m0))
req(m0 < 1.0, "no unit vector fits: the axis is MAXIMAL (margin %.1f per cent)"
    % (100*(1 - m0)))
mp = [maxnorm(A, e)[0] for e in (-1e-6, 1e-6, 1e-4)]
req(max(mp) < 1.0, "and still maximal with the right-hand sides moved by +-1e-6 and +1e-4 "
                   "(%s), so this is not a floating-point decision"
    % ", ".join("%.6f" % x for x in mp))

rng = np.random.default_rng(0)
D = rng.normal(size=(200000, k))
D /= np.linalg.norm(D, axis=1)[:, None]
t = 1.0 / np.maximum((D @ C.T) / b, 1e-300).max(1)
req(t.max() <= m0 + 1e-9, "an independent 200000-direction scan of the same maximum, needing no "
                          "enumeration, reaches %.9f and never exceeds it" % t.max())

worst_extra = 0
recovered = []
for j in range(len(A)):
    Cc = np.vstack([np.delete(A, j, 0), Z, W, -W])
    bb = np.concatenate([np.full(len(A)-1, 0.5), np.full(len(Z), np.sqrt(3)/2),
                         np.full(2*len(W), 1/np.sqrt(2))])
    V = HalfspaceIntersection(np.hstack([Cc, -bb[:, None]]), np.zeros(k)).intersections
    u = V[np.linalg.norm(V, axis=1) > 1 - 1e-9]
    recovered.append(len(u))
    if len(u):
        u = u / np.linalg.norm(u, axis=1)[:, None]
        worst_extra = max(worst_extra, int((np.abs(u @ A[j]) < 1 - 1e-9).sum()))
req(min(recovered) >= 1, "CONTROL: removing any one axis point makes a unit vector fit again, so "
                         "the test can fail (recovered %d..%d)" % (min(recovered), max(recovered)))
req(worst_extra == 0, "and the recovered vector is ALWAYS the point removed -- 0 others over all "
                      "%d removals, so no one-for-one swap opens the axis either" % len(A))

# --------------------------------------------------------------------------- [3] what is open
print("\n[3] what is NOT at a ceiling")
L0 = int(SIZES.max())
shape = ("all %d of them at %d" % (len(SIZES), L0) if SIZES.min() == SIZES.max() else
         "the %d here run %d..%d, the repeated lines having been dropped -- a subset of a class "
         "is again a class" % (len(SIZES), int(SIZES.min()), int(SIZES.max())))
print("   the CLASS.  The template class has %d lines (%s); one more line in the template is "
      "worth +%d in dimension %d.  248 is maximal over all 48576 type-B lines and rigid, but the "
      "Delsarte LP on the Leech line scheme allows 425.45 lines, and on the 84 octads "
      "transverse to one trio -- where every KNOWN class lives, though no proof says every "
      "class must -- the Hoffman ratio bound allows 373.33.  So this is open."
      % (L0, shape, 2*w, DIM))
if int(SIZES.sum()) < len(GRP)*L0:
    print("   the PACKING.  %d owner lines against the %d that %d full classes would give: "
          "%d repeats, worth up to +%d.  Two neighbourhoods around this family are CLOSED "
          "exhaustively, rather than by a search stopping.  (a) The SIGN WORDS: the overlap of "
          "two images depends on them only through c_i xor c_j, every one of the 861 exact pair "
          "tables reaches 0 somewhere, and annealing still cannot leave %d.  (b) The RE-SLOTTING: "
          "one trio admits exactly 2016 supports and its stabiliser has order 64512, so 256 "
          "images realise each support and they differ in the 4-slot taken inside each light "
          "coset -- which no sign word can reach, since a sign word TRANSLATES a slot.  All "
          "42 x 256 of them over all 4096 sign words is 44 million placements and zero "
          "improvements (research/collab2531/dim31/reslot31.py).  What is open is search DEPTH: "
          "on octad coverage, on the bit budget against the 42 x 12 = 504 bits of sign freedom, "
          "on the pair probability and on the code agreement alike, five fresh builds match this "
          "family and lose twice as many lines (KNOWLEDGE.md section 149)."
          % (int(SIZES.sum()), len(GRP)*L0, len(GRP),
             len(GRP)*L0 - int(SIZES.sum()), 4*(len(GRP)*L0 - int(SIZES.sum())),
             len(GRP)*L0 - int(SIZES.sum())))
else:
    print("   the PACKING is closed: %d x %d = %d owner lines, every class full and disjoint."
          % (len(GRP), L0, int(SIZES.sum())))
if TAU_PROVED[k] != TAU_KNOWN[k]:
    print("   tau(%d) itself, which is not proved (see [1])." % k)

print("\n%s" % ("*** ALL CHECKS PASS ***" if not FAILS else "*** %d CHECKS FAILED ***" % len(FAILS)))
sys.exit(1 if FAILS else 0)
