#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ceilings.py --- which factors of this configuration are at a ceiling, and on what evidence

    python ceilings.py            (a few seconds)

verify28.py certifies the LOWER BOUND.  This script asks the other question -- how much of the
configuration is finished -- and every line of the answer says whether it is a PROOF, an
EXHAUSTIVE search, or a measurement.  Nothing here is a sample.  Dimension 28 is the one member
of the family 28-31 where BOTH answers are unconditional, because tau(4) = 24 is a theorem
(Musin) where tau(5), tau(6) and tau(7) are not.

    K(28)  =  196560  +  2 w L  +  |axis|  +  96 k ,     k = 4, w = 16, L = 248, |axis| = 16

[1] THE DIRECTION WEIGHT.  `w = sum_i (|G_i| - 1)`: a group of three directions is worth 2 per
    owner line, a group of two 1, a singleton 0.  Two owner classes sharing a direction `z`
    would meet at `<u,u'> <= 1` everywhere -- the cap constraint `(2/3)<u,u'> + (4/3) <= 2` --
    and so would be ONE class, so the groups are disjoint; inside a group the directions are
    pairwise at `<= -1/2` and `|sum z|^2 >= 0` caps a group at THREE; across groups they are at
    `<= 1/2`.  So the groups' union is a 60-degree code in R^k, `3t + 2q <= tau(k)`, and

        w  <=  floor(2 tau(k)/3)  =  floor(2 x 24/3)  =  16 ,   attained.

    tau(4) = 24 is PROVED, so this ceiling is unconditional -- the only one of the four.

[2] THE AXIS IS MAXIMAL, exhaustively.  A 17th axis point is a UNIT vector `a` of R^4 with
    `<a,a_i> <= 1/2`, `<a,z> <= sqrt3/2` and `|<a,w>| <= 1/sqrt2`.  Every right-hand side is
    positive, so `P = {a : C a <= b}` is a polytope with the origin interior, and `|a|` is
    convex, so `max_P |a|` is attained at a VERTEX.  Enumerate them: the largest is 0.765, so
    P contains no unit vector at all.  CONTROL: drop one of the 16 and the enumeration finds it
    again, and nothing else.

    THE SCOPE is this axis code against this layer frame and these cap directions.
"""
import numpy as np, sys, itertools
from scipy.spatial import HalfspaceIntersection

TAU_KNOWN = {4: 24}          # H. Cohn, Table of kissing number bounds; common/published.py
TAU_PROVED = {4: 24}         # dimension 4 is solved (Musin 2003), so the two agree
FAILS = []


def req(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg, flush=True)
    if not ok:
        FAILS.append(msg)


k = 4
C1 = np.array([v for v in itertools.product((1, 0, -1), repeat=4)
               if sorted(np.abs(v)) == [0, 0, 1, 1]], float)          # the 24-cell, norm^2 = 2
Z = C1 / np.sqrt(2.0)                                                 # cap directions, unit
A = np.array(list(itertools.product((0.5, -0.5), repeat=4)))          # axis: the 16 half-vectors
W = np.eye(4)                                                         # layer frame; +-W is the
L0 = 248                                                              # cross-polytope
req(len(Z) == 24, "the 24-cell has 24 directions")
req(len(A) == 16, "the axis is the 16 half-vectors of the dual 24-cell")

print("== dimension 28: k = 4, 8 classes of %d, %d axis points, %d cap directions ==" %
      (L0, len(A), len(Z)))

# --------------------------------------------------------------------------- [1] the weight
print("\n[1] the direction weight")
tris = [t for t in itertools.combinations(range(24), 3)
        if np.abs(C1[list(t)].sum(0)).max() == 0]
byp = {}
for t in tris:
    for p in t:
        byp.setdefault(p, []).append(t)


def cover(rem, acc):
    if not rem:
        return acc
    p = min(rem)
    for t in byp[p]:
        if set(t) <= rem:
            r = cover(rem - set(t), acc + [t])
            if r is not None:
                return r
    return None


GRP = cover(set(range(24)), [])
req(GRP is not None and len(GRP) == 8, "the 24 directions partition into 8 zero-sum triangles")
COS = Z @ Z.T
req(all(COS[a, b] <= -0.5 + 1e-12 for g in GRP for a, b in itertools.combinations(g, 2)),
    "inside a group the directions are at cosine <= -1/2")
req(all(COS[a, b] <= 0.5 + 1e-12 for i, g in enumerate(GRP) for h in GRP[i+1:]
        for a in g for b in h),
    "across groups they are at cosine <= 1/2, so the union of the groups is a 60-degree code")
t3 = sum(1 for g in GRP if len(g) == 3)
q2 = sum(1 for g in GRP if len(g) == 2)
w = 2*t3 + q2
ck, cp = 2*TAU_KNOWN[k]//3, 2*TAU_PROVED[k]//3
req(w == 16 and w == ck, "weight w = 2t + q = %d IS the ceiling floor(2 tau(4)/3) = %d" % (w, ck))
req(TAU_PROVED[k] == TAU_KNOWN[k],
    "tau(4) = 24 is PROVED, so dimension 28's direction weight is at an UNCONDITIONAL ceiling "
    "-- unlike dimensions 29, 30 and 31, whose tau(k) is only the best known")

# --------------------------------------------------------------------------- [2] the axis
print("\n[2] the axis, exhaustively")
b = np.concatenate([np.full(len(A), 0.5), np.full(len(Z), np.sqrt(3)/2),
                    np.full(2*len(W), 1/np.sqrt(2))])
C = np.vstack([A, Z, W, -W])
req(bool((b > 0).all()), "every constraint has a positive right-hand side, so the origin is "
                         "interior and the polytope is the right object")


def maxnorm(Ax, eps=0.0):
    Cc = np.vstack([Ax, Z, W, -W])
    bb = np.concatenate([np.full(len(Ax), 0.5), np.full(len(Z), np.sqrt(3)/2),
                         np.full(2*len(W), 1/np.sqrt(2))]) * (1.0 + eps)
    V = HalfspaceIntersection(np.hstack([Cc, -bb[:, None]]), np.zeros(k)).intersections
    return float(np.linalg.norm(V, axis=1).max()), V


m0, V = maxnorm(A)
print("   %d vertices;  max |a| over the admissible polytope = %.9f" % (len(V), m0))
req(m0 < 1.0, "no unit vector fits: the axis is MAXIMAL (margin %.1f per cent)" % (100*(1 - m0)))
mp = [maxnorm(A, e)[0] for e in (-1e-6, 1e-6, 1e-4)]
req(max(mp) < 1.0, "and still maximal with the right-hand sides moved by +-1e-6 and +1e-4 (%s), "
                   "so this is not a floating-point decision"
    % ", ".join("%.6f" % x for x in mp))
rng = np.random.default_rng(0)
D = rng.normal(size=(200000, k))
D /= np.linalg.norm(D, axis=1)[:, None]
t = 1.0 / np.maximum((D @ C.T) / b, 1e-300).max(1)
req(t.max() <= m0 + 1e-9, "an independent 200000-direction scan of the same maximum, needing no "
                          "enumeration, reaches %.9f and never exceeds it" % t.max())
extra, rec = 0, []
for j in range(len(A)):
    _, V = maxnorm(np.delete(A, j, 0))
    u = V[np.linalg.norm(V, axis=1) > 1 - 1e-9]
    rec.append(len(u))
    if len(u):
        u = u / np.linalg.norm(u, axis=1)[:, None]
        extra = max(extra, int((np.abs(u @ A[j]) < 1 - 1e-9).sum()))
req(min(rec) >= 1, "CONTROL: removing any one axis point makes a unit vector fit again, so the "
                   "test can fail (recovered %d..%d)" % (min(rec), max(rec)))
req(extra == 0, "and the recovered vector is ALWAYS the point removed -- 0 others over all %d "
                "removals, so no one-for-one swap opens the axis either" % len(A))

# --------------------------------------------------------------------------- [3] what is open
print("\n[3] what is NOT at a ceiling")
print("   the CLASS, and nothing else.  Every class here has %d lines; one more line in the "
      "template is worth +%d in dimension 28.  248 is maximal over all 48576 type-B lines and "
      "rigid under 850000 plateau moves, but the Delsarte LP on the Leech line scheme allows "
      "425.45 lines, and on the 84 octads transverse to one trio -- where every KNOWN class "
      "lives, though no proof says every class must -- the Hoffman ratio bound allows 373.33.  "
      "So it is open." % (L0, 2*w))
print("   The packing is closed (8 disjoint full classes), the layer is 96k = 384 by the frame "
      "argument, the axis is maximal above, and tau(4) = 24 is proved.")

print("\n%s" % ("*** ALL CHECKS PASS ***" if not FAILS else "*** %d CHECKS FAILED ***" % len(FAILS)))
sys.exit(1 if FAILS else 0)
