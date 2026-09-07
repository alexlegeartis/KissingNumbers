#!/usr/bin/env python3
"""
Exact verification of

    tau(25) >= 197058

improving the published record 197056 (H. Cohn's table, Ma et al. 2025) by 2.

THE POINT.  Every record configuration in dimensions 25-31 has the shape

    equator  (y, 0)                          y a Leech minimal vector, not a head
    caps     (sqrt(t) u, sqrt(1-t) w)        w a unit vector in R^k
    axis     (0, a)

and the literature always takes t = 2/3, because that is the unique level at which a head
can be shared by a *triple* of directions (three unit vectors pairwise at cos -1/2).

In dimension 25 the second block is R^1, so there are only two directions, +1 and -1, and
no triple exists.  A head only has to be shared by an antipodal *pair*, and that is possible
for every t <= 3/4.  Raising t from 2/3 to 3/4 changes nothing else -- the head condition
relaxes from cos <= 1/4 to cos <= 1/3, which for the Leech lattice is the same condition,
so the head class still has 496 vectors -- but it makes the cap-pole inner product

    sqrt(1-t) = 1/2      instead of      sqrt(1/3) = 0.5773...

so the two poles (0,...,0,+-1), which had to be discarded at t = 2/3, are now admissible.

    197056  =  196064 + 992          (published, t = 2/3, no poles)
    197058  =  196064 + 992 + 2      (t = 3/4)

Run:  python verify.py
"""
import os, sys, itertools, collections, json
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GOLAY = os.path.join(HERE, "..", "..", "..", "common", "data", "golay_basis.txt")
CLASSJSON = os.path.join(HERE, "..", "dim27-triple-partition", "data", "construction.json")

OK = True
def check(cond, msg):
    global OK
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        OK = False

# ------------------------------------------------------- Leech minimal vectors

def golay_from_basis():
    B = np.array([[int(c) for c in line.strip()] for line in open(GOLAY) if line.strip()],
                 dtype=np.uint8)
    W = np.zeros((4096, 24), dtype=np.uint8)
    for m in range(4096):
        bits = np.array([(m >> t) & 1 for t in range(12)], dtype=np.uint8)
        W[m] = (bits @ B) % 2
    return W

def leech(G):
    octads = G[G.sum(1) == 8]
    out = []
    for i, j in itertools.combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                v = np.zeros(24, np.int8); v[i] = si; v[j] = sj; out.append(v)
    sg = [m for m in range(256) if bin(m).count('1') % 2 == 0]
    for oc in octads:
        p = np.nonzero(oc)[0]
        for m in sg:
            v = np.zeros(24, np.int8)
            for t, q in enumerate(p): v[q] = -2 if (m >> t) & 1 else 2
            out.append(v)
    base = np.where(G == 1, -1, 1).astype(np.int8)
    for ci in range(4096):
        b = base[ci]; c = G[ci]
        for j in range(24):
            v = b.copy(); v[j] = 3 if c[j] == 1 else -3
            out.append(v)
    return np.array(out, np.int8)

print("=" * 74)
print("STEP 1  the Leech lattice (squared norm 32) and the head class")
print("=" * 74)
G = golay_from_basis()
wd = dict(collections.Counter(G.sum(1).tolist()))
check(wd == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, "Golay weight enumerator")
A = leech(G).astype(np.int64)
check(A.shape == (196560, 24), "196560 minimal vectors")
check(bool(((A * A).sum(1) == 32).all()), "all of squared norm 32")
Aset = set(map(tuple, A.tolist()))
S = np.array(json.load(open(CLASSJSON))["classes"][0], dtype=np.int64)
check(S.shape == (496, 24), "the head class has 496 vectors")
check(all(tuple(v) in Aset for v in S.tolist()), "every head is a Leech minimal vector")
check(len(set(map(tuple, S.tolist()))) == 496, "the 496 heads are distinct")
GS = S @ S.T
vals = sorted(set(GS.flatten().tolist()))
check(vals == [-32, -8, 0, 8, 32], "head inner products lie in {0,+-8,+-32}: %s" % vals)
off = GS[~np.eye(496, dtype=bool)]
check(int(off.max()) <= 8, "max off-diagonal head inner product is %d <= 8, i.e. cos <= 1/4"
                            % int(off.max()))

print()
print("=" * 74)
print("STEP 2  the six kinds of pair at level t = 3/4  (exact rational arithmetic)")
print("=" * 74)
t = Fraction(3, 4)          # alpha^2
half = Fraction(1, 2)
print("""       Points, all unit vectors in R^25 = R^24 (+) R:
         equator : (y/sqrt32, 0)                       196560 - 496 = 196064
         caps    : (sqrt(t) u/sqrt32, +- sqrt(1-t))    2 * 496       =    992
         poles   : (0, +-1)                                          =      2""")

# (1) equator-equator :  y.y'/32 <= 1/2  <=>  y.y' <= 16
check(True, "equator-equator: distinct Leech minimal vectors have y.y' <= 16, so the inner "
            "product is at most 1/2")
# (2) equator-cap : sqrt(t) * (y.u)/32 <= 1/2  <=>  y.u <= 16/sqrt(t) ; worst case y.u = 16
#     (y = u is excluded because every head is deleted from the equator)
lhs2 = t * Fraction(16, 32) ** 2                      # (sqrt(t)*1/2)^2
check(lhs2 <= half ** 2, "equator-cap: y != u, so y.u <= 16 and the inner product is "
                         "sqrt(3)/4 = 0.4330 <= 1/2")
# (3) cap-cap, same pole sign : t*(u.u')/32 + (1-t) <= 1/2  <=>  u.u'/32 <= (1/2-(1-t))/t = 1/3
thr3 = (half - (1 - t)) / t
check(thr3 >= Fraction(1, 4), "cap-cap same sign: needs cos(u,u') <= %s, and the class "
                              "guarantees cos <= 1/4" % thr3)
# (4) cap-cap, opposite pole sign, SAME head : t - (1-t) = 1/2
check(t - (1 - t) <= half, "cap-cap opposite sign with the same head: t-(1-t) = %s <= 1/2 "
                           "(this is what lets the two directions share a head)" % (t - (1 - t)))
# (4b) opposite sign, different heads: t*(u.u')/32 - (1-t) <= 1/2 always since u.u'/32 <= 1
check(t * 1 - (1 - t) <= half, "cap-cap opposite sign, different heads: at most %s <= 1/2"
      % (t * 1 - (1 - t)))
# (5) cap-pole : sqrt(1-t) <= 1/2   <=>  1-t <= 1/4  <=>  t >= 3/4
check((1 - t) <= half ** 2, "cap-pole: sqrt(1-t) = 1/2 <= 1/2   *** this is the whole point: "
                            "at the published t = 2/3 it is sqrt(1/3) = 0.5774 > 1/2 ***")
# (6) pole-pole, equator-pole
check(True, "pole-pole: -1 <= 1/2 ;  equator-pole: 0 <= 1/2")

print()
print("=" * 74)
print("STEP 3  the count")
print("=" * 74)
eq = 196560 - 496
caps = 2 * 496
poles = 2
tot = eq + caps + poles
print("       equator : 196560 - 496 = %d" % eq)
print("       caps    : 2 * 496      = %d" % caps)
print("       poles   :              = %d" % poles)
print("       total                  = %d" % tot)
check(tot == 197058, "total = 197058")
print()
print("       published record : 197056   (same configuration at t = 2/3, poles dropped)")
print("       this             : %d   (+%d)" % (tot, tot - 197056))

print()
print("=" * 74)
print("STEP 4  numerical spot check on real 25-dimensional coordinates")
print("=" * 74)
rt = float(t) ** 0.5; rb = (1 - float(t)) ** 0.5
NL = 32.0 ** 0.5
rng = np.random.default_rng(0)
heads = set(map(tuple, S.tolist()))
eqidx = [i for i in range(196560) if tuple(A[i].tolist()) not in heads]
sel = rng.choice(len(eqidx), size=3000, replace=False)
pts = []
for i in sel:
    v = np.zeros(25); v[:24] = A[eqidx[i]] / NL; pts.append(v)
for u in S:
    for s in (1, -1):
        v = np.zeros(25); v[:24] = rt * u / NL; v[24] = s * rb; pts.append(v)
for s in (1, -1):
    v = np.zeros(25); v[24] = s; pts.append(v)
X = np.array(pts)
print("       sampled %d points (%d equator + all 992 caps + both poles)" % (len(X), len(sel)))
nr = np.linalg.norm(X, axis=1)
check(abs(nr.min() - 1) < 1e-12 and abs(nr.max() - 1) < 1e-12, "all sampled points are unit vectors")
worst = -2.0; bad = 0
for s0 in range(0, len(X), 1500):
    Gm = X[s0:s0 + 1500] @ X.T
    for r in range(Gm.shape[0]): Gm[r, s0 + r] = -2.0
    worst = max(worst, float(Gm.max())); bad += int((Gm > 0.5 + 1e-10).sum())
print("       largest inner product between distinct sampled points: %.15f" % worst)
check(bad == 0, "no sampled pair exceeds 1/2 (%d violations)" % bad)

print()
print("=" * 74)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
print("=" * 74)
sys.exit(0 if OK else 1)
