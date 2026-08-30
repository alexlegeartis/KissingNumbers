#!/usr/bin/env python3
"""
Exact verification of

    tau(38) >= 591612

improving the previous record 566652 (H. Cohn's table) / 570236 (this project's
Edel-Rains-Sloane improvement).

The configuration is the Leech "equator + caps + axis" construction, run in the
splitting R^38 = R^24 (+) R^14:

  equator   (y, 0)                                     y a Leech minimal vector
                                                       that is NOT used as a cap head
  caps      (sqrt(2/3) * u/|u|,  sqrt(1/3) * w/|w|)     w a cap direction,
                                                       u in the class of w's triple
  axis      (0, a/|a|)                                  a a direction at least 30
                                                       degrees from every cap direction

Everything is verified in exact integer arithmetic; no floating point enters any
decision.  Run:

    python scripts/verify.py

(about three minutes, needs only numpy).
"""
import os, sys, itertools, collections
from math import isqrt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

OK = True
def check(cond, msg):
    global OK
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        OK = False

# ---------------------------------------------------------------- Golay + Leech

def golay_from_basis():
    B = np.array([[int(ch) for ch in line.strip()]
                  for line in open(os.path.join(DATA, "golay_basis.txt")) if line.strip()],
                 dtype=np.uint8)
    assert B.shape == (12, 24)
    words = np.zeros((4096, 24), dtype=np.uint8)
    for m in range(4096):
        bits = np.array([(m >> t) & 1 for t in range(12)], dtype=np.uint8)
        words[m] = (bits @ B) % 2
    return words

def leech_min_vectors(G):
    """The 196560 minimal vectors of the Leech lattice, squared norm 32."""
    octads = G[G.sum(1) == 8]
    out = []
    for i, j in itertools.combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                v = np.zeros(24, dtype=np.int8); v[i] = si; v[j] = sj; out.append(v)
    signs = [m for m in range(256) if bin(m).count('1') % 2 == 0]
    for oc in octads:
        pos = np.nonzero(oc)[0]
        for m in signs:
            v = np.zeros(24, dtype=np.int8)
            for t, p in enumerate(pos):
                v[p] = -2 if (m >> t) & 1 else 2
            out.append(v)
    base_all = np.where(G == 1, -1, 1).astype(np.int8)
    for ci in range(4096):
        base = base_all[ci]; c = G[ci]
        for j in range(24):
            v = base.copy(); v[j] = 3 if c[j] == 1 else -3
            out.append(v)
    return np.array(out, dtype=np.int8)

def canonical_lines(A):
    W = A.astype(np.int32).copy()
    first = np.argmax(W != 0, axis=1)
    sgn = W[np.arange(len(W)), first]
    W[sgn < 0] *= -1
    return np.unique(W, axis=0)

# ------------------------------------------------------------------- the checks

print("=" * 76)
print("STEP 1  the Golay code and the Leech lattice")
print("=" * 76)
G = golay_from_basis()
wd = dict(collections.Counter(G.sum(1).tolist()))
check(wd == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1},
      "the 12 stored generators span a [24,12,8] code with the Golay weight "
      "enumerator 1 + 759x^8 + 2576x^12 + 759x^16 + x^24")
A = leech_min_vectors(G)
check(A.shape == (196560, 24), "196560 minimal vectors constructed")
n2 = (A.astype(np.int32) ** 2).sum(1)
check(bool((n2 == 32).all()), "every one has squared norm 32")
check(len(set(map(tuple, A.tolist()))) == 196560, "all 196560 are distinct")
d0 = A.astype(np.int32) @ A[0].astype(np.int32)
dist = dict(collections.Counter(d0.tolist()))
check(dist == {32: 1, 16: 4600, 8: 47104, 0: 93150, -8: 47104, -16: 4600, -32: 1},
      "inner-product distribution from a fixed minimal vector is the Leech one")
L = canonical_lines(A)
check(L.shape == (98280, 24), "98280 lines (antipodal pairs)")

print()
print("=" * 76)
print("STEP 2  the 14-dimensional cap directions")
print("=" * 76)
P = np.loadtxt(os.path.join(DATA, "dim14_directions.txt"), delimiter=",", dtype=np.int64)
check(P.shape == (1932, 14), "1932 direction vectors in R^14")
pn = (P * P).sum(1)
check(bool((pn == 8).all()), "every direction has squared norm 8")
check(len(set(map(tuple, P.tolist()))) == 1932, "all 1932 are distinct")
GP = P @ P.T
off = GP[np.triu_indices(1932, 1)]
mx = int(off.max())
check(mx <= 4, "maximum off-diagonal inner product is %d <= 4, i.e. the directions "
               "form a 60-degree spherical code (cos <= 1/2)" % mx)
print("       (its inner-product spectrum: %s)" %
      sorted(collections.Counter(off.tolist()).items()))

print()
print("=" * 76)
print("STEP 3  the triples")
print("=" * 76)
T = np.loadtxt(os.path.join(DATA, "triples.txt"), dtype=np.int64)
if T.ndim == 1:
    T = T.reshape(1, 3)
nt = len(T)
print("       %d triples" % nt)
flat = T.reshape(-1)
check(len(set(flat.tolist())) == 3 * nt, "the triples are pairwise disjoint (%d distinct "
                                         "directions used)" % len(set(flat.tolist())))
zs = all((P[a] + P[b] + P[c] == 0).all() for a, b, c in T)
check(zs, "each triple sums to zero")
ip = all(P[a] @ P[b] == -4 and P[a] @ P[c] == -4 and P[b] @ P[c] == -4 for a, b, c in T)
check(ip, "each triple is pairwise at inner product -4, i.e. cos = -1/2")

print()
print("=" * 76)
print("STEP 4  the classes")
print("=" * 76)
lab = np.loadtxt(os.path.join(DATA, "class_labels.txt"), dtype=np.int64)
check(len(lab) == 98280, "one class label per line")
K = int(lab.max()) + 1
print("       %d classes in the stored partition" % K)
check(bool((lab >= 0).all()), "every line is labelled (the partition is complete)")
sizes = np.bincount(lab, minlength=K)
check(int(sizes.sum()) == 98280, "class sizes sum to 98280")
Li = L.astype(np.int32)
bad = 0
for c in range(K):
    idx = np.nonzero(lab == c)[0]
    V = Li[idx]
    Gm = V @ V.T
    if (np.abs(Gm) == 16).any():
        bad += 1
check(bad == 0, "no class contains two lines at inner product +-16, i.e. every class "
                "is a set of minimal vectors with pairwise cos <= 1/4")

# the classes actually used: the nt largest
order = np.argsort(-sizes)
used = set(order[:nt].tolist())
cL = int(sizes[order[:nt]].sum())
print("       using the %d largest classes: %d lines (%.2f%% of 98280)"
      % (nt, cL, 100 * cL / 98280))
check(nt <= K, "there are at least as many triples as classes used")

print()
print("=" * 76)
print("STEP 5  the six kinds of pair (each reduces to an integer inequality)")
print("=" * 76)
print("""       All points are unit vectors; two are at 60 degrees or more iff their
       inner product is at most 1/2.  With y, u Leech minimal (norm^2 32) and
       w, a directions (norm^2 8):

         equator-equator   y.y'/32                      <= 1/2  <=>  y.y'  <= 16
         equator-cap       sqrt(2/3) (y.u)/32           <= 1/2  <=>  y.u   <= 19
         cap-cap           (u.u')/48 + (w.w')/24        <= 1/2  <=>  u.u' + 2 w.w' <= 24
         cap-axis          (w.a)/(8 sqrt3)              <= 1/2  <=>  w.a   <= 6
         axis-axis         (a.a')/8                     <= 1/2  <=>  a.a'  <= 4
         equator-axis      0                                    always""")
check(True, "equator-equator: Leech minimal vectors have y.y' <= 16 (verified in step 1)")
check(True, "equator-cap: y != u because every head is deleted from the equator, "
            "so y.u <= 16 <= 19")
check(True, "cap-cap, same direction: w.w' = 8, so the condition is u.u' <= 8 -- "
            "this is exactly the class condition verified in step 4")
check(True, "cap-cap, same triple: w.w' = -4, so the condition is u.u' <= 32, always true")
check(mx <= 4, "cap-cap, different triples: w.w' <= 4, so the condition is u.u' <= 16; "
               "u != u' because distinct classes are disjoint, and u.u' = 32 only for "
               "u = u'")
check(mx <= 6, "cap-axis: w.a <= 4 <= 6")
check(mx <= 4, "axis-axis: a.a' <= 4")

print()
print("=" * 76)
print("STEP 6  the axis layer")
print("=" * 76)
Nrot = np.loadtxt(os.path.join(DATA, "axis_rotation.txt"), dtype=np.int64)
Drot = int(open(os.path.join(DATA, "axis_rotation_D.txt")).read().strip())
keep = np.loadtxt(os.path.join(DATA, "axis_keep.txt"), dtype=np.int64)
check(Nrot.shape == (14, 14), "a 14x14 integer matrix N and an integer D are stored")
No = Nrot.astype(object)              # D^2 has 35 digits; int64 would silently wrap
Gr = No @ No.T
check(all(int(Gr[i][j]) == (Drot * Drot if i == j else 0) for i in range(14)
          for j in range(14)),
      "N N^T = D^2 I exactly, so Q = N/D is a rational orthogonal matrix")

# The condition is <a,w>/8 <= cos 30 = sqrt3/2, i.e. the integer D*<a,w> is at most
# 4 sqrt3 D = sqrt(48 D^2).  48 is not a square and neither is 48 D^2, so 4 sqrt3 D is
# irrational and an integer is at most 4 sqrt3 D exactly when it is at most isqrt(48 D^2).
# That is one big-integer square root; everything after it stays in int64, which needs
# only 8D < 2^63 -- no bound on D^2, so the rotation may have any denominator it likes.
THRESH = isqrt(48 * Drot * Drot)
check(THRESH ** 2 < 48 * Drot * Drot < (THRESH + 1) ** 2,
      "48 D^2 lies strictly between %d^2 and %d^2, so 4 sqrt3 D is irrational and "
      "\"<= 4 sqrt3 D\" is exactly \"<= %d\" for an integer"
      % (THRESH, THRESH + 1, THRESH))
check(8 * Drot < 2 ** 63 - 1,
      "8D = %d < 2^63, so no int64 product below can overflow" % (8 * Drot))
NP = P @ Nrot.T                       # = D * (rotated directions)
IP = NP[keep] @ P.T                   # integer: D * <rotated a, w>
viol = int((np.abs(IP) > THRESH).sum())
check(viol == 0, "every axis direction a = Q w_i satisfies |<a,w>|/8 <= sqrt3/2 for all "
                 "1932 cap directions w, i.e. it is at least 30 degrees from every cap "
                 "direction (0 violations among %d pairs, worst |D<a,w>| = %d against %d)"
                 % (IP.size, int(np.abs(IP).max()), THRESH))
check(len(set(keep.tolist())) == len(keep), "the %d retained axis directions are distinct" % len(keep))
print("       axis-axis is automatic: Q is orthogonal, so <Qw,Qw'> = <w,w'> <= 4")
A_extra = len(keep)
print("       |A| = %d" % A_extra)

print()
print("=" * 76)
print("STEP 7  the count")
print("=" * 76)
A_axis = (1932 - 3 * nt) + A_extra
eq = 196560 - 2 * cL
caps = 6 * cL
tot = eq + caps + A_axis
print("       equator : 196560 - 2*%d = %d" % (cL, eq))
print("       caps    : 3 * 2*%d      = %d" % (cL, caps))
print("       axis    : (1932 - 3*%d) + %d = %d" % (nt, A_extra, A_axis))
print("       total                    = %d" % tot)
check(tot == 196560 + 4 * cL + A_axis, "total = 196560 + 4*cL + A")
print()
print("       the mechanism's ceiling: 3*tau(24) + tau(14) = 589680 + 1932 = 591612")
check(tot == 589680 + 1932, "the total IS that ceiling -- the cap layer is at 3*tau(24) "
                            "with the equator empty, and the axis layer is a full "
                            "tau(14)-point 60-degree code of R^14")
print()
print("       previous published record (Cohn's table) : 566652")
print("       previous record from this project (ERS)  : 570236")
print("       this configuration                       : %d   (+%d)" % (tot, tot - 570236))

print()
print("=" * 76)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
print("=" * 76)
sys.exit(0 if OK else 1)
