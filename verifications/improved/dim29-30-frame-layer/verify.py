#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""verify.py --- exact verification of K(29) >= 209968 and K(30) >= 220948

    python verify.py 29        python verify.py 30        (add --fast to skip the two sweeps)

Reads only `data/` and `lib/golay.py`.  The Leech shell is rebuilt from the Golay code and its
inner-product signature is checked before anything is trusted, so nothing depends on the code
that produced the configuration.

THE CONFIGURATION.  Norm-4 units: a point of R^{24+k} is (x, y) with |x|^2 + |y|^2 = 4, and two
points are compatible iff their inner product is at most 2.

    equator  (u, 0)                        u a minimal vector that is not an owner
    cap      (sqrt(2/3) u, (2/sqrt3) z)    u an owner, z a direction of its group
    axis     (0, 2a)                       a a unit vector of R^k
    layer    (v/2, sqrt2 w)                v one of the 48 vectors of a Leech FRAME,
                                           w one of the 2k directions of a cross-polytope

The layer is the new part.  A frame vector has |v|^2 = 8, so |v/2|^2 = 2 and the layer sits at
height sqrt2.  Three of its four constraints are free:

  * it deletes no equator point, because a norm-8 lattice vector meets a minimal one at
    <v,u> <= 4, hence <v/2, u> <= 2 exactly;
  * one head carries every direction with <w,w'> <= 0, so a whole cross-polytope, 2k of them;
  * two heads share a direction as soon as <v,v'> <= 0, so all 48 vectors of a frame carry all
    2k directions: 96k layer points.

The fourth is the whole content.  Against a cap the constraint is

    <u,v>/sqrt6 + (2 sqrt2/sqrt3) <z,w>  <=  2,

so if every owner is TYPE B for the frame -- |<u,v_i>| = 2 for all 24 heads, which for the frame
8e_i means the owner has shape (+-2^8, 0^16) -- the threshold on the directions is

    |<z,w>|  <=  sqrt(3/2) - 1/2  =  0.724744871...,

and the tau(k) cap directions are a root system, which has an orthonormal frame in which every
root has maximum coordinate 1/sqrt2 = 0.707106781...  The margin is the single inequality

    1/sqrt2 <= sqrt(3/2) - 1/2   <=>   1 + 1/sqrt2 <= sqrt3   <=>   sqrt2 < 3/2.

Exit code 0 iff every check passes.
"""
import numpy as np, sys, os, json, itertools, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from golay import golay24

DIM = int(sys.argv[1]) if len(sys.argv) > 1 else 29
if DIM not in (29, 30):
    sys.exit("usage: python verify.py 29|30 [--fast]")
FULL = '--fast' not in sys.argv

FAILS = []
T0 = time.time()


def req(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg, flush=True)
    if not ok:
        FAILS.append(msg)


# --------------------------------------------------------------------------------------------
# exact arithmetic in Z[sqrt2, sqrt3]: a quadruple (a,b,c,d) is (a + b r2 + c r3 + d r6)/DEN
# --------------------------------------------------------------------------------------------
def qmul(x, y):
    a, b, c, d = x
    p, q, r, s = y
    return (a*p + 2*b*q + 3*c*r + 6*d*s,
            a*q + b*p + 3*c*s + 3*d*r,
            a*r + c*p + 2*b*s + 2*d*q,
            a*s + d*p + b*r + c*q)


def qdot(u, v):
    t = (0, 0, 0, 0)
    for x, y in zip(u, v):
        m = qmul(x, y)
        t = (t[0]+m[0], t[1]+m[1], t[2]+m[2], t[3]+m[3])
    return t


def qsign(t):
    """sign of a + b sqrt2 + c sqrt3 + d sqrt6, exactly"""
    import sympy as sp
    a, b, c, d = (int(v) for v in t)
    return int(sp.sign(a + b*sp.sqrt(2) + c*sp.sqrt(3) + d*sp.sqrt(6)))


def qle(t, num, den, scale):
    """decide  (a + b r2 + c r3 + d r6)/scale  <=  num/den   exactly"""
    a, b, c, d = (int(v) for v in t)
    return qsign((a*den - num*scale, b*den, c*den, d*den)) <= 0


def qval(t, scale):
    a, b, c, d = (float(v) for v in t)
    return (a + b*np.sqrt(2) + c*np.sqrt(3) + d*np.sqrt(6)) / scale


# --------------------------------------------------------------------------------------------
# the Leech shell, rebuilt from the Golay code, at minimal norm 32
# --------------------------------------------------------------------------------------------
print("== the Leech shell, from lib/golay.py ==")
G = golay24().astype(np.int8)
OCT = G[G.sum(1) == 8]
req(len(G) == 4096 and len(OCT) == 759, "Golay code: 4096 words, 759 octads")

rows = []
for i, j in itertools.combinations(range(24), 2):            # (+-4^2, 0^22)
    for si in (4, -4):
        for sj in (4, -4):
            v = np.zeros(24, np.int8); v[i] = si; v[j] = sj; rows.append(v)
ev = np.array([s for s in itertools.product((1, -1), repeat=8) if np.prod(s) == 1], np.int8)
for o in OCT:                                                # (+-2^8, 0^16) on an octad
    p = np.flatnonzero(o)
    for s in ev:
        v = np.zeros(24, np.int8); v[p] = 2*s; rows.append(v)
V = 1 - 2*G.astype(np.int8)
for j in range(24):                                          # (-+3, +-1^23)
    X = V.copy(); X[:, j] = -3*V[:, j]; rows.append(X)
LS = np.vstack([r if r.ndim == 2 else r[None, :] for r in rows]).astype(np.int8)
req(len(LS) == 196560, "196560 vectors built (%d)" % len(LS))
N = (LS.astype(np.int64)**2).sum(1)
req((N == 32).all(), "every vector has norm 32 (= 4 in norm-4 units)")
req(np.abs(LS).max() == 4, "no coordinate exceeds 4, so <u, 8e_i>/8 <= 4 for every head")


def canon(X):
    X = np.asarray(X, np.int64)
    f = (X != 0).argmax(1)
    return X * np.sign(X[np.arange(len(X)), f])[:, None]


def key(X):
    """one bytes object per line, sign-canonicalised; coordinates fit in an int8"""
    return [bytes(r) for r in np.ascontiguousarray(canon(X).astype(np.int8))]


LN = np.unique(canon(LS), axis=0)
req(len(LN) == 98280, "98280 distinct lines (%d)" % len(LN))
KEY = {k: i for i, k in enumerate(key(LN))}

if FULL:
    print("   inner-product signature of the 98280 lines ...", flush=True)
    # float32 is exact here: coordinates are at most 4 and a dot product at most 32
    prof = np.zeros(5, np.int64)
    L32 = LN.astype(np.float32)
    for i in range(0, 98280, 4096):
        B = np.abs(L32[i:i+4096] @ L32.T)
        for t in range(5):
            prof[t] += int((B == 8*t).sum())
    req(prof.sum() == 98280*98280,
        "every inner product is 0, 8, 16, 24 or 32 -- i.e. 0, 1, 2, 3 or 4 in norm-4 units")
    prof[4] -= 98280                                          # each line against itself
    per = prof // 98280
    req(prof.sum() == 98280*98280 - 98280, "the profile covers every ordered pair of lines")
    req(list(per) == [46575, 47104, 4600, 0, 0],
        "signature 46575 / 47104 / 4600 / 0 / 0 per line at |<u,u'>| = 0, 1, 2, 3, 4 (got %s)"
        % list(per))
    print("   -> distinct minimal lines meet at |<u,u'>| <= 2, which settles equator x equator")

FR = 8*np.eye(24, dtype=np.int64)
req(((FR @ FR.T) == 64*np.eye(24)).all(), "the 24 heads 8e_i are orthogonal of norm 64 (= 8)")


def in_leech(x):
    x = np.asarray(x, np.int64)
    m = int(x[0] % 2)
    if not ((x % 2) == m).all():
        return False
    S = ((x % 4) == (2 if m == 0 else 1)).astype(np.int8)
    if not any(bool((S == c).all()) for c in G):
        return False
    return int(x.sum()) % 8 == (4*m) % 8


req(all(in_leech(f) for f in FR), "the heads are Leech vectors (Conway-Sloane criterion)")
req(all(in_leech(LN[i]) for i in (0, 1000, 50000, 98279)), "sampled lines pass the same criterion")
req(not in_leech(LN[0] + np.eye(24, dtype=np.int64)[0]),
    "   (negative control: the criterion rejects a perturbed vector)")

# --------------------------------------------------------------------------------------------
# the owner classes
# --------------------------------------------------------------------------------------------
print("\n== owner classes ==")
OWN = np.load(os.path.join(HERE, 'data', 'owners%d.npy' % DIM)).astype(np.int64)
BND = np.load(os.path.join(HERE, 'data', 'bounds%d.npy' % DIM))
CLS = [OWN[BND[i]:BND[i+1]] for i in range(len(BND)-1)]
req(all(k in KEY for k in key(OWN)), "every owner line is a Leech minimal line")
req(len(set(key(OWN))) == len(OWN), "the classes are pairwise disjoint")
req(np.abs(OWN).max() == 2, "every owner is TYPE B for this frame: |<u, 8e_i>|/8 = 2, never 3 or 4")
worst = max(int(np.abs(c @ c.T - 32*np.eye(len(c), dtype=np.int64)).max()) for c in CLS)
req(worst <= 8, "every class is a 75.52-degree code: |<u,u'>| <= 1 in norm-4 units (max %g)"
    % (worst/8.0))
print("   %d classes, sizes %d..%d, %d owner lines in total"
      % (len(CLS), min(map(len, CLS)), max(map(len, CLS)), len(OWN)))

# --------------------------------------------------------------------------------------------
# directions, axis, frame -- all exact
# --------------------------------------------------------------------------------------------
print("\n== directions, axis and layer, exactly ==")
M = json.load(open(os.path.join(HERE, 'data', 'geom%d.json' % DIM)))
k, DEN = M['k'], M['denominator']
Z, W, A, GRP = M['directions'], M['frame'], M['axis'], M['groups']
req(k == DIM - 24, "k = %d" % k)
req(len(GRP) == len(CLS), "one direction group per class (%d)" % len(GRP))
# A group of three directions is worth 2 per owner vector and a group of two only 1, so the
# pairing must give the larger classes the larger groups.  Both lists are stored largest first.
req(all(len(CLS[i]) >= len(CLS[i+1]) for i in range(len(CLS)-1))
    and all(len(GRP[i]) >= len(GRP[i+1]) for i in range(len(GRP)-1)),
    "classes and groups are both stored largest first, so the pairing is the best one")
flat = [i for g in GRP for i in g]
req(len(set(flat)) == len(flat) <= len(Z), "the groups use %d distinct directions of %d"
    % (len(set(flat)), len(Z)))

S2 = DEN*DEN
req(all(qdot(z, z) == (S2, 0, 0, 0) for z in Z), "cap directions are unit vectors")
req(all(qdot(a, a) == (S2, 0, 0, 0) for a in A), "axis points are unit vectors")
req(all(qdot(w, w) == (S2, 0, 0, 0) for w in W)
    and all(qdot(W[i], W[j]) == (0, 0, 0, 0) for i in range(k) for j in range(i+1, k)),
    "the layer frame is orthonormal")

req(all(qle(qdot(Z[a], Z[b]), -1, 2, S2) for g in GRP for a, b in itertools.combinations(g, 2)),
    "inside a group the directions are at cosine <= -1/2 (so a cap carries them all)")
req(all(qle(qdot(Z[a], Z[b]), 1, 2, S2)
        for i, g in enumerate(GRP) for h in GRP[i+1:] for a in g for b in h),
    "across groups the directions are at cosine <= 1/2")
# cap x layer:  <u,v>/sqrt6 + (2 sqrt2/sqrt3) <z,w> <= 2 with <u,v> = 2 gives
# <z,w> <= sqrt(3/2) - 1/2 = (sqrt6 - 1)/2 = (-12 + 12 sqrt6)/DEN.  qdot carries <z,w>*DEN^2, so the
# comparison <z,w> <= THR/DEN is  qdot*DEN - THR*DEN^2 <= 0.
THR = (-12, 0, 0, 12)
RHS = tuple(S2*t for t in THR)
ok = all(qsign(tuple(s*x*DEN - y for x, y in zip(qdot(Z[i], W[j]), RHS))) <= 0
         for i in set(flat) for j in range(k) for s in (1, -1))
mx = max(abs(qval(qdot(Z[i], W[j]), S2)) for i in set(flat) for j in range(k))
req(ok, "cap x layer: |<z,w>| <= sqrt(3/2) - 1/2 for every direction and frame vector "
        "(max %.9f vs %.9f)" % (mx, np.sqrt(1.5) - 0.5))
req(all(qle(qdot(A[i], A[j]), 1, 2, S2) for i in range(len(A)) for j in range(len(A)) if i != j),
    "axis x axis: the axis is a 60-degree code")
req(all(qle(qdot(a, Z[i]), 0, 1, 1) or qsign(tuple(4*x - y for x, y in
        zip(qmul(qdot(a, Z[i]), qdot(a, Z[i])), (3*S2*S2, 0, 0, 0)))) <= 0
        for a in A for i in set(flat)),
    "axis x caps: <a,z> <= sqrt3/2")
req(all(qsign(tuple(2*x - y for x, y in zip(qmul(qdot(a, w), qdot(a, w)), (S2*S2, 0, 0, 0)))) <= 0
        for a in A for w in W),
    "axis x layer: |<a,w>| <= 1/sqrt2")
req(qsign((-3, 2, 0, 0)) < 0, "cap x layer binds at (sqrt6 + 2 sqrt3)/3 < 2, i.e. sqrt2 < 3/2")

# --------------------------------------------------------------------------------------------
# the count
# --------------------------------------------------------------------------------------------
print("\n== the count ==")
own_vec = 2*len(OWN)
caps = sum(2*len(c)*len(g) for c, g in zip(CLS, GRP))
layer = 96*k
total = 196560 - own_vec + caps + len(A) + layer
print("   equator %6d + caps %6d + axis %3d + layer %3d = %d"
      % (196560 - own_vec, caps, len(A), layer, total))
req(total == M['total'], "the count agrees with the shipped total (%d)" % M['total'])

# --------------------------------------------------------------------------------------------
# an independent numerical sweep: every pair except equator x equator
# --------------------------------------------------------------------------------------------
if FULL:
    print("\n== independent sweep of every pair except equator x equator ==")
    r8 = np.sqrt(8.0)
    Zf = np.array([[qval(t, DEN) for t in z] for z in Z])
    Af = np.array([[qval(t, DEN) for t in a] for a in A])
    Wf = np.array([[qval(t, DEN) for t in w] for w in W])
    ownk = set(key(OWN))          # key() canonicalises the sign, so this drops both +-u
    keep = np.array([kk not in ownk for kk in key(LS)])
    EQ = np.hstack([LS[keep].astype(np.float64)/r8, np.zeros((int(keep.sum()), k))])
    cap = []
    for c, g in zip(CLS, GRP):
        U = np.vstack([c, -c]).astype(np.float64)/r8
        for zi in g:
            cap.append(np.hstack([np.sqrt(2/3)*U, np.tile((2/np.sqrt(3))*Zf[zi], (len(U), 1))]))
    cap = np.vstack(cap)
    ax = np.hstack([np.zeros((len(Af), 24)), 2*Af])
    Wpm = np.vstack([Wf, -Wf])
    I24 = np.eye(24)
    lay = np.array([np.hstack([np.sqrt(2)*s*I24[i], np.sqrt(2)*w])
                    for i in range(24) for s in (1, -1) for w in Wpm])
    pts = np.vstack([EQ, cap, ax, lay])
    req(len(pts) == total, "the assembled point set has %d points" % total)
    req(np.abs((pts*pts).sum(1) - 4).max() < 1e-9, "every point has squared length 4")
    req(len(np.unique(np.round(pts, 8), axis=0)) == len(pts), "all points are distinct")
    rest = np.vstack([cap, ax, lay])
    worstf = -9.0
    for i in range(0, len(rest), 4000):
        B = rest[i:i+4000]
        Gm = B @ rest.T
        Gm[np.arange(len(B)), np.arange(i, i+len(B))] = -9
        worstf = max(worstf, Gm.max(), (B @ EQ.T).max())
    print("   worst inner product found: %.12f" % worstf)
    req(worstf <= 2 + 1e-9, "every pair outside equator x equator is at inner product <= 2")

print("\n%s   K(%d) >= %d      (%.0f s)"
      % ("*** ALL CHECKS PASS ***" if not FAILS else "*** %d CHECKS FAILED ***" % len(FAILS),
         DIM, total, time.time() - T0))
sys.exit(1 if FAILS else 0)
