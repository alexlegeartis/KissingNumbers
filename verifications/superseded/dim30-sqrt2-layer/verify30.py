# -*- coding: utf-8 -*-
"""Exact verification of  K(30) >= 220440 + (layer)   (python verify30.py)

Reads data/classes.npy (24 x 496 x 24, Cohn coordinates, minimal vectors of norm 32: the owner
classes, one per zero-sum triangle of the cap E6), data/heads.npy (N x 24, Leech vectors of norm
48: the heads of the new layer) and data/directions.json (unit vectors of R^6 as sympy strings:
the 24 triangles, the 72 axis directions, and the three E6* directions of the layer), rebuilds
the Leech minimal vectors from the Golay code, and checks every pair of points exactly: integer
inner products for the x-parts, sympy in Q(sqrt2, sqrt3) for the y-parts.

In norm-4 units the configuration is
    equator  (u/sqrt8, 0)                u a minimal vector that is not an owner
    caps     (u/sqrt12, (2/sqrt3) z)     u in class i, z in triangle i          (3 per owner)
    axis     (0, 2a)                     a one of the 72 axis directions
    layer    (+-f/sqrt24, sqrt2 w_j)     f a head line, w_j one of its own directions
and two points are compatible iff their inner product is at most 2.
"""
from __future__ import print_function
import os, sys, json, time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build, in_lattice        # noqa: E402

t0 = time.time()


def log(s):
    print('  ' + s); sys.stdout.flush()


def ok(expr):
    """expr <= 0, decided exactly (zero is allowed: touching spheres)"""
    g = sp.simplify(expr)
    if g == 0:
        return True
    return bool(sp.N(g, 60) < 0)


# ---------------------------------------------------------------- data
C = np.load(os.path.join(HERE, 'data', 'classes.npy')).astype(np.int64)
F = np.load(os.path.join(HERE, 'data', 'heads.npy')).astype(np.int64)
D = json.load(open(os.path.join(HERE, 'data', 'directions.json')))
vec = lambda ss: sp.Matrix([sp.sympify(s) for s in ss])
TRI = [[vec(z) for z in T] for T in D['triangles']]
AX = [vec(a) for a in D['axis']]
WS = [vec(w) for w in D['w']]
K = len(WS)
HD = [list(h) for h in D['head_dirs']]                # the directions carried by each head line
N = len(F)
assert C.shape == (24, 496, 24) and F.shape == (N, 24) and len(TRI) == 24 and len(AX) == 72
VAL = (C != 0).any(2)                                 # a zero row is an owner given up to another class
NOWN = int(VAL.sum())
assert len(HD) == N and all(0 <= j < K for hd in HD for j in hd)
LAYER = sum(2 * len(hd) for hd in HD)
log('data: 24 classes, %d owners in all (%d of the 24 x 496 slots are empty), %d head lines (both signs used), %d layer directions, %d layer points'
    % (NOWN, 24 * 496 - NOWN, N, K, LAYER))

# ---------------------------------------------------------------- 0. directions
Z = [z for T in TRI for z in T]
for T in TRI:
    assert len(T) == 3 and all(sp.simplify(z.dot(z) - 1) == 0 for z in T), 'cap directions are unit'
    assert all(sp.simplify(c) == 0 for c in (T[0] + T[1] + T[2])), 'each triangle sums to zero'
assert len({tuple(sp.nsimplify(c) for c in z) for z in Z}) == 72, '72 distinct cap directions'
assert all(sp.simplify(a.dot(a) - 1) == 0 for a in AX) and len({tuple(sp.nsimplify(c) for c in a) for a in AX}) == 72
assert all(sp.simplify(w.dot(w) - 1) == 0 for w in WS), 'layer directions are unit'
cosWW = [[sp.nsimplify(sp.simplify(WS[a].dot(WS[b]))) for b in range(K)] for a in range(K)]
assert all(cosWW[a][b] <= 0 for hd in HD for a in hd for b in hd if a != b), (
    'the directions carried by one head are pairwise non-acute')
cosZZ = [[sp.nsimplify(sp.simplify(Z[p].dot(Z[q]))) for q in range(72)] for p in range(72)]
assert all(cosZZ[p][q] in (0, sp.Rational(1, 2), -sp.Rational(1, 2), 1, -1) for p in range(72) for q in range(72)), \
    'cap-direction cosines are 0, +-1/2, +-1 (an E6 root system)'
assert all(cosZZ[p][q] <= sp.Rational(1, 2) for p in range(72) for q in range(72) if p != q), 'cap directions pairwise at >= 60 degrees'
log('directions: 24 zero-sum unit triangles, cosines in {0, +-1/2, +-1}; 72 unit axis directions; %d unit layer directions; on one head the cosines are %s'
    % (K, sorted({str(cosWW[a][b]) for hd in HD for a in hd for b in hd if a != b})))

# ---------------------------------------------------------------- 1. classes
M = leech_build().astype(np.int64)
pw = 9 ** np.arange(24, dtype=np.int64)
key = lambda V: ((V + 4) * pw).sum(-1)
SK = np.sort(key(M))
kC = key(C)[VAL]
pos = np.searchsorted(SK, kC)
assert (SK[np.minimum(pos, len(SK) - 1)] == kC).all(), 'every class vector is a Leech minimal vector'
assert len(np.unique(kC)) == NOWN, 'the classes are pairwise disjoint'
assert ((C ** 2).sum(2)[VAL] == 32).all(), 'every owner has norm 32'
Mij = np.zeros((24, 24), dtype=np.int64)         # max inner product between distinct vectors of classes i and j
for i in range(24):
    for j in range(24):
        G = C[i] @ C[j].T
        if i == j:
            np.fill_diagonal(G, -10 ** 6)
        G = np.where(VAL[i][:, None] & VAL[j][None, :], G, -10 ** 6)
        Mij[i, j] = G.max()
assert (np.diag(Mij) <= 8).all(), 'inside a class the inner product is at most 8 (angle >= arccos(1/4) = 75.52 deg)'
log('classes: %d distinct minimal vectors; inside a class the largest inner product between distinct vectors is %d (bound 8: the owners of a class share all three cap directions, so 60-degree-freeness is not enough)' % (NOWN, np.diag(Mij).max()))
owners = C[VAL]
kO = key(owners)

# ---------------------------------------------------------------- 2. heads
assert ((F ** 2).sum(1) == 48).all(), 'heads have norm 48'
assert all(in_lattice(f) for f in F), 'heads are Leech vectors'
GF = F @ F.T
AGF = np.abs(GF)
np.fill_diagonal(AGF, 0)
assert (AGF < 48).all(), 'the head lines are distinct and none is the negative of another'
log("heads: %d Leech vectors of norm 48; largest |<f,f'>| between distinct head lines is %d" % (N, AGF.max() if N > 1 else 0))

# ---------------------------------------------------------------- 3. cap vs cap
# (u/sqrt12, (2/sqrt3) z).(u'/sqrt12, (2/sqrt3) z') = <u,u'>/12 + (4/3) cos <= 2   <=>   <u,u'> <= 24 - 16 cos
thr = {c: int(sp.nsimplify(24 - 16 * c)) for c in (0, sp.Rational(1, 2), -sp.Rational(1, 2), 1, -1)}
tight = 0
for i in range(24):
    for j in range(24):
        for p in range(3):
            for q in range(3):
                c = cosZZ[3 * i + p][3 * j + q]
                m = Mij[i, j]
                if i == j and p != q:
                    m = max(m, 32)                     # the same owner on two directions of its triangle
                assert m <= thr[c], 'cap-cap (classes %d,%d directions %d,%d)' % (i, j, p, q)
                tight += int(m == thr[c])
log('cap-cap: all class pairs and direction pairs admissible (thresholds %s); %d tight cases' % (thr, tight))

# ---------------------------------------------------------------- 4. cap vs axis, axis vs axis
h = 2 / sp.sqrt(3)
worst = None
for z in Z:
    for a in AX:
        v = sp.simplify(2 * h * z.dot(a))
        assert ok(v - 2), 'cap-axis'
        if worst is None or sp.N(v, 30) > sp.N(worst, 30):
            worst = v
log('cap-axis: maximum inner product %s = %.6f (bound 2)' % (worst, float(worst)))
for i, a in enumerate(AX):
    for b in AX[i + 1:]:
        assert ok(4 * a.dot(b) - 2), 'axis-axis'
log('axis-axis: 72 points of norm 4, pairwise inner product <= 2')

# ---------------------------------------------------------------- 5. equator vs caps and vs layer
mx_eq_owner = -10 ** 6
mx_eq_head = -10 ** 6
B = 4096
for s in range(0, len(M), B):
    blk = M[s:s + B]
    G = np.rint(blk.astype(np.float64) @ owners.T.astype(np.float64)).astype(np.int64)   # exact: |entries| <= 4, 24 terms
    same = (key(blk)[:, None] == kO[None, :])
    G[same] = -10 ** 6
    mx_eq_owner = max(mx_eq_owner, int(G.max()))
    mx_eq_head = max(mx_eq_head, int(np.abs(np.rint(blk.astype(np.float64) @ F.T.astype(np.float64))).max()))
assert mx_eq_owner <= 16, 'equator-cap: <u,u_owner> <= 16 for u != owner'
assert mx_eq_head <= 24, 'equator-layer: |<u,f>| <= 24 (both signs of every head are used)'
assert ok(sp.Rational(16, 1) / sp.sqrt(96) - 2) and ok(sp.Rational(24, 1) / sp.sqrt(192) - 2)
log('equator: max <u,u_owner> = %d (cap inner product 16/sqrt96 = %.4f), max |<u,f>| = %d (layer inner product 24/sqrt192 = %.4f)'
    % (mx_eq_owner, 16 / np.sqrt(96), mx_eq_head, 24 / np.sqrt(192)))

# ---------------------------------------------------------------- 6. cap vs layer  (the constraint that decides the layer)
# (u/sqrt12, (2/sqrt3) z).(f/sqrt24, s sqrt2 w) = <u,f>/sqrt288 + s (2 sqrt2/sqrt3) <z,w> <= 2
cosZW = [[sp.nsimplify(sp.simplify(z.dot(w))) for z in Z] for w in WS]
mx_cf = np.where(VAL[:, :, None], C @ F.T, 0)     # 24 x 496 x N, empty slots ignored
mcf = mx_cf.max(1)                                # 24 x N : max <u,f> over class i
amcf = np.maximum(mcf, -mx_cf.min(1))            # both signs of a head: max |<u,f>| over the class
CAPLAYER = {}                                    # (m, cos) -> the exact inner product, decided once
dirty = sorted({i for hh in range(N) for a in HD[hh] for i in range(24) if any(sp.N(cosZW[a][3 * i + p], 30) > 0 for p in range(3))})
worst = None
for i in range(24):
    for hh in range(N):
        for a in HD[hh]:
            for p in range(3):
                ky = (int(amcf[i, hh]), cosZW[a][3 * i + p])
                if ky not in CAPLAYER:
                    v = sp.simplify(sp.Integer(ky[0]) / sp.sqrt(288) + 2 * sp.sqrt(2) / sp.sqrt(3) * ky[1])
                    assert ok(v - 2), 'cap-layer (max <u,f> = %s, cos = %s)' % ky
                    CAPLAYER[ky] = v
                v = CAPLAYER[ky]
                if worst is None or sp.N(v, 30) > sp.N(worst, 30):
                    worst = v
log('cap-layer: cos(z,w_j) in %s; %d of the 24 triangles are dirty for some head; over every (dirty class, head) pair max |<u,f>| = %d (a conflict would be 24); worst inner product %.6f (bound 2)'
    % (sorted({str(c) for cw in cosZW for c in cw}), len(dirty),
       max([int(amcf[i, hh]) for hh in range(N) for a in HD[hh] for i in range(24)
            if any(sp.N(cosZW[a][3 * i + p], 30) > 0 for p in range(3))]), float(worst)))

# ---------------------------------------------------------------- 7. axis vs layer, layer vs layer
USED = sorted({a for hd in HD for a in hd})
for a in AX:
    for j in USED:
        assert ok(2 * sp.sqrt(2) * a.dot(WS[j]) - 2), 'axis-layer'
worst = max((sp.simplify(2 * sp.sqrt(2) * a.dot(WS[j])) for a in AX for j in USED), key=lambda e: sp.N(e, 30))
log('axis-layer: maximum inner product %s = %.6f (bound 2)' % (worst, float(worst)))
# layer-layer: (f/sqrt24, sqrt2 w_a).(f'/sqrt24, sqrt2 w_b) = <f,f'>/24 + 2 cos(w_a, w_b) <= 2
LAYERLAYER = {}
for i in range(N):
    for j in range(N):
        for a in HD[i]:
            for b in HD[j]:
                if i == j and a == b:
                    m = -48                        # the two signs of one head line at one direction
                elif i == j:
                    m = 48                         # one head vector at two of its directions
                else:
                    m = int(AGF[i, j])             # either sign of two distinct head lines
                ky = (m, cosWW[a][b])
                if ky not in LAYERLAYER:
                    assert ok(sp.Integer(ky[0]) / 24 + 2 * ky[1] - 2), 'layer-layer (<f,f_> = %s, cos = %s)' % ky
                    LAYERLAYER[ky] = True
log("layer-layer: all %d layer points pairwise admissible (largest |<f,f'>| between head lines %d)" % (LAYER, AGF.max() if N > 1 else 0))

# ---------------------------------------------------------------- 8. count
count = 196560 - NOWN + 3 * NOWN + 72 + LAYER
log('count = 196560 - %d + 3 * %d + 72 + %d = %d   (%.0f s)' % (NOWN, NOWN, LAYER, count, time.time() - t0))
print('ALL CHECKS PASS   K(30) >= %d' % count)
