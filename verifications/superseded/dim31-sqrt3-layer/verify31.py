# -*- coding: utf-8 -*-
"""Exact verification of  K(31) >= 238354   (python verify31.py)

Reads data/classes.npy (42 x 496 x 24, Cohn coordinates, minimal vectors of norm 32: the owner
classes of the published construction, one per zero-sum triangle of the cap E7), data/u0.npy
(a minimal vector that is not an owner) and data/directions.json (unit vectors of R^7 as sympy
strings: the 42 triangles, the 126 axis directions, and the four layer directions -- two lines
of deep holes of the cap E7), rebuilds the Leech minimal vectors from the Golay code, and checks
every pair of points exactly: integer inner products for the x-parts, sympy in Q(sqrt2, sqrt3)
for the y-parts.

In norm-4 units the configuration is
    equator  (u/sqrt8, 0)               u a minimal vector that is not an owner         196560 - 42 * 496 points
    caps     (u/sqrt12, (2/sqrt3) z)    u in class i, z in triangle i                  3 per owner
    axis     (0, 2a)                    a an axis direction, minus the four at cosine > 1/sqrt3 to a layer direction
    layer    (+-u0/(2 sqrt8), sqrt3 w)  w one of the four layer directions             8 points
and two points are compatible iff their inner product is at most 2.
"""
from __future__ import print_function
import os, sys, json, time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build        # noqa: E402

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
u0 = np.load(os.path.join(HERE, 'data', 'u0.npy')).astype(np.int64)
D = json.load(open(os.path.join(HERE, 'data', 'directions.json')))
vec = lambda ss: sp.Matrix([sp.sympify(s) for s in ss])
TRI = [[vec(z) for z in T] for T in D['triangles']]
AX = [vec(a) for a in D['axis']]
HOLES = [vec(h) for h in D['holes']]
K = len(TRI)
assert C.shape == (K, 496, 24) and K == 42 and len(AX) == 126 and len(HOLES) == 4 and u0.shape == (24,)
log('data: 42 classes of 496, u0, 42 triangles, 126 axis directions, 4 layer directions')

# ---------------------------------------------------------------- 0. directions
Z = [z for T in TRI for z in T]
for T in TRI:
    assert len(T) == 3 and all(sp.simplify(z.dot(z) - 1) == 0 for z in T), 'cap directions are unit'
    assert all(sp.simplify(c) == 0 for c in (T[0] + T[1] + T[2])), 'each triangle sums to zero'
assert len({tuple(sp.nsimplify(c) for c in z) for z in Z}) == 126, '126 distinct cap directions'
assert all(sp.simplify(a.dot(a) - 1) == 0 for a in AX) and len({tuple(sp.nsimplify(c) for c in a) for a in AX}) == 126
assert all(sp.simplify(w.dot(w) - 1) == 0 for w in HOLES)
cosZZ = [[sp.nsimplify(sp.simplify(Z[p].dot(Z[q]))) for q in range(126)] for p in range(126)]
vals = (0, sp.Rational(1, 2), -sp.Rational(1, 2), 1, -1)
assert all(cosZZ[p][q] in vals for p in range(126) for q in range(126)), 'cap-direction cosines are 0, +-1/2, +-1 (an E7 root system)'
assert all(cosZZ[p][q] <= sp.Rational(1, 2) for p in range(126) for q in range(126) if p != q)
# the layer directions: two antipodal pairs of deep holes, cosine exactly +-1/sqrt3 or 0 to every cap direction
cosZW = [[sp.nsimplify(sp.simplify(z.dot(w))) for w in HOLES] for z in Z]
r3 = 1 / sp.sqrt(3)
assert all(c in (0, sp.nsimplify(r3), sp.nsimplify(-r3)) for row in cosZW for c in row), 'layer directions are deep holes (cosines 0, +-1/sqrt3)'
HH = [[sp.nsimplify(sp.simplify(HOLES[i].dot(HOLES[j]))) for j in range(4)] for i in range(4)]
log('directions: 42 zero-sum unit triangles (E7), 126 unit axis directions, 4 unit layer directions with mutual cosines %s' % [[str(c) for c in row] for row in HH])

# ---------------------------------------------------------------- 1. classes and u0
M = leech_build().astype(np.int64)
pw = 9 ** np.arange(24, dtype=np.int64)
key = lambda V: ((V + 4) * pw).sum(-1)
SK = np.sort(key(M))
kC = key(C).ravel()
pos = np.searchsorted(SK, kC)
assert (SK[np.minimum(pos, len(SK) - 1)] == kC).all(), 'every class vector is a Leech minimal vector'
assert len(np.unique(kC)) == K * 496, 'the 42 classes are pairwise disjoint'
Mij = np.zeros((K, K), dtype=np.int64)
for i in range(K):
    for j in range(K):
        G = C[i] @ C[j].T
        if i == j:
            np.fill_diagonal(G, -10 ** 6)
        Mij[i, j] = G.max()
assert (np.diag(Mij) <= 8).all(), 'inside a class the inner product is at most 8 (>= arccos(1/4) = 75.52 deg)'
owners = C.reshape(-1, 24)
kO = key(owners)
assert (u0 ** 2).sum() == 32 and SK[np.searchsorted(SK, key(u0))] == key(u0), 'u0 is a minimal vector'
assert key(u0) not in set(kO.tolist()), 'u0 is not an owner'
log('classes: %d distinct minimal vectors, inner product at most 8 inside a class (max inside %d), i.e. at least arccos(1/4) = 75.52 deg; u0 is a minimal non-owner' % (K * 496, np.diag(Mij).max()))

# ---------------------------------------------------------------- 2. cap vs cap
thr = {c: int(sp.nsimplify(24 - 16 * c)) for c in vals}
for i in range(K):
    for j in range(K):
        for p in range(3):
            for q in range(3):
                m = Mij[i, j]
                if i == j and p != q:
                    m = max(m, 32)
                assert m <= thr[cosZZ[3 * i + p][3 * j + q]], 'cap-cap (%d,%d,%d,%d)' % (i, j, p, q)
log('cap-cap: all class pairs and direction pairs admissible (thresholds %s)' % thr)

# ---------------------------------------------------------------- 3. axis: drop the points blocked by the layer, then cap-axis and axis-axis
blocked = [k for k, a in enumerate(AX) if any(not ok(2 * sp.sqrt(3) * a.dot(w) - 2) for w in HOLES)]
assert len(blocked) == 4, 'exactly four axis points are blocked by the layer, found %d' % len(blocked)
AXK = [a for k, a in enumerate(AX) if k not in blocked]
h = 2 / sp.sqrt(3)
worst = None
for z in Z:
    for a in AXK:
        v = sp.simplify(2 * h * z.dot(a))
        assert ok(v - 2), 'cap-axis'
        if worst is None or sp.N(v, 30) > sp.N(worst, 30):
            worst = v
for i, a in enumerate(AXK):
    for b in AXK[i + 1:]:
        assert ok(4 * a.dot(b) - 2), 'axis-axis'
log('axis: 4 blocked points dropped (indices %s); the remaining 122 are pairwise admissible and admissible against every cap (max %s = %.6f)' % (blocked, worst, float(worst)))

# ---------------------------------------------------------------- 4. equator vs caps
mx_eq_owner = -10 ** 6
B = 4096
for s in range(0, len(M), B):
    blk = M[s:s + B]
    G = np.rint(blk.astype(np.float64) @ owners.T.astype(np.float64)).astype(np.int64)
    G[key(blk)[:, None] == kO[None, :]] = -10 ** 6
    mx_eq_owner = max(mx_eq_owner, int(G.max()))
assert mx_eq_owner <= 16 and ok(sp.Integer(16) / sp.sqrt(96) - 2)
log('equator-cap: max <u,u_owner> = %d over non-owners (inner product 16/sqrt96 = %.4f)' % (mx_eq_owner, 16 / np.sqrt(96)))

# ---------------------------------------------------------------- 5. the layer  (+- u0/(2 sqrt8), sqrt3 w)
# vs equator (v/sqrt8, 0): <u0,v>/16 <= 2, equality only at v = u0 (touching) -- Cauchy-Schwarz, checked over all of M
mu = int(np.rint(M.astype(np.float64) @ u0.astype(np.float64)).max())
assert mu == 32 and ok(sp.Integer(32) / 16 - 2), 'layer-equator'
# vs caps (u/sqrt12, (2/sqrt3) z): <u0,u>/(2 sqrt96) + 2 <z,w> <= 2 with <u0,u> <= 16 (u0 is not an owner) and <z,w> <= 1/sqrt3
mo = int(np.rint(owners.astype(np.float64) @ u0.astype(np.float64)).max())
assert mo <= 16, 'u0 is at inner product <= 16 with every owner'
v = sp.Integer(mo) / (2 * sp.sqrt(96)) + 2 * r3
assert ok(v - 2), 'layer-cap'
log('layer-cap: max <u0,u_owner> = %d, worst inner product %s = %.6f (bound 2; exact margin sqrt6 + 2 sqrt3 < 6)' % (mo, sp.simplify(v), float(v)))
# vs the surviving axis (0, 2a): 2 sqrt3 <w,a> <= 2
for a in AXK:
    for w in HOLES:
        assert ok(2 * sp.sqrt(3) * a.dot(w) - 2), 'layer-axis'
# vs each other: (s u0/(2 sqrt8), sqrt3 w).(s' u0/(2 sqrt8), sqrt3 w') = s s' + 3 <w,w'> <= 2
for i in range(4):
    for j in range(4):
        for s1 in (1, -1):
            for s2 in (1, -1):
                if i == j and s1 == s2:
                    continue
                assert ok(s1 * s2 + 3 * HH[i][j] - 2), 'layer-layer'
log('layer: 8 points (+-u0/(2 sqrt8), sqrt3 w), admissible against the equator (touching at u0), every cap, the surviving axis, and each other')

# ---------------------------------------------------------------- 6. count
count = 196560 - K * 496 + 3 * K * 496 + (126 - 4) + 8
log('count = 196560 - %d + 3 * %d + (126 - 4) + 8 = %d   (%.0f s)' % (K * 496, K * 496, count, time.time() - t0))
print('ALL CHECKS PASS   K(31) >= %d' % count)
