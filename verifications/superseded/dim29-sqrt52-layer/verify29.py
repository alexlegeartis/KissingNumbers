# -*- coding: utf-8 -*-
"""Exact verification of  K(29) >= 209594   (python verify29.py)

Reads data/classes.npy (14 x 496 x 24, Cohn coordinates, minimal vectors of norm 32: the owner
classes, one per direction group of the cap D5), data/heads.npy (4 x 24, Leech vectors of norm 48:
the heads of the new layer) and data/directions.json (unit vectors of R^5 as sympy strings, in the
orthonormal D5 frame of Cohn's dimension-29 block: the 14 direction groups -- 12 zero-sum triangles
and 2 antipodal pairs -- the 10 axis directions, and the 32 deep holes of D5 with their sign
parity), rebuilds the Leech minimal vectors from the Golay code, and checks every pair of points
exactly: integer inner products for the x-parts, sympy in Q(sqrt2, sqrt3, sqrt5) for the y-parts.

In norm-4 units the configuration is

    equator  (u/sqrt8, 0)                     u a minimal vector that is not an owner   189616 points
    caps     (u/sqrt12, (2/sqrt3) z)          u in class g, z in group g                 19840 points
    axis     (0, 2a)                          a one of the 10 axis directions +-E_i          10 points
    layer    (s_j v_i/sqrt32, sqrt(5/2) w_j)  v_i a head, w_j a deep hole,                  128 points
                                              s_j = +1 on the 16 even holes, -1 on the odd

and two points are compatible iff their inner product is at most 2.

The layer is the whole content.  Its height h = sqrt(5/2) is chosen so that three thresholds
coincide for D5:

  * |x|^2 = 4 - h^2 = 3/2 is exactly the norm-6 scale, so the head is v/2 with v a norm-6 Leech
    vector and <v/2, u> <= 3/2 < 2 against every equator point: the layer deletes NOTHING
    (section 5);
  * the axis threshold 1/h = 2/sqrt10 is exactly D5's covering cosine on S^4, attained at the 32
    deep holes, so a head may use deep holes as directions at all (section 7);
  * the directions-per-head threshold 1 - 2/h^2 = 1/5 is exactly the cosine of two sign patterns
    at Hamming distance 2, so one head carries 16 directions instead of 4 (section 7).

The published dimension-29 configuration (Ma et al. 2025, arXiv:2511.13391, the PackingStar
configurations; the value in H. Cohn's table, whose coordinate file dimensions25-31.txt supplies
the block) has 209496 points: the same equator and caps, and a 40-point axis.  Here the axis is
replaced by the 10 points +-E_i, of which the layer blocks none, and the layer adds 128, so the
gain is 128 - 30 = 98.
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
V = np.load(os.path.join(HERE, 'data', 'heads.npy')).astype(np.int64)
D = json.load(open(os.path.join(HERE, 'data', 'directions.json')))
vec = lambda ss: sp.Matrix([sp.sympify(s) for s in ss])
GRP = [[vec(z) for z in g] for g in D['groups']]
AX = [vec(a) for a in D['axis']]
WS = [vec(w) for w in D['holes']]
PAR = list(D['parity'])
K = len(GRP)
NH = len(V)
assert C.shape == (14, 496, 24) and V.shape == (4, 24) and K == 14
assert len(AX) == 10 and len(WS) == 32 and len(PAR) == 32
NOWN = 14 * 496
NCAP = sum(len(g) for g in GRP) * 496
LAYER = NH * len(WS)
log('data: 14 classes of 496 (%d owners), %d heads of norm 48, %d direction groups of sizes %s, '
    '%d axis directions, %d deep holes, %d layer points'
    % (NOWN, NH, K, [len(g) for g in GRP], len(AX), len(WS), LAYER))

# ---------------------------------------------------------------- 0. the direction system
Z = [z for g in GRP for z in g]
assert len(Z) == 40, '40 cap directions'
for g in GRP:
    assert all(sp.simplify(z.dot(z) - 1) == 0 for z in g), 'cap directions are unit'
    if len(g) == 3:
        assert all(sp.simplify(c) == 0 for c in (g[0] + g[1] + g[2])), 'a triangle sums to zero'
    else:
        assert len(g) == 2 and all(sp.simplify(c) == 0 for c in (g[0] + g[1])), 'a pair is antipodal'
key5 = lambda v: tuple(sp.nsimplify(c) for c in v)
assert len({key5(z) for z in Z}) == 40, 'the 14 groups cover 40 distinct cap directions once each'
cosZZ = [[sp.nsimplify(sp.simplify(Z[p].dot(Z[q]))) for q in range(40)] for p in range(40)]
VALS = (0, sp.Rational(1, 2), -sp.Rational(1, 2), 1, -1)
assert all(cosZZ[p][q] in VALS for p in range(40) for q in range(40)), \
    'cap-direction cosines are 0, +-1/2, +-1 (a D5 root system)'
assert all(cosZZ[p][q] <= sp.Rational(1, 2) for p in range(40) for q in range(40) if p != q), \
    'cap directions pairwise at >= 60 degrees'
assert all(sp.simplify(a.dot(a) - 1) == 0 for a in AX) and len({key5(a) for a in AX}) == 10
assert all(sp.simplify(w.dot(w) - 1) == 0 for w in WS) and len({key5(w) for w in WS}) == 32
# the deep holes: at D5's covering cosine 2/sqrt10 from every cap direction, and nothing closer
cov = sp.nsimplify(2 / sp.sqrt(10))
cosZW = [[sp.nsimplify(sp.simplify(z.dot(w))) for w in WS] for z in Z]
seenZW = {c for row in cosZW for c in row}
assert seenZW == {0, cov, -cov}, \
    'every deep hole is at cosine 0 or +-2/sqrt10 from every cap direction, found %s' % seenZW
cosWW = [[sp.nsimplify(sp.simplify(WS[i].dot(WS[j]))) for j in range(32)] for i in range(32)]
assert {c for row in cosWW for c in row} == {sp.Rational(k, 5) for k in (-5, -3, -1, 1, 3, 5)}, \
    'deep-hole cosines are the sign-pattern distances (5-2d)/5'
assert all(PAR[i] == PAR[j] for i in range(32) for j in range(32)
           if cosWW[i][j] == sp.Rational(1, 5)), \
    'cosine 1/5 (Hamming distance 2) happens only inside one parity class'
log('directions: %d groups (%d zero-sum triangles, %d antipodal pairs) covering 40 unit cap '
    'directions with cosines in {0,+-1/2,+-1} (a D5 root system); %d unit axis directions; '
    '32 deep holes at cosine exactly 2/sqrt10 = %.9f from every cap direction (D5 covering radius)'
    % (K, sum(1 for g in GRP if len(g) == 3), sum(1 for g in GRP if len(g) == 2), len(AX), float(cov)))

# ---------------------------------------------------------------- 1. the classes
M = leech_build().astype(np.int64)
pw = 9 ** np.arange(24, dtype=np.int64)
key = lambda X: ((X + 4) * pw).sum(-1)
SK = np.sort(key(M))
kC = key(C).ravel()
pos = np.searchsorted(SK, kC)
assert (SK[np.minimum(pos, len(SK) - 1)] == kC).all(), 'every class vector is a Leech minimal vector'
assert len(np.unique(kC)) == NOWN, 'the 14 classes are pairwise disjoint'
assert ((C ** 2).sum(2) == 32).all(), 'every owner has norm 32'
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
log('classes: %d distinct minimal vectors, pairwise disjoint; inside a class the largest inner '
    'product between distinct vectors is %d (bound 8: the owners of a class share every direction '
    'of their group, so 60-degree-freeness is not enough)' % (NOWN, np.diag(Mij).max()))

# ---------------------------------------------------------------- 2. the heads
assert ((V ** 2).sum(1) == 48).all(), 'heads have norm 48 (norm 6 in the standard scaling)'
assert all(in_lattice(v) for v in V), 'heads are Leech vectors'
GV = V @ V.T
assert (GV - np.diag(np.diag(GV)) == -16 * (1 - np.eye(NH, dtype=np.int64))).all(), \
    'the heads are pairwise at inner product -16'
assert np.abs(V.sum(0)).max() == 0, 'the four heads sum to zero'
# the simplex bound at |x|^2 = 3/2: k heads may share one direction only if <v,v'>/32 + 5/2 <= 2,
# that is <v,v'> <= -16, i.e. cosine <= -1/3, and k unit vectors at pairwise cosine <= -1/3 force
# k <= 4 with equality only for a zero sum.  Both are met exactly.
assert NH == 4 and sp.simplify(sp.Integer(-16) / 32 + sp.Rational(5, 2) - 2) == 0, \
    'four heads per direction is the simplex bound at |x|^2 = 3/2, and it is tight'
log('heads: %d Leech vectors of norm 48, Gram 48 on the diagonal and -16 off it, summing to zero '
    '-- the simplex bound k <= 2/(2 - |x|^2) = 4 at |x|^2 = 3/2, attained and tight' % NH)

# ---------------------------------------------------------------- 3. cap vs cap
# (u/sqrt12, (2/sqrt3) z).(u'/sqrt12, (2/sqrt3) z') = <u,u'>/12 + (4/3) cos <= 2  <=>  <u,u'> <= 24 - 16 cos
thr = {c: int(sp.nsimplify(24 - 16 * c)) for c in VALS}
tight = 0
off = np.cumsum([0] + [len(g) for g in GRP])
for i in range(K):
    for j in range(K):
        for p in range(len(GRP[i])):
            for q in range(len(GRP[j])):
                c = cosZZ[off[i] + p][off[j] + q]
                m = Mij[i, j]
                if i == j and p != q:
                    m = max(m, 32)                     # the same owner on two directions of its group
                assert m <= thr[c], 'cap-cap (classes %d,%d directions %d,%d)' % (i, j, p, q)
                tight += int(m == thr[c])
log('cap-cap: all class pairs and direction pairs admissible (thresholds %s); %d tight cases'
    % (thr, tight))

# ---------------------------------------------------------------- 4. cap vs axis, axis vs axis
h = 2 / sp.sqrt(3)
worst = max((sp.simplify(2 * h * z.dot(a)) for z in Z for a in AX), key=lambda e: sp.N(e, 30))
assert ok(worst - 2), 'cap-axis'
for i, a in enumerate(AX):
    for b in AX[i + 1:]:
        assert ok(4 * a.dot(b) - 2), 'axis-axis'
log('cap-axis: maximum inner product %s = %.6f (bound 2).  axis-axis: %d points of norm 4, '
    'pairwise inner product <= 2' % (worst, float(worst), len(AX)))

# ---------------------------------------------------------------- 5. equator vs caps, equator vs layer
# the section that proves the layer is deletion-free
mx_eq_owner = -10 ** 6
mx_eq_head = -10 ** 6
B = 1024
for s in range(0, len(M), B):
    blk = M[s:s + B]
    G = np.rint(blk.astype(np.float64) @ owners.T.astype(np.float64)).astype(np.int64)  # exact: |entries| <= 4
    G[key(blk)[:, None] == kO[None, :]] = -10 ** 6
    mx_eq_owner = max(mx_eq_owner, int(G.max()))
    mx_eq_head = max(mx_eq_head, int(np.abs(np.rint(blk.astype(np.float64) @ V.T.astype(np.float64))).max()))
assert mx_eq_owner <= 16 and ok(sp.Integer(16) / sp.sqrt(96) - 2), 'equator-cap'
assert mx_eq_head <= 24 and ok(sp.Integer(24) / 16 - 2), 'equator-layer'
log('equator: max <u,u_owner> = %d over non-owners (cap inner product 16/sqrt96 = %.4f); '
    'max |<u,v_i>| = %d over the WHOLE shell (layer inner product 24/16 = %.4f) -- the layer '
    'deletes no equator point at all' % (mx_eq_owner, 16 / np.sqrt(96), mx_eq_head, 24 / 16.0))

# ---------------------------------------------------------------- 6. cap vs layer (the constraint that decides the layer)
# (u/sqrt12, (2/sqrt3) z).(s v/sqrt32, sqrt(5/2) w) = <s v,u>/sqrt384 + (2/sqrt3) sqrt(5/2) <z,w> <= 2
H = sp.sqrt(sp.Rational(5, 2))
mvu = int(np.abs(owners @ V.T).max())
assert mvu <= 16, 'every owner avoids the forbidden set |<v_i,u>| = 24'
CAPLAYER = {}
worst = None
for j in range(32):
    for gi in range(K):
        for p in range(len(GRP[gi])):
            c = cosZW[off[gi] + p][j]
            if c not in CAPLAYER:
                e = sp.simplify(sp.Integer(mvu) / sp.sqrt(384) + h * H * c)
                assert ok(e - 2), 'cap-layer (max |<v,u>| = %d, cos = %s)' % (mvu, c)
                CAPLAYER[c] = e
            e = CAPLAYER[c]
            if worst is None or sp.N(e, 30) > sp.N(worst, 30):
                worst = e
dirty = sum(1 for gi in range(K)
            if any(sp.N(cosZW[off[gi] + p][j], 30) > 0
                   for p in range(len(GRP[gi])) for j in range(32)))
assert dirty == K, 'every group is dirty, so every class had to avoid the forbidden set'
log('cap-layer: max |<v_i,u>| over all %d owners and all %d heads is %d -- 24 would be a conflict '
    '(the forbidden set U = {u : |<v_i,u>| = 24}, |U| = 4320, which every class avoids); all %d '
    'groups are dirty for some deep hole; worst inner product %s = %.9f (bound 2, headroom %.4f)'
    % (NOWN, NH, mvu, dirty, sp.simplify(worst), float(worst), 2 - float(worst)))

# ---------------------------------------------------------------- 7. axis vs layer, layer vs layer
# (0,2a).(x, sqrt(5/2) w) = 2 sqrt(5/2) <a,w> = sqrt10 <a,w> <= 2  <=>  <a,w> <= 2/sqrt10
worst = max((sp.simplify(2 * H * a.dot(w)) for a in AX for w in WS), key=lambda e: sp.N(e, 30))
assert ok(worst - 2), 'axis-layer'
assert all(ok(2 * H * a.dot(w) - 2) for a in AX for w in WS), 'axis-layer: not one of the 10 is blocked'
log('axis-layer: maximum inner product %s = %.6f (bound 2; the threshold is cos <= 2/sqrt10, '
    'exactly D5\'s covering cosine) -- none of the %d axis points is blocked'
    % (sp.simplify(worst), float(worst), len(AX)))
# layer-layer: (s v_i/sqrt32, sqrt(5/2) w).(s' v_j/sqrt32, sqrt(5/2) w') = s s' <v_i,v_j>/32 + (5/2) cos
SGN = [1 if PAR[j] == 0 else -1 for j in range(32)]
LL = {}
touch = 0
for i in range(NH):
    for j in range(NH):
        for a in range(32):
            for b in range(32):
                if i == j and a == b:
                    continue
                ky = (int(SGN[a] * SGN[b] * GV[i, j]), cosWW[a][b])
                if ky not in LL:
                    e = sp.simplify(sp.Integer(ky[0]) / 32 + sp.Rational(5, 2) * ky[1])
                    assert ok(e - 2), 'layer-layer (<v,v_> = %s, cos = %s)' % ky
                    LL[ky] = e
                touch += int(LL[ky] == 2)
log('layer-layer: all %d layer points pairwise admissible over %d distinct (inner product, cosine) '
    'cases; %d ordered pairs touch at exactly 2' % (LAYER, len(LL), touch))
log('   the touching families: %s'
    % '; '.join('<v,v\'> = %+d at cos %s' % (k[0], k[1]) for k, e in sorted(LL.items(), key=str) if e == 2))

# ---------------------------------------------------------------- 8. the four families are disjoint
# Everything above is about ANGLES; the count also needs the four families to be disjoint as point
# sets, which is a different kind of statement.  Each family has a constant (|x|^2, |y|^2), and the
# four are distinct, so no point of one is a point of another:
SIG = {'equator': (sp.Integer(4), sp.Integer(0)),
       'caps': (sp.Integer(32) / 12, (2 / sp.sqrt(3)) ** 2),
       'axis': (sp.Integer(0), sp.Integer(4)),
       'layer': (sp.Integer(48) / 32, H ** 2)}
for nm, (a, b) in SIG.items():
    assert sp.simplify(a + b - 4) == 0, '%s is not on the unit sphere of radius 2' % nm
assert len({(sp.nsimplify(a), sp.nsimplify(b)) for a, b in SIG.values()}) == 4,     'two families share a norm signature, so the count could double-count'
# and inside each family: the equator is the minimal vectors minus the owners (both counted by
# distinct keys above), the caps are (owner, direction of its own group) pairs with the owners
# pairwise distinct and the groups covering each direction once, the axis is 10 distinct unit
# vectors, and two layer points that coincided would appear as an inner product of 4 in section 7.
log('disjointness: the four families have distinct (|x|^2, |y|^2) = %s, so the count adds no point '
    'twice' % sorted((str(sp.nsimplify(a)), str(sp.nsimplify(b))) for a, b in SIG.values()))

# ---------------------------------------------------------------- 9. count
count = 196560 - NOWN + NCAP + len(AX) + LAYER
assert count == 209594
log('count = 196560 - %d + %d + %d + %d = %d   (%.0f s)'
    % (NOWN, NCAP, len(AX), LAYER, count, time.time() - t0))
print('ALL CHECKS PASS   K(29) >= %d   (published 209496 -- Ma et al. 2025, arXiv:2511.13391 --'
      ' gain +%d)' % (count, count - 209496))
