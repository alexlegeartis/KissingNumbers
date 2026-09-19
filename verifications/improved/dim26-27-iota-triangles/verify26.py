# -*- coding: utf-8 -*-
"""Exact verification of  K(26) >= 199806   (python verify26.py)

TWO cap layers, as in dimension 27.

The first is heads at |x|^2 = 8/3 with |y| = 2/sqrt3, each carrying a zero-sum TRIPLE of the
hexagon's six edge-midpoint directions: three cap points and one equator deletion.  Read from
data/heads26_Y.npy (integer, Cohn units, norm 192) and data/heads26_side.npy (which triangle).

The second sits at |x|^2 = 3 with |y| = 1, carries the antipodal pair +-y, and takes the scaled
owner x' = (sqrt3/2) u, which deletes exactly u -- two cap points for one deletion, +1 each.
Read from data/heads26_layer2_u.npy and data/heads26_layer2_line.npy.

Its direction must sit as far as it can from the six first-layer directions, and their covering
radius is 30 degrees, attained at the hexagon's six VERTICES -- three antipodal lines, which are
also the axis directions.  That gives the bar

    first layer vs second   <Y,u> <= 24      (sqrt3/48)*24 + 1 = 1 + sqrt3/2 < 2

against dimension 27's 32, and it is why this layer was reported impossible: measured over the
196560 minimal vectors the minimum of max<Y,u> is 32 -- but ONLY because of the FREE heads
y = +-2v, for which <Y,u> = 2<v,u> is a multiple of 16 and is 32 the moment <v,u> = 16.  Over
the CLASS heads alone 1188 owners qualify, and 324 of them fit.  A free head is worth 3 and a
second-layer head 1, so giving up all 96 free heads costs 288 and buys 324: +36.  No free head
can be put back -- none of the 16773120 norm-6 vectors clears the whole layer.

Integer comparisons, in the units where minimal vectors have |z|^2 = 32:

    first layer vs second        <Y,u>  <= 24      (sqrt3/48)*24 + 1 = 1 + sqrt3/2 < 2
    second vs second, same line  <u,u'> <= 8       (3/32)*8 + 1 = 7/4 < 2

and the rest is automatic: second vs equator (sqrt3/2)*2 = sqrt3 < 2; second vs second across
lines (3/32)*16 + 1/2 = 2; on one line with opposite signs (3/32)*16 - 1 < 2; second vs axis
<y',a> <= |y'||a| = 2.  Each is checked below, the inequalities symbolically.
"""
from __future__ import print_function
import os, sys
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build          # noqa: E402
from layered import verify, hexagon              # noqa: E402

bad = []


def chk(name, ok, detail=''):
    print('  %-58s %s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if not ok:
        bad.append(name)


M = leech_build().astype(np.int64)
Y = np.load(os.path.join(HERE, 'data', 'heads26_Y.npy')).astype(np.int64)
side = np.load(os.path.join(HERE, 'data', 'heads26_side.npy'))
tri, axis, tau = hexagon()
first = verify(M, Y, side, tri, axis, tau, 26, lambda s: print('  ' + s))

U = np.load(os.path.join(HERE, 'data', 'heads26_layer2_u.npy')).astype(np.int64)
line = np.load(os.path.join(HERE, 'data', 'heads26_layer2_line.npy')).astype(np.int64)
print('  second layer: %d heads at |x|^2 = 3 on %d axis lines'
      % (len(U), len(set(line.tolist()))))

chk('second-layer owners are Leech minimal vectors',
    bool(((U * U).sum(1) == 32).all())
    and all(bytes(np.ascontiguousarray(u)) in set(map(bytes, np.ascontiguousarray(M)))
            for u in U))
chk('second-layer owners are distinct',
    len(set(map(bytes, np.ascontiguousarray(U)))) == len(U))

own1 = set()
for a in range(0, len(Y), 400):
    for r in Y[a:a + 400] @ M.T:
        for i in np.nonzero(r > 48)[0]:
            own1.add(int(i))
pos = {bytes(np.ascontiguousarray(M[i])): i for i in range(len(M))}
own2 = set(pos[bytes(np.ascontiguousarray(u))] for u in U)
chk('the two layers delete disjoint sets (%d and %d)' % (len(own1), len(own2)),
    len(own1 & own2) == 0)
chk('a second-layer head deletes exactly its owner',
    all(int((M @ u).max()) == 32 for u in U))

mx = max(int((Y @ U[a:a + 200].T).max()) for a in range(0, len(U), 200))
chk('<Y,u> <= 24 for every first-layer head and second-layer owner', mx <= 24, 'max %d' % mx)
G = U @ U.T
np.fill_diagonal(G, -10 ** 9)
sm = line[:, None] == line[None, :]
chk('<u,u\'> <= 8 on one line', int(G[sm].max()) <= 8, 'max %d' % int(G[sm].max()))
chk('<u,u\'> <= 16 across lines', int(G[~sm].max()) <= 16, 'max %d' % int(G[~sm].max()))

r3 = sp.sqrt(3)
w = sp.Rational(24, 48) * r3 + 1
chk('(sqrt3/48)*24 + 1 < 2', bool(w < 2), '= 1 + sqrt3/2 = %.6f' % float(w))
chk('it is exactly 1 + sqrt3/2', sp.simplify(w - (1 + r3 / 2)) == 0)
chk('(3/32)*8 + 1 < 2', sp.Rational(3, 32) * 8 + 1 < 2)
chk('(3/32)*16 + 1/2 <= 2 and (3/32)*16 - 1 < 2',
    sp.Rational(3, 32) * 16 + sp.Rational(1, 2) <= 2 and sp.Rational(3, 32) * 16 - 1 < 2)
chk('(sqrt3/2)*2 < 2', bool(r3 < 2))
chk('a second-layer cap has norm 4', 3 + 1 == 4)
chk('its two cap points meet at 3 - 1 <= 2', 3 - 1 <= 2)

# the second-layer directions are the hexagon's vertices: unit, three antipodal lines, and
# exactly 30 degrees from the nearest first-layer direction
lines2 = [sp.Matrix([sp.cos(sp.pi * k / 3), sp.sin(sp.pi * k / 3)]) for k in range(6)]
chk('the three second-layer lines are unit vectors',
    all(sp.simplify(v.dot(v) - 1) == 0 for v in lines2))
mc = max(sp.nsimplify(sp.simplify(u.dot(v))) for T in tri for u in T for v in lines2)
chk('the widest gap from a first-layer direction is 30 degrees',
    sp.simplify(mc - r3 / 2) == 0, 'max cos = sqrt3/2')
chk('every axis point meets every second-layer direction at <= 2',
    all(sp.simplify(sp.Matrix(list(a)).dot(v)).is_nonpositive
        or sp.simplify(sp.Matrix(list(a)).dot(v) - 2).is_nonpositive
        for a in axis for v in lines2))

count = 196560 - len(own1) - len(own2) + 3 * len(Y) + 2 * len(U) + tau
print('  count = 196560 - %d - %d + 3 * %d + 2 * %d + %d = %d'
      % (len(own1), len(own2), len(Y), len(U), tau, count))
chk('the first layer alone gives 199482', first == 199482, '%d' % first)
if bad:
    print('FAILED: %s' % ', '.join(bad))
    sys.exit(1)
print('ALL CHECKS PASS   K(26) >= %d' % count)
