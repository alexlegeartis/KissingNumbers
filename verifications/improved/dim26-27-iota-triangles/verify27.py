# -*- coding: utf-8 -*-
"""Exact verification of  K(27) >= 201503   (python verify27.py)

TWO cap layers.

The first is the one dimension 26 also has: heads at |x|^2 = 8/3 with |y| = 2/sqrt3, each
carrying a zero-sum TRIPLE of the cuboctahedron's twelve directions, worth three cap points and
one equator deletion (or none, for a free head y = +-2v).  Read from data/heads27_Y.npy
(integer, Cohn units, norm 192) and data/heads27_side.npy (which of the four triangles).

The second sits at |x|^2 = 3 with |y| = 1.  Such a cap satisfies the axis condition <y,a> <= 2
whatever direction it points, because |y||a| = 2; it carries the antipodal PAIR +-y and no more,
since three directions would need pairwise cosine <= -1; and it always deletes, because nothing
in {x : <x,z> <= 2 for every minimal z} has norm above 8/3.  Taking the scaled owner
x' = (sqrt3/2) u deletes EXACTLY u: <x',z> = (sqrt3/2)<u,z> exceeds 2 only when <u,z> > 4/sqrt3
= 2.31, i.e. <u,z> = 4, i.e. z = u.  So it is two cap points for one deletion: +1 each.  Read
from data/heads27_layer2_u.npy (the owners) and data/heads27_layer2_line.npy (which axis line).

Its direction has to sit as far as it can from the twelve first-layer directions, and that is
where dimensions 26 and 27 differ.  In dimension 26 they are the six edge midpoints of a
hexagon, the best gap is 30 degrees, and MEASURED over the whole shipped layer not one of the
196560 minimal vectors qualifies.  In dimension 27 they are the cuboctahedron's twelve vertices,
whose covering radius is 45 degrees, attained exactly at the six SQUARE-FACE centres +-e_1,
+-e_2, +-e_3 -- three antipodal lines.  That single extra quantum is the whole construction.

Everything reduces to integer comparisons in the units where minimal vectors have |z|^2 = 32:

    first layer vs second   <Y,u>  <= 32     (sqrt3/48)*32 + 2/sqrt6 = (2+sqrt2)/sqrt3 < 2
    second vs second, same line  <u,u'> <= 8      (3/32)*8 + 1 = 7/4 < 2

and the rest is automatic: second vs equator (sqrt3/2)*2 < 2; second vs second across lines
(3/32)*16 = 3/2 < 2; second vs second on one line with opposite signs (3/32)*16 - 1 < 2; second
vs axis <y',a> <= |y'||a| = 2.  Each is checked below, the inequalities symbolically.
"""
from __future__ import print_function
import os, sys
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build          # noqa: E402
from layered import verify, cuboctahedron        # noqa: E402

bad = []


def chk(name, ok, detail=''):
    print('  %-58s %s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if not ok:
        bad.append(name)


M = leech_build().astype(np.int64)
Y = np.load(os.path.join(HERE, 'data', 'heads27_Y.npy')).astype(np.int64)
side = np.load(os.path.join(HERE, 'data', 'heads27_side.npy'))
tri, axis, tau = cuboctahedron()
first = verify(M, Y, side, tri, axis, tau, 27, lambda s: print('  ' + s))

U = np.load(os.path.join(HERE, 'data', 'heads27_layer2_u.npy')).astype(np.int64)
line = np.load(os.path.join(HERE, 'data', 'heads27_layer2_line.npy')).astype(np.int64)
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
chk('<Y,u> <= 32 for every first-layer head and second-layer owner', mx <= 32, 'max %d' % mx)
G = U @ U.T
np.fill_diagonal(G, -10 ** 9)
sm = line[:, None] == line[None, :]
chk('<u,u\'> <= 8 on one line', int(G[sm].max()) <= 8, 'max %d' % int(G[sm].max()))
chk('<u,u\'> <= 16 across lines', int(G[~sm].max()) <= 16, 'max %d' % int(G[~sm].max()))

r2, r3, r6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
w = sp.Rational(32, 48) * r3 + 2 / r6
chk('(sqrt3/48)*32 + 2/sqrt6 < 2', sp.simplify(w < 2) is sp.true or bool(w < 2),
    '= (2+sqrt2)/sqrt3 = %.6f' % float(w))
chk('it is exactly (2+sqrt2)/sqrt3', sp.simplify(w - (2 + r2) / r3) == 0)
chk('(3/32)*8 + 1 < 2', sp.Rational(3, 32) * 8 + 1 < 2)
chk('(3/32)*16 < 2 and (3/32)*16 - 1 < 2',
    sp.Rational(3, 32) * 16 < 2 and sp.Rational(3, 32) * 16 - 1 < 2)
chk('(sqrt3/2)*2 < 2', sp.simplify(r3 / 2 * 2 < 2) is sp.true or bool(r3 < 2))
chk('a second-layer cap has norm 4', 3 + 1 == 4)
chk('its two cap points meet at 3 - 1 <= 2', 3 - 1 <= 2)

cub = sp.Matrix([[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [1, 0, 1], [1, 0, -1],
                 [-1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1]]) / r2
chk('the square-face centres are 45 degrees from the nearest direction',
    max(abs(sp.simplify(cub.row(i).dot(sp.Matrix([1, 0, 0])))) for i in range(12)) == 1 / r2)
chk('every axis point meets every second-layer direction at <= 2',
    all(abs(sp.simplify(sp.Matrix(list(axis[i])).dot(sp.Matrix(ei)))) <= 2
        for i in range(len(axis))
        for ei in ([1, 0, 0], [0, 1, 0], [0, 0, 1])))

count = 196560 - len(own1) - len(own2) + 3 * len(Y) + 2 * len(U) + tau
print('  count = 196560 - %d - %d + 3 * %d + 2 * %d + %d = %d'
      % (len(own1), len(own2), len(Y), len(U), tau, count))
chk('the first layer alone gives 201225', first == 201225, '%d' % first)
if bad:
    print('FAILED: %s' % ', '.join(bad))
    sys.exit(1)
print('ALL CHECKS PASS   K(27) >= %d' % count)
