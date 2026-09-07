# -*- coding: utf-8 -*-
"""Exact verification of  K(27) >= 201010   (python verify27.py)

Reads data/heads27_Y.npy (integer heads, Cohn units, norm 192) and data/heads27_side.npy
(index of the direction triangle, 0..3, of the cuboctahedron split into four zero-sum
triangles), rebuilds the Leech minimal vectors from the Golay code, and checks every pair of
points exactly: integer for head-equator and head-head, symbolic in Q(sqrt2, sqrt3) for the
axis, which is the cuboctahedron rotated by 45 degrees about the third coordinate axis.
"""
from __future__ import print_function
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build          # noqa: E402
from layered import verify, cuboctahedron        # noqa: E402

M = leech_build().astype(np.int64)
Y = np.load(os.path.join(HERE, 'data', 'heads27_Y.npy')).astype(np.int64)
side = np.load(os.path.join(HERE, 'data', 'heads27_side.npy'))
tri, axis, tau = cuboctahedron()
count = verify(M, Y, side, tri, axis, tau, 27, lambda s: print('  ' + s))
print('ALL CHECKS PASS   K(27) >= %d' % count)
