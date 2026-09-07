# -*- coding: utf-8 -*-
"""Exact verification of  K(26) >= 199632   (python verify26.py)

Reads data/heads26_Y.npy (integer heads, Cohn units, norm 192) and data/heads26_side.npy
(0 = directions 30/150/270 degrees, 1 = directions 90/210/330), rebuilds the Leech minimal
vectors from the Golay code, and checks every pair of points in exact arithmetic: integer
for head-equator and head-head, symbolic in Q(sqrt3) for the hexagon.  See lib/layered.py.
"""
from __future__ import print_function
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build          # noqa: E402
from layered import verify, hexagon              # noqa: E402

M = leech_build().astype(np.int64)
Y = np.load(os.path.join(HERE, 'data', 'heads26_Y.npy')).astype(np.int64)
side = np.load(os.path.join(HERE, 'data', 'heads26_side.npy'))
tri, axis, tau = hexagon()
count = verify(M, Y, side, tri, axis, tau, 26, lambda s: print('  ' + s))
print('ALL CHECKS PASS   K(26) >= %d' % count)
