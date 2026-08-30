#!/usr/bin/env python3
"""Does the PRESENTATION of the Gram matrix change the bound?  Measured, at k = 3.

    python scripts/presentation_invariance.py        # ~1 min

WHY THIS EXISTS.  scan_k3.py reports a different value for almost every one of the sixteen
realizable determinant-256 Gram matrices -- 625 956 480 up to 627 822 180 -- and scan_k4.py
reports 98 distinct values across 98 determinant-1024 Grams, a spread of 0.7%.  That was read
as the bound being basis-dependent among isometric Grams.  It is not.  The sixteen Grams are
sixteen bases of ONE lattice, A_3 scaled to minimum 8, and the programme they define is the
same programme; the spread is entirely the quality of the dual that the numeric solve hands
to the certifier.

WHAT GOES WRONG.  HiGHS's presolve, on these systems, returns points that violate the
equality constraints by tens of vectors.  A dual read off such a solve is still repaired into
exact feasibility by the certifier -- so no bound was ever WRONG, every one of those numbers
is a valid lower bound -- but it is a weak witness, and how weak depends on the basis, which
is what made the spread look structural.

WHAT THIS SHOWS.  Running the same engine with extra_solves=True, which adds the
presolve=False attempts and keeps every dual instead of the first, all sixteen presentations
certify the identical value.  Since the certificate is the maximum over candidate duals and
each is verified exactly, adding candidates can only raise a bound; the default path is left
alone so that no shipped number moves.

The consequence for the claims is that "which Gram?" is settled by the DETERMINANT alone --
Hermite's floor -- and the choice of basis at that determinant costs nothing.
"""
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from kpoint_lp import run3   # noqa: E402  the one shared engine

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
M1G = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}

G72 = np.load(os.path.join(PKG, '..', 'dim73-95-gamma72-caps', 'data',
                           'gamma72_gram.npy')).astype(np.int64)
U = np.load(os.path.join(PKG, 'data', 'gamma72_d4_k4.npy')).astype(np.int64)
A4 = np.array([a for a in itertools.product(range(-2, 3), repeat=4) if any(a)], np.int64)
V = A4 @ U
V = V[np.einsum('ij,jk,ik->i', V, G72, V) == 8]
IP = V @ G72 @ V.T
assert len(V) == 24, len(V)

grams = {}
for t in itertools.combinations(range(24), 3):
    M = IP[np.ix_(t, t)]
    if not (np.diag(M) == 8).all():
        continue
    if int(round(np.linalg.det(M.astype(float)))) != 256:
        continue
    grams.setdefault(tuple(map(tuple, M.tolist())), t)

print(__doc__)
print("=" * 78)
print("   %d realizable determinant-256 Gram matrices, all bases of one scaled A_3" % len(grams))
print()
print("   Gram (off-diagonals)   default path   every dual kept   agree")
rows = []
t0 = time.time()
for g in sorted(grams):
    G = [list(r) for r in g]
    a = run3(72, 6218175600, 8, [1, 2, 3, 4], M1G, G)
    b = run3(72, 6218175600, 8, [1, 2, 3, 4], M1G, G, extra_solves=True)
    ca = a[2] if a and a[0] == 'cert' else None
    cb = b[2] if b and b[0] == 'cert' else None
    rows.append((ca, cb))
    print("   (%3d,%3d,%3d)            %-13s  %-16s  %s"
          % (g[0][1], g[0][2], g[1][2], ca, cb, "yes" if cb == rows[0][1] else "NO"))

best = {cb for _, cb in rows}
print()
print("   distinct certificates on the default path : %d" % len({ca for ca, _ in rows}))
print("   distinct certificates with every dual kept: %d" % len(best))
ok = len(best) == 1
print()
if ok:
    print("   ALL %d PRESENTATIONS CERTIFY %d.  The presentation costs nothing; the spread"
          % (len(rows), best.pop()))
    print("   scan_k3.py and scan_k4.py report is the numeric dual, not the Gram matrix.")
else:
    print("   *** the presentations do NOT agree: %s" % sorted(best))
print("   (%.0fs)" % (time.time() - t0))
sys.exit(0 if ok else 1)
