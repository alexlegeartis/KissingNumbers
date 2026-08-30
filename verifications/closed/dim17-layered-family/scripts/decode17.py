"""Decode the dimension-17 record's layer structure from Cohn's coordinates.

The arithmetic 4320 + 2*704 + 2 = 5730 is suggestive, but suggestive is not decoded.  This
looks for an actual layering axis: a direction z such that the heights <z,x> over the 5730
points take few values, with an equator at height 0 and symmetric cap layers.  Every point of
the configuration is tried as an axis, plus the axis is polished by averaging the extreme layer.
"""
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from extract import blocks, SRC

B = [b for b in blocks(SRC) if b['dim'] == 17 and b['n'] == 5730]
assert B, "no 5730-point dimension-17 block"
b = B[0]
X = np.array(b['pts'], dtype=float) * np.sqrt(np.array(b['c'], dtype=float))[None, :]
U = X / np.linalg.norm(X, axis=1)[:, None]
print(__doc__)
print("dimension 17, %d points; max inner product %.9f"
      % (len(U), (U @ U.T - 2 * np.eye(len(U))).max()))

best = None
for i in range(len(U)):
    h = U @ U[i]
    q = np.round(h * 720).astype(np.int64)
    lay = collections.Counter(q.tolist())
    if best is None or len(lay) < best[0]:
        best = (len(lay), i, lay)
nl, i0, lay = best
print("\nfewest distinct heights over all %d candidate axes: %d, at point %d" % (len(U), nl, i0))
for v in sorted(lay):
    print("    height %+8.5f   count %6d" % (v / 720.0, lay[v]))

h = U @ U[i0]
eq = np.abs(h) < 1e-9
print("\nequator (height 0): %d points" % eq.sum())
if eq.sum():
    E = U[eq]
    print("   equator max inner product: %.9f" % (E @ E.T - 2 * np.eye(eq.sum())).max())
for lvl in sorted({round(float(x), 6) for x in h if abs(x) > 1e-9}):
    sel = np.abs(h - lvl) < 1e-7
    C = U[sel]
    if sel.sum() < 2:
        continue
    Cp = C - np.outer(C @ U[i0], U[i0])
    nz = np.linalg.norm(Cp, axis=1)
    if nz.min() < 1e-9:
        continue
    Cp = Cp / nz[:, None]
    print("   layer at %+8.5f : %5d points, projected max inner product %.9f"
          % (lvl, sel.sum(), (Cp @ Cp.T - 2 * np.eye(sel.sum())).max()))
