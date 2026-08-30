#!/usr/bin/env python3
"""Dimension 96 = 72 + 24: the Gamma_72 cap construction at k = 24, and what it is worth.

k = 24 is where the ladder of dimensions 73-95 SATURATES.  For k <= 23 the number of parts
the direction set W admits is the binding constraint; at k = 24, W = the 196560 minimal
vectors of the Leech lattice admits floor(196560/3) = 65520 zero-sum triples, far more than
the 31000-32000 pairwise-disjoint Gamma_72 classes that exist, so from k = 24 upward the
CLASS FAMILY alone decides the bound and every dimension 96, 97, ... inherits the same gain.

The surplus is not wasted: every direction of W not used by a triple is a legal POLE, since
a pole meets a cap at sqrt(1-t)*<z,w> = <z,w>/sqrt(3) <= 1/2 iff <z,w> <= sqrt(3)/2 = 0.866,
and distinct non-parallel Leech minimal vectors have |cos| <= 1/2.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DELIV = os.path.join(HERE, '..', '..', '..')
sys.path.insert(0, os.path.join(DELIV, 'common'))
from leech import load

N72 = 6218175600                       # tau(72), Nebe's Gamma_72
TAU24 = 196560                         # tau(24), the Leech lattice

# ---------------------------------------------------------------- 1. the direction set
A = load().astype(np.int64)            # Cohn's scaling: norm^2 = 32
T = np.load(os.path.join(HERE, 'data', 'triples24.npy')).astype(np.int64)
assert A.shape == (TAU24, 24)
assert ((A ** 2).sum(1) == 32).all(), "some direction is not a minimal vector"

# every triple: three distinct minimal vectors summing to zero, pairwise cos = -1/2
a, b, c = A[T[:, 0]], A[T[:, 1]], A[T[:, 2]]
assert (a + b + c == 0).all(), "a triple does not sum to zero"
for x, y in ((a, b), (a, c), (b, c)):
    ip = (x * y).sum(1)
    assert (ip == -16).all(), "a triple is not at 120 degrees"      # -16/32 = -1/2
flat = T.ravel()
assert len(np.unique(flat)) == len(flat), "the triples are not disjoint"
print("zero-sum triples verified exactly: %d, covering %d of %d directions"
      % (len(T), 3 * len(T), TAU24))

# ---------------------------------------------------------------- 2. the class families
FAMILIES = []
for nm, fn in (('images(Aut Gamma_72)', 'class_sizes_32000.npy'),
               ('greedy(spread pool)', 'class_sizes_gpu.npy')):
    fp = os.path.join(HERE, '..', '..', 'improved', 'dim73-95-gamma72-caps', 'data', fn)
    s = np.sort(np.load(fp).astype(np.int64))[::-1]
    FAMILIES.append((nm, s, np.concatenate([[0], np.cumsum(s)])))

# ---------------------------------------------------------------- 3. the bound
print("\n  family                     classes  triples used   poles     gain          tau(96) >=")
best = None
for nm, s, cum in FAMILIES:
    n_tri = min(len(s), len(T))                 # one triple per class; both are the binder
    assert 3 * n_tri <= TAU24
    gain = 4 * int(cum[n_tri])                  # 2*|C_i|*(|Z_i|-1) with |Z_i| = 3
    poles = TAU24 - 3 * n_tri                   # every unused direction is a legal pole
    tot = N72 + gain + poles
    print("  %-24s %7d  %11d  %6d  %11d   %d" % (nm, len(s), n_tri, poles, gain + poles, tot))
    if best is None or tot > best[0]:
        best = (tot, nm, n_tri, gain, poles)

tot, nm, n_tri, gain, poles = best
print("\n  best: %s, %d triples, %d poles" % (nm, n_tri, poles))
print("  tau(96) >= %d + %d + %d = %d" % (N72, gain, poles, tot))

# ---------------------------------------------------------------- 4. what it has to beat
floor_ds = N72 + TAU24                                   # Gamma_72 (+) Leech, min 8
ers_lvl0 = 1 << 33                                       # A(96,24) >= 2^33, [96,33,24]
print("\n  direct-sum floor  Gamma_72 (+) Leech      : %d" % floor_ds)
print("  Edel-Rains-Sloane level 0 alone, A(96,24) : %d   ([96,33,24], codetables.de)" % ers_lvl0)
print("\n  cap construction beats the direct sum by %d (x%.4f)"
      % (tot - floor_ds, tot / floor_ds))
print("  cap construction is BELOW Edel-Rains-Sloane by %d (x%.4f)"
      % (ers_lvl0 - tot, tot / ers_lvl0))

# ---------------------------------------------------------------- 5. what k = 24 hands back
# At k = 23 the scheme needs T + P = 30990 + 82 = 31072 parts.  This section used to price the
# 72 classes the greedy family was SHORT of that, because dimension 95 was falling back on the
# images family as a result.  The family was extended to 32000 classes on 2026-08-26 and the
# shortfall is gone, so the section now reports the realised comparison instead of a forecast.
# The old forecast is kept as a line, because it turned out to be 0.8% low and that is worth
# knowing about this kind of estimate: it assumed continued classes no smaller than the
# smallest already seen (1938), and the realised ones ran 1932-2118 with mean 2049.
T95, P95, POLES95 = 30990, 82, 248
def gain95(s):
    c = np.concatenate([[0], np.cumsum(s)])
    return 4 * int(c[T95]) + 2 * int(c[T95 + P95] - c[T95]) if len(s) >= T95 + P95 else None
img = FAMILIES[0][1]
fre = FAMILIES[1][1]
print("")
print("  dimension 95, for comparison:")
print("    images(Aut Gamma_72), %d classes, mean %.0f : tau(95) >= %d"
      % (len(img), img.mean(), N72 + gain95(img) + POLES95))
if len(fre) >= T95 + P95:
    g_img, g_fre = gain95(img), gain95(fre)
    print("    greedy(spread pool), %d classes, mean %.0f : tau(95) >= %d   (+%d, x%.4f)"
          % (len(fre), fre.mean(), N72 + g_fre + POLES95, g_fre - g_img, g_fre / g_img))
    print("    the gap is the CLASS SIZE, not the count: %.0f lines a class against %.0f."
          % (fre.mean(), img.mean()))
    print("    (forecast before the extension: +69289824, which was 0.8% low)")
else:
    short = T95 + P95 - len(fre)
    ext = np.concatenate([fre, np.full(short, fre.min(), dtype=np.int64)])
    print("    greedy(spread pool) has %d, needs %d -- %d short, so it cannot serve dim 95"
          % (len(fre), T95 + P95, short))
    print("    were it extended by %d classes of its own minimum size %d:"
          % (short, fre.min()))
    print("                                       tau(95) >= %d   (+%d)"
          % (N72 + gain95(ext) + POLES95, gain95(ext) - gain95(img)))
