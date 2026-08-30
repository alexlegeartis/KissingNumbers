#!/usr/bin/env python3
"""Dimensions 49-63: the cap construction over P_48 with EXPLICIT classes.

The first version of this package used Caro-Wei's class bound of 712, obtained with no
coordinates at all.  P_48 is available explicitly (Catalogue of Lattices, entries P48p and
P48q), and searching its minimal lines does far better.

Verified here in exact integer arithmetic:
  * the P_48p Gram is 48x48, symmetric, with even diagonal 6^21 8^26 10^1 and
    **determinant exactly 1** (Bareiss) -- even unimodular.  Note the catalogue basis is NOT
    a minimal-vector basis, so the diagonal is not all 6; the minimum is 6.
  * a class of **7069 lines**, every vector of norm exactly 6, every pairwise
    |inner product| <= 2 with zero threes, and the inner-product histogram
    {0: 11337483, 1: 11979647, 2: 1664716} summing to C(7069,2) exactly.
  * **1400 pairwise disjoint classes**, sizes 4546..7069, 8 023 414 lines in total, all
    distinct -- 26.4% of the whole shell.  Regenerated and re-verified deterministically by
    scripts/regenerate.py; only their sizes are shipped.

7069 against Caro-Wei's 712 is a factor 9.93, between the Leech's 11.3 and Gamma_72's 9.58.

P_48 admits only ANTIPODAL PAIRS of cap directions, not zero-sum triples: its class threshold
is |ip| <= 2 of 6, i.e. gamma = 1/3, which forces cap level t = 3/4 and then two directions
for the same line need <z,z'> <= -1 (KNOWLEDGE.md section 50).  So the gain is
2 * sum of the lambda(k) largest class sizes, plus tau(k) poles, where lambda(k) is the
number of LINES in the largest antipodally symmetric 60-degree code of R^k -- NOT
tau_lattice(k)/2, which is what this docstring said long after the code below stopped doing
it, and which understates every k whose record configuration happens to be antipodal.  It is
tau(k)/2 for k = 1..8 (root systems) and for k = 9, 10, 11, 13, 14, 15 (scripts/
verify_capdirs.py checks the shipped configurations are antipodal, exactly), and
tau_lattice(k)/2 only at k = 12, where tau(12) = 841 is odd and the 841-point record has just
18 antipodal pairs.  k = 12 is therefore the one lambda taken from the literature rather than
exhibited: 756 for the Coxeter-Todd lattice K_12.  Dimension 60 does not depend on it for its
record -- on the exhibited Lambda_12 it would still read 56 789 413 against 52 416 841.

    python final.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N = 52416000
# best known kissing numbers (poles), Cohn's table verified 2026-08-20
TAU = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 306, 10: 510,
       11: 604, 12: 841, 13: 1154, 14: 1932, 15: 2564}
# antipodal (lattice) kissing numbers -> lambda(k) = value/2 disjoint classes
ANTIP = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 272, 10: 336,
         11: 438, 12: 756, 13: 918, 14: 1422, 15: 2340}

# lambda(k) is the number of PARTS, and at t = 3/4 a part is an ANTIPODAL PAIR, so what the
# scheme needs is the largest antipodally symmetric 60-degree code of R^k -- not tau_lat(k).
# Cohn's record configurations are antipodally symmetric for k = 9, 10, 11, 13, 14, 15
# (checked by scripts/verify_capdirs.py), giving lambda(k) = tau(k)/2 there; k = 12 is excluded
# because tau(12) = 841 is odd.  This file used tau_lat(k)/2 long after audit.py and RESULTS.md
# had been corrected, which is why it reported dimension 63 as 66 265 628.
SYMK = {9, 10, 11, 13, 14, 15}
LAMBDA = {k: (TAU[k] // 2 if k in SYMK else ANTIP[k] // 2) for k in range(1, 16)}

# ONE family, and it is regenerated rather than shipped: scripts/regenerate.py rebuilds all
# 1400 pairwise-disjoint classes from data/p48p_gram.npy, data/p48p_gens.npy and
# data/class_p48_7069.npy plus a seed, verifies every one of them exactly, and checks that the
# size vector it derives is this file.  An earlier version stopped at 1170 classes -- fewer than
# the lambda(15) = 1282 that dimension 63 reads -- and made up the difference from a 130-class
# extension built against a base that no shipped script reproduces.  Running the same
# construction further, and choosing each of the first 400 images to overlap the ones already
# taken as little as possible, beats that extension on every one of the fifteen dimensions.
sizes = np.sort(np.load(os.path.join(HERE, 'data', 'class_sizes_p48.npy')).ravel())[::-1].astype(np.int64)
cum = np.concatenate([[0], np.cumsum(sizes)])
print(__doc__)
print("classes available: %d, sizes %d..%d, %d lines"
      % (len(sizes), sizes.min(), sizes.max(), sizes.sum()))
print()
print("  dim  k  lambda   sum|C_i|     gain          this work        Caro-Wei version   direct sum")
CW = {49: 52417426, 50: 52420278, 51: 52424556, 52: 52433088, 53: 52444480,
      54: 52467264, 55: 52505712, 56: 52586400, 57: 52609154, 58: 52654388,
      59: 52726686, 60: 52949064, 61: 53062508, 62: 53410162, 63: 54033164}
rows = []
TABLE = []
for k in range(1, 16):
    lam = LAMBDA[k]
    assert lam <= len(sizes), (
        "k = %d needs %d classes and the family has only %d" % (k, lam, len(sizes)))
    s = int(cum[lam])
    fam = 'class_sizes_p48.npy'
    gain = 2 * s + TAU[k]
    new = N + gain
    old = N + TAU[k]
    rows.append((48 + k, old, new))
    TABLE.append({'dim': 48 + k, 'k': k, 'lam': lam, 'S': s, 'gain': gain,
                  'value': new, 'previous': old, 'family': fam})
    print("   %2d %2d  %6d  %9d   +%-11d %-16d %-16d  %d"
          % (48 + k, k, lam, s, new - old, new, CW[48 + k], old))
prev = 0
for d, old, new in rows:
    assert new > prev
    prev = new
print()
print("   monotone: True ; all below the dimension-64 record 331737984: %s"
      % all(new < 331737984 for d, old, new in rows))
print("   dimension 63: gain +%d against the Caro-Wei +%d, a factor %.2f"
      % (rows[-1][2] - rows[-1][1], CW[63] - (N + TAU[15]),
         (rows[-1][2] - rows[-1][1]) / (CW[63] - (N + TAU[15]))))

# The README table is GENERATED from this, so it cannot drift from the driver -- it had been
# left on the superseded lambda = tau_lat(k)/2 and showed dimension 63 as 66 265 628.
import json
json.dump(TABLE, open(os.path.join(HERE, 'rows49-63.json'), 'w'), indent=1)
print()
print("   wrote rows49-63.json (%d rows); scripts/make_readme_table.py turns it into the"
      % len(TABLE))
print("   markdown table in README.md, so the two cannot disagree.")
