#!/usr/bin/env python3
"""Dimension 68: the four-point cross-section of Gamma_72, with an exact integer certificate.

    tau(68) >= 361 275 480      (floor 331 737 984, from dimension 64 by monotonicity)

The four-point moment LP with lattice-range constraints certifies that many minimal vectors of
Gamma_72 orthogonal to a four-tuple whose Gram has determinant 1024.  They span a
68-dimensional subspace and are pairwise at >= 60 degrees.

WHY D_4 AND NOTHING ELSE.  The four vectors generate a rank-4 sublattice L of Gamma_72, so
min(L) >= 8 and det(Gram) = det(L).  Hermite's constant gamma_4 = 2^(1/2) is attained only by
D_4, so det(L) >= (8/2^(1/2))^4 = 1024 with equality exactly there -- and the bound falls as
the determinant rises.  So the Gram is forced, not chosen:

    scripts/k4_gram_min_det.py   det >= 1024 by brute force over all 339393 candidate Grams
    scripts/leech_truth.py       the LP is TIGHT on the Leech through k = 4, against counts
    scripts/scan_k4.py           all 104 realizable det-1024 Grams; the certificate is
                                 basis-dependent and ranges over 358774871 .. 361275383
    scripts/verify68.py          the winning four-tuple is exhibited in Gamma_72, exactly

The lattice is forced but the BASIS is not, so the Gram below is the best of the 104, not the
textbook D_4 one; the plain D_4 pattern gives 361184461, and 359925453 was what a single numeric dual
source reached before the certificate took the best of two.

A_4 has determinant 1280 and therefore cannot compete; on the Leech it gives 15540 against
D_4's 17400.  It is kept behind --also-a4 as a comparison, not as a claim.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from kpoint_lp import run3   # the one shared engine, ../../../../common/

M1 = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
N = 6218175600
FLOOR = 331737984
CLAIM = 361275383   # run3's default-path certificate; the CLAIM is
                    # exact_vertex's 361275480, checked by verify68.py
BEST = [[8, 0, -4, 4], [0, 8, -4, 4], [-4, -4, 8, -4], [4, 4, -4, 8]]   # scan_k4.py winner
D4 = [[8, -4, -4, -4], [-4, 8, 0, 0], [-4, 0, 8, 0], [-4, 0, 0, 8]]
A4 = [[8, -4, 0, 0], [-4, 8, -4, 0], [0, -4, 8, -4], [0, 0, -4, 8]]

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--tdes', type=int, default=11, help='design degree (Gamma_72 is an 11-design)')
ap.add_argument('--also-a4', action='store_true', help='also run the A_4 Gram, for comparison')
args = ap.parse_args()

print(__doc__)
rows = ([(BEST, 'best basis (det 1024)', CLAIM), (D4, 'D_4 pattern (det 1024)', 361184461)]
        + ([(A4, 'A_4  (det 1280)', None)] if args.also_a4 else []))
fail = False
for G, name, expect in rows:
    t0 = time.time()
    r = run3(72, N, 8, [1, 2, 3, 4], M1, G, tdes=args.tdes)
    dt = time.time() - t0
    if r is None:
        print("   %-16s NO LP SOLUTION   (%.0fs)" % (name, dt), flush=True)
        fail = fail or expect is not None
        continue
    if r[0] != 'cert':
        print("   %-16s LP %.0f but NO CERTIFICATE   (%.0fs)" % (name, r[1], dt), flush=True)
        fail = fail or expect is not None
        continue
    v, m = r[2], r[3]
    note = ("factor %.2f over the floor" % (v / FLOOR)) if v > FLOOR else "below the floor"
    print("   %-16s certified %-12d %-28s %5d cells   (%.0fs)" % (name, v, note, m, dt),
          flush=True)
    if expect is not None and v != expect:
        print("      MISMATCH: expected %d" % expect)
        fail = True

print()
if fail:
    print("FAILED")
else:
    print("tau(68) >= %d   -- the Gram is at the Hermite floor, so only a larger k can help"
          % CLAIM)
sys.exit(1 if fail else 0)
