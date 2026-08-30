#!/usr/bin/env python3
"""Verification of the dimension-68 claim.

    tau(68) >= 361 275 480     (floor 331 737 984, from dimension 64 by monotonicity)

The four-point moment LP certifies that many minimal vectors of Gamma_72 orthogonal to a
four-tuple of minimal vectors whose Gram has determinant 1024.  Those vectors span a
68-dimensional subspace and are pairwise at >= 60 degrees, so they are a kissing configuration
of R^68.

WHICH FOUR-TUPLE, AND WHY.  The bound falls as det(Gram) rises.  The four vectors generate a
rank-4 sublattice L of Gamma_72, so min(L) >= 8 and det(Gram) = det(L); Hermite's constant in
dimension 4 is gamma_4 = 2^(1/2), attained only by D_4, so

    8 <= min(L) <= gamma_4 * det(L)^(1/4)   ==>   det(L) >= (8/2^(1/2))^4 = 1024,

with equality exactly when L is a scaled D_4.  scripts/k4_gram_min_det.py confirms the floor by
brute force over all 339393 integer Grams a four-tuple of minimal lines could have.  So the
LATTICE is forced, and so is the BOUND: the basis costs nothing.  data/scan_k4.json records a
different certificate for almost every one of the 104 determinant-1024 Grams, 358 774 871 up
to 361 275 383, and that 0.7% spread used to be read here as basis-dependence among isometric
Grams.  It is not.  It is how good a dual the numeric solve hands to the certifier:
scripts/presentation_invariance.py runs the k = 3 analogue, where the sixteen realizable
determinant-256 Grams are sixteen bases of one scaled A_3, and the thirteen distinct
certificates of the default path collapse to ONE, 627 822 180, as soon as more than the first
dual is kept.  Every one of those numbers is a valid lower bound either way -- the certifier
repairs any dual into exact feasibility -- so the weak ones are weak witnesses, not errors.
scripts/scan_k4.py runs all 104 and this script checks the best of them.

REALIZABILITY is by exhibition, and needs no pool of minimal vectors.  data/gamma72_d4_k4.npy
is a four-tuple of Gamma_72 minimal vectors with the standard D_4 Gram; the sublattice it
generates is D_4 scaled to minimum 8, whose 24 minimal vectors are all minimal vectors of
Gamma_72, and all 104 determinant-1024 Grams occur among their four-tuples.  This script
rebuilds those 24 and exhibits the winning four-tuple among them."""
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G = np.load(os.path.join(HERE, '..', '..', 'dim73-95-gamma72-caps', 'data',
                         'gamma72_gram.npy')).astype(np.int64)
U = np.load(os.path.join(HERE, '..', 'data', 'gamma72_d4_k4.npy')).astype(np.int64)
print(__doc__)
assert G.shape == (72, 72) and (G == G.T).all() and (np.diag(G) == 8).all()

D4 = np.array([[8, -4, -4, -4], [-4, 8, 0, 0], [-4, 0, 8, 0], [-4, 0, 0, 8]], dtype=np.int64)
M0 = U @ G @ U.T
print("the shipped four-tuple, shape %s, has the standard D_4 Gram: %s"
      % (str(U.shape), bool((M0 == D4).all())))
assert (M0 == D4).all()

# the 24 minimal vectors of the sublattice it generates
A = np.array([a for a in itertools.product(range(-2, 3), repeat=4) if any(a)], dtype=np.int64)
V = A @ U
V = V[np.einsum('ij,jk,ik->i', V, G, V) == 8]
print("minimal vectors of that sublattice: %d   (D_4 scaled to minimum 8 has 24)" % len(V))
assert len(V) == 24

# the winning basis, from scripts/scan_k4.py
BEST = (0, -4, 4, -4, 4, -4)
TARGET = np.array([[8, BEST[0], BEST[1], BEST[2]],
                   [BEST[0], 8, BEST[3], BEST[4]],
                   [BEST[1], BEST[3], 8, BEST[5]],
                   [BEST[2], BEST[4], BEST[5], 8]], dtype=np.int64)
IP = V @ G @ V.T
quad = next(t for t in itertools.permutations(range(24), 4)
            if (IP[np.ix_(t, t)] == TARGET).all())
W = V[list(quad)]
M = W @ G @ W.T
print("winning four-tuple exhibited at sublattice indices %s, shape %s"
      % (str(quad), str(W.shape)))
print("exact Gram under the catalogue Gram of Gamma_72:")
for r in M.tolist():
    print("   ", r)

ok_norm = bool((np.diag(M) == 8).all())
ok_gram = bool((M == TARGET).all())
det = int(round(np.linalg.det(M.astype(float))))
rank = int(np.linalg.matrix_rank(M.astype(float)))
floor = round((8 / 2 ** 0.5) ** 4)
print("all four are minimal vectors (norm 8): %s" % ok_norm)
print("Gram is the scanned winner exactly: %s      determinant %d" % (ok_gram, det))
print("Hermite floor (8/gamma_4)^4 for a rank-4 sublattice: %d      attained: %s"
      % (floor, det == floor))
print("   so no four-tuple can do better -- scripts/k4_gram_min_det.py checks this by brute")
print("   force over all 339393 integer Grams a four-tuple of minimal lines could have")
print("rank: %d   (so the orthogonal complement is 72 - 4 = 68 dimensional)" % rank)
assert ok_norm and ok_gram and det == 1024 and rank == 4 and det == floor

# the sublattice minimum really is 8, which is what the Hermite argument needs
combos = np.array([a for a in itertools.product(range(-3, 4), repeat=4) if any(a)],
                  dtype=np.int64)
mn = int(np.einsum('ij,jk,ik->i', combos, M, combos).min())
print("minimum of the rank-4 sublattice: %d   (must be 8)" % mn)
assert mn == 8

# the scan, as recorded
sf = os.path.join(HERE, '..', 'data', 'scan_k4.json')
if os.path.exists(sf):
    sc = json.load(open(sf))
    vals = [v for v in sc.values() if v is not None]
    none = sum(1 for v in sc.values() if v is None)
    print("scan_k4.json: %d Grams, %d certified, %d that HiGHS does not solve (all 104 have"
          % (len(sc), len(vals), none))
    print("   the same 1681 cells and 3473 rows, so those six are hard programmes, not big ones;")
    print("   they return status 1, time limit reached, at 900 s and at 4000 s)")
    print("   range %d .. %d, best %d" % (min(vals), max(vals), max(vals)))
    assert max(vals) == 361275383, max(vals)
    assert sc.get(json.dumps(list(BEST))) == 361275383
    print("   those are run3's DEFAULT-PATH certificates, from a rounded floating-point")
    print("   dual; the claim below does not come from them -- see the exact vertex.")

print()
print("THE EXACT OPTIMUM OF THE PROGRAMME (../../common/exact_vertex.py).  run3 proposes a")
print("dual with a floating-point solve, rounds it and repairs it; that is sound but lossy,")
print("and 361 275 383 is what it repairs into.  Solving for the dual VERTEX instead -- the")
print("identifications substituted out exactly, the active system solved in Fractions, and")
print("A^T y <= c then verified over every column -- gives the programme's true optimum,")
print("which the numeric primal matches, so no dual can do better:")
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'common'))
from exact_vertex import exact_optimum   # noqa: E402
M1G = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
_t = time.time()
_v, _rep = exact_optimum(72, 6218175600, 8, [1, 2, 3, 4], M1G, D4.tolist())
print("   %s" % _rep)
print("   certified n[0,0,0,0] >= %d      (%.0f s)" % (_v, time.time() - _t))
assert _v == 361275480, _v
print("   floor for dimension 68: 331737984      factor %.2f" % (361275480 / 331737984.))
print("   monotonicity: tau(68) <= tau(69) = 627822180 : %s" % (361275480 <= 627822180))
print("   the default path certifies 361275383 for the same programme, 97 lower; both are")
print("   valid lower bounds and the exact one is claimed.")
print()
print("ALL CHECKS PASSED -- tau(68) >= 361275480")
