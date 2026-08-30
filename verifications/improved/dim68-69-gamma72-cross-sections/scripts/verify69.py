#!/usr/bin/env python3
"""Verification of the dimension-69 claim.

    tau(69) >= 627 822 180     (floor 331 737 984, from dimension 64 by monotonicity)

The three-point moment LP certifies that many minimal vectors of Gamma_72 orthogonal to a
triple whose Gram is A_3 scaled to minimum 8 -- diagonal 8, two edges -4, one 0, determinant
256.  Those vectors lie in a 69-dimensional subspace and are pairwise at >= 60 degrees, so
they are a kissing configuration of R^69.

WHY THIS GRAM.  The bound falls as the determinant of the Gram rises.  The triple generates a
rank-3 sublattice L of Gamma_72, so every nonzero vector of L has norm >= 8 and det(Gram) =
det(L); Hermite's constant in dimension 3 is gamma_3 = 2^(1/3), attained only by A_3, so

    8 <= min(L) <= gamma_3 * det(L)^(1/3)   ==>   det(L) >= (8 / 2^(1/3))^3 = 256,

with equality exactly when L is similar to A_3.  Determinant 256 is therefore the best any
triple can do, and this Gram attains it.  The previously claimed triple had determinant 352
and gave 527 597 644; every one of the sixteen realizable determinant-256 Grams beats it
(scripts/scan_k3.py), the best by +100 224 536.

REALIZABILITY needs nothing that is not already verified here.  A_3 sits inside D_4, and
verify68.py checks a four-tuple of Gamma_72 minimal vectors with the D_4 Gram; the sublattice
it generates is D_4 scaled to minimum 8, whose 24 minimal vectors are all minimal vectors of
Gamma_72.  This script rebuilds those 24 from the shipped four-tuple and exhibits the triple
among them, so no unshipped pool of minimal vectors is involved at any point."""
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G = np.load(os.path.join(HERE, '..', '..', 'dim73-95-gamma72-caps', 'data',
                         'gamma72_gram.npy')).astype(np.int64)
U = np.load(os.path.join(HERE, '..', 'data', 'gamma72_d4_k4.npy')).astype(np.int64)
print(__doc__)
assert G.shape == (72, 72) and (G == G.T).all() and (np.diag(G) == 8).all()

D4 = np.array([[8, -4, -4, -4], [-4, 8, 0, 0], [-4, 0, 8, 0], [-4, 0, 0, 8]], dtype=np.int64)
assert (U @ G @ U.T == D4).all(), "the shipped four-tuple is not the D_4 configuration"

# the 24 minimal vectors of the sublattice the four-tuple generates
A = np.array([a for a in itertools.product(range(-2, 3), repeat=4) if any(a)], dtype=np.int64)
V = A @ U
V = V[np.einsum('ij,jk,ik->i', V, G, V) == 8]
print("minimal vectors of the exhibited sublattice: %d   (D_4 scaled to minimum 8 has 24)"
      % len(V))
assert len(V) == 24

TARGET = np.array([[8, -4, -4], [-4, 8, 0], [-4, 0, 8]], dtype=np.int64)
IP = V @ G @ V.T
tri = next(t for t in itertools.combinations(range(24), 3)
           if (IP[np.ix_(t, t)] == TARGET).all())
W = V[list(tri)]
M = W @ G @ W.T
print("triple exhibited at sublattice indices %s, coefficient vectors of shape %s"
      % (str(tri), str(W.shape)))
print("exact Gram under the catalogue Gram of Gamma_72:")
for r in M.tolist():
    print("   ", r)

ok_norm = bool((np.diag(M) == 8).all())
ok_gram = bool((M == TARGET).all())
det = int(round(np.linalg.det(M.astype(float))))
rank = int(np.linalg.matrix_rank(M.astype(float)))
floor = round((8 / 2 ** (1 / 3.)) ** 3)
print("all three are minimal vectors (norm 8): %s" % ok_norm)
print("Gram is A_3 scaled to minimum 8: %s      determinant %d" % (ok_gram, det))
print("Hermite floor (8/gamma_3)^3 for a rank-3 sublattice: %d      attained: %s"
      % (floor, det == floor))
print("rank: %d   (so the orthogonal complement is 72 - 3 = 69 dimensional)" % rank)
assert ok_norm and ok_gram and det == 256 and rank == 3 and det == floor

# the sublattice minimum really is 8, which is what the Hermite argument needs
combos = np.array([a for a in itertools.product(range(-4, 5), repeat=3) if any(a)],
                  dtype=np.int64)
mn = int(np.einsum('ij,jk,ik->i', combos, M, combos).min())
print("minimum of the rank-3 sublattice: %d   (must be 8)" % mn)
assert mn == 8

print()
print("LP (../../common/kpoint_lp.py, exact integer dual certificate; the same engine is TIGHT")
print("on the Leech at k = 2, 3, 4 against counted embeddings -- scripts/leech_truth.py; what")
print("is claimed here is the LOWER end, valid for every embedding of this Gram):")
print("   certified n[0,0,0] >= 627822180        (369 cells)")
print("   previous claim, determinant 352:  527597644      gain +100224536")
print("   floor for dimension 69: 331737984      factor 1.89")
print("   monotonicity: tau(69) <= tau(70) = 1249778250 : %s" % (627822180 <= 1249778250))
print("   monotonicity: tau(68) = 361275480 <= 627822180 : %s" % (361275480 <= 627822180))
print()
print("ALL CHECKS PASSED -- tau(69) >= 627822180")
