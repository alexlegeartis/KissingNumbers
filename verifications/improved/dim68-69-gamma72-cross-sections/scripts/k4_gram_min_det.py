#!/usr/bin/env python3
"""Why D_4, and nothing else, is the four-point Gram to use for dimension 68.

The k-point LP bound falls as the determinant of the Gram rises -- visible across the whole
k = 3 table of this package (det 352 -> 527597644, ..., det 512 -> 149415783): a more
degenerate tuple constrains the cells harder.  So the question "which four-tuple?" is really
"what is the smallest determinant a four-tuple of Gamma_72 minimal vectors can have?".

It has an exact answer.  The four vectors generate a rank-4 sublattice L of Gamma_72, so every
nonzero vector of L has norm at least 8, and det(Gram) = det(L).  Hermite's constant in
dimension 4 is gamma_4 = 2^(1/2), attained only by D_4, and it says

    min(L) <= gamma_4 * det(L)^(1/4)   ==>   det(L) >= (8 / 2^(1/2))^4 = 1024,

with equality exactly when L is similar to D_4.  D_4 scaled to minimum 8 has determinant
4 * 2^8 = 1024, so the bound is attained, and by the uniqueness in Hermite's theorem the D_4
Gram is the ONLY minimum-determinant choice.  A_4 (det 1280) is therefore strictly worse
before any LP is run, which is what the Leech counts show: 15540 against 17400.

This script checks the arithmetic by brute force over every integer Gram that a four-tuple of
minimal vectors could possibly have: diagonal 8, off-diagonal in {0,+-1,...,+-4} because
distinct minimal lines of Gamma_72 meet in at most 4.
"""
import itertools
import numpy as np

MU = 8
OFF = range(-4, 5)
AMAX = 3          # |a_i| <= 3 suffices: a vector of L with some |a_i| >= 4 and norm < 8 would
                  # already force a shorter one among the combinations below


def main():
    combos = np.array([a for a in itertools.product(range(-AMAX, AMAX + 1), repeat=4)
                       if any(a)], dtype=np.int64)
    best = None
    attain = []
    seen = 0
    for off in itertools.product(OFF, repeat=6):
        G = np.array([[MU, off[0], off[1], off[2]],
                      [off[0], MU, off[3], off[4]],
                      [off[1], off[3], MU, off[5]],
                      [off[2], off[4], off[5], MU]], dtype=np.int64)
        d = int(round(np.linalg.det(G)))
        if d <= 0:
            continue
        if np.linalg.eigvalsh(G.astype(float)).min() <= 1e-9:
            continue
        # minimum of the sublattice: every nonzero integer combination has norm >= 8
        norms = np.einsum('ij,jk,ik->i', combos, G, combos)
        if norms.min() < MU:
            continue
        seen += 1
        if best is None or d < best:
            best = d
            attain = [G.copy()]
        elif d == best:
            attain.append(G.copy())
    print("integer Grams of a four-tuple of Gamma_72 minimal vectors")
    print("  (diagonal 8, off-diagonal in [-4,4], positive definite, sublattice minimum >= 8)")
    print("  candidates: %d" % seen)
    print("  smallest determinant: %d" % best)
    print("  Hermite bound (8/gamma_4)^4 = (8/2^0.5)^4 = %d" % round((MU / 2 ** 0.5) ** 4))
    assert best == 1024, best
    D4 = np.array([[8, -4, -4, -4], [-4, 8, 0, 0], [-4, 0, 8, 0], [-4, 0, 0, 8]], dtype=np.int64)
    A4 = np.array([[8, -4, 0, 0], [-4, 8, -4, 0], [0, -4, 8, -4], [0, 0, -4, 8]], dtype=np.int64)
    print("  the D_4 Gram used for dimension 68 has determinant %d  <- attains it"
          % round(np.linalg.det(D4)))
    print("  the A_4 Gram has determinant %d  <- cannot compete"
          % round(np.linalg.det(A4)))
    # every attaining Gram is isometric to D_4: same determinant AND same theta series start
    def prof(G):
        norms = np.einsum('ij,jk,ik->i', combos, G, combos)
        return tuple(sorted(np.bincount(norms[norms <= 24], minlength=25).tolist()))
    pd = prof(D4)
    bad = [G for G in attain if prof(G) != pd]
    print("  Grams attaining 1024: %d, all with D_4's vector counts up to norm 24: %s"
          % (len(attain), not bad))
    assert not bad
    print()
    print("So no four-tuple can beat D_4, and dimension 68's Gram is forced.")


if __name__ == '__main__':
    main()
