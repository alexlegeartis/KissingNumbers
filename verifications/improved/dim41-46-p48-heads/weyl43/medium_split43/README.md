# 43D medium-body splitting continuation

This exact continuation gives **2,919,803** points, four more than the retained 2,919,799-point Weyl construction in the parent directory.

The certificate removes six individual medium-body points and three old explicit tiny points. It inserts twelve medium-body points with new rational heads and one new tiny-body point. Thus the change is 13 - 6 - 3 = 4. The older certificate is retained unchanged.

Run from this directory with Python 3.9+, NumPy and SymPy:

- `python ../verify.py` replays the inherited configuration.
- `python verify.py` checks the new construction against that base.
- `python verify_independent.py` independently checks all changed interactions and also freshly replays the old independent verifier. Its output directory can be set with `--output`.

`candidate.json` is the exact integer certificate. A medium deletion ID is **495 times its family index plus its row index**, not an index into a flattened list of only medium families. Tiny deletion IDs index the parent's `explicit_points_43.npz` array. Heads have the shared denominator in the certificate.

`medium_families/` supplies all 495 distinct raw body vectors in each of the 60 medium families. Both verifiers check their norm, complement-lattice membership and exact gluing class against the parent's replayed geometry. Completeness uses the inherited population 495. All retained-medium/new-point and retained-tiny/new-point products, all new pairs, and all pure-head interactions are checked directly with integers. The large families use the inherited universal shell bound.

For a raw body h of norm squared 288r and a head z/D, the physical vector before normalization is (h/sqrt(864), z/(D sqrt(3))). Its scaled squared norm is r + |z/D|². The verifier squares positive doubled inner-product inequalities and checks nonpositive products directly. The total dimension is 40 + 3 = 43.

The fresh certificate does not enumerate the large shells. P48n minimum/design, large-family populations, and the stated completeness identities retain exactly the inherited scope of the parent certificate and `../../POPULATIONS.txt`. The parent geometry replay and medium-body checks are fresh. This is a verified repository improvement, not an exhaustive worldwide-priority claim. The manuscript is unchanged.

`verification.json` and `independent_audit.json` retain the successful development receipts. The scripts generate new receipts; the outer package runs them in a temporary copy to preserve shipped hashes.
