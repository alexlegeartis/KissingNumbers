# Polynomial-image support search in 32D

This is an actual construction search, not a new bound. It starts from the 352,800-point structural seed in `seed32.json`, which `verify_seed.py` independently checks via support intersections and the minimum distance of its binary sign code. The generated vectors have rational coordinates before normalization.

The finite-field engine enumerates coefficient grids in GF(32), rejects non-permutation polynomials, maps affine 3-flats under the remaining maps, translates and deduplicates their images, and retains supports meeting a chosen conflict budget. The exchange search chooses individual signed antipodal lines, including both sign parities, and may remove individual lines from full seed families. It does not assume whole support orbits must be kept or deleted. A numerical MILP proposes exchanges; it is not an upper-bound proof.

From the repository root:

    python research/gf32_polynomial_octads_20261105/verify_seed.py
    python research/gf32_polynomial_octads_20261105/run.py --support-limit 11 --seconds 120
    python research/gf32_polynomial_octads_20261105/verify_exchange.py research/gf32_polynomial_octads_20261105/proposal.json

The verifier proves any resulting count by checking the unchanged source, every proposed new line, and all necessary source deletions. It does not infer a record from a solver status. A run over all 11 supports in this particular coefficient grid found no gain; this is a solver result, not an exclusion proof. Other polynomial families, partial non-antipodal signs, and other source constructions are not covered.
