# A 611,638-point kissing configuration in dimension 38

This construction improves the repository's 611,522-point configuration by 116 points. The original cap/equator/axis proof and both baseline verifiers have been replayed successfully. The new construction and all cross-layer inequalities pass exact integer checks. A separate independent audit also passes, including fresh Golay/Leech reconstruction and all retained grid and equator pairs.

## Construction

Keep the original 6,272 Fano-grid points, 13,710 equator points and 1,596 axes. Replace the 384 cube-edge latitude points by 512 points as follows. Indices below are zero-based.

Choose the eight orthogonal Leech root lines

- 4(e14 − e15), 4(e14 + e15),
- 4(e16 − e17), 4(e16 + e17),
- 4(e18 − e19), 4(e18 + e19),
- 4(e20 − e21), 4(e20 + e21).

For each of the two seven-coordinate tail halves and each sign word y in {−1,+1}^7, let its color be the XOR of the labels i+1 for which y_i is negative. Choose the corresponding root u from the ordered list above. Include both points

(±2u/√3, 16y/√21),

with zeroes in the other tail half. Each has squared norm 128. There are 2 × 128 × 2 = 512 distinct points.

Retain all original caps except the caps on the 14 grid body lines and the eight new body lines. This removes 22 lines, or 132 caps. Relative to the preceding construction, four old cube body lines are restored and six new lines are removed.

The count is

589,548 + 13,710 + 1,596 + 6,272 + 512 = 611,638.

## Exact checks

All pairwise inner products must be at most 64 at squared norm 128.

- New/new: multiplying products by 21 gives 28 u·u′ + 256 y·y′, an integer. All 130,816 pairs pass. The maximum is 1,280/21 < 64. Equivalently, adjacent cube words have distinct orthogonal bodies, while words of the same syndrome differ in at least three coordinates.
- New/grid: body products vanish because the grid uses coordinates 0–13. The maximum doubled-head product is 6. Squaring the remaining radical gives 1,024 × 6² ≤ 4,096 × 21.
- New/retained caps: all retained cap lines differ from the eight new body lines. Their maximum absolute Leech body product is 16, and the maximum sign-head product is 4. The required bound reduces to 2/√3 + 4√(2/21) ≤ √6. After two sign-preserving squarings, this is the integer inequality 25,088 ≤ 30,492.
- New/flat axes: the maximum head product is 4, and 4² ≤ 21.
- New/coordinate axes: 256 ≤ 32 × 21.
- New/equator: the maximum absolute root/equator product is 24. The required bound is 24² ≤ 1,152.
- All old/old interactions are inherited from the replayed source proof. Restored caps were part of that source; their interactions with retained grid points satisfy the original different-body bound because their lines are not grid owners. Their interactions with all new points are explicitly included in the retained-cap check above.

## Files and replay

- `certificate.json` pins source hashes, body roots, cap-line accounting and the point file.
- `new_points.npz` stores integer bodies, sign heads and construction metadata.
- `construct.py --check` regenerates every point and cap owner without search.
- `verify.py` checks every new/new and cross-layer interaction with exact integers.
- `verify_independent.py` independently reconstructs the Golay/Leech inputs, checks retained cap classes, equator/grid pairs and all new interactions. It does not import the first verifier.
- `negative_controls.py` requires rejection of a duplicate point, wrong head norm, missing cap owner and wrong total. The point-file hash is refreshed in geometric corruptions so rejection cannot rely only on checksum mismatch.
- `verification.json`, `independent_verification.json` and `negative_control_results.json` archive the receipts.

From the package root, `python verify.py` checks the complete SHA-256 manifest and runs the historical source, both current proof routes, reconstruction and four negative controls in a temporary directory. `python final.py` derives the count only.

The inherited cap-cap argument uses the source's elementary Golay/Leech minimum-norm proof. Neither route enumerates all pairs of the 98,280 Leech lines or supplies a new formal proof of that lattice theorem. Source attribution and historical records are retained in `../source/` and `../refill/`. This is an improvement over the preceding repository certificate, not a global optimality claim or a fresh comprehensive priority survey. The paper is unchanged.
