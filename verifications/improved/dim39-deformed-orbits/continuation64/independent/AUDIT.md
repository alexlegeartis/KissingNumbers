# Independent exact audit: 39 dimensions, 769,372 points

PASS, 6 November 2026. Certificate SHA256: bfea4089709e11e7b810da37f751868552ed0c39d34256c410e386531f5bc401.

A newly written verifier, not importing or executing the primary or copied direct verifier, checked the complete construction. Run `python audit.py [source-directory]`; Python standard library and g++ are required. The sibling certificate source directory is the default. `relative.cpp` is compiled afresh. Mathematical calculations use exact integers; execution time alone uses floating point.

## Dense layer

All 5,127 integer representatives have 39 nonzero coordinates and squared norms between 99,999,981,539,255 and 100,000,022,283,800. Their norms obey L=9,999,996 and U=10,000,004. The 64 distinct bit masks include identity and are closed under XOR, so they form a coordinate-sign group. Nonzero coordinates make every orbit free. Primitive absolute coordinate tuples and canonical sign masks modulo the group independently exclude duplicate normalized rays among different orbits.

All 841,315,065 relative products were checked: 64 n(n+1)/2 minus n identity diagonal entries. Closure reduces arbitrary group pairs to a relative mask, and symmetry of these real diagonal sign operators permits i<=j. The maximum integer dot is 49,999,015,430,245, at representative indices 228 and 2900 and group index 56. It is below L²/2=49,999,960,000,008. Therefore normalized dot products are at most 1/2. Python checks a conservative bound 2*39*M²<2^63, where maximum coordinate magnitude M=3,268,139, before C++ signed 64-bit accumulation.

## Sparse layer and cross products

The four support edits were applied only after checking their original rows. All 3,424 supports have size eight, are distinct, and their pairwise intersections are at most four, checked directly for all support pairs. Each contributes 128 even-sign vectors. Equal-support distinct sign patterns differ on at least two coordinates. Retained roots are explicitly constructed after deleting the 14 center-leaf pairs: 2,908 distinct roots. All 64 heavy vectors are explicitly constructed and checked distinct. All sparse vectors have squared norm 32. Middle/root/heavy shapes distinguish their families. Heavy-heavy, heavy-root and exact parity-maximized heavy-middle dot maxima are all 16. Root-root and root-middle upper bounds are 16 from their coordinate shapes.

Absolute-value rearrangement bounds check each dense representative against every possible sparse vector, and remain valid under every group sign change. The minimum exact slacks are:
- middle: 2L²-(sum of eight largest magnitudes)² = 2,601,206,588;
- roots: L²-2(sum of two largest magnitudes)² = 53,504,845,470,216;
- heavy: 2L²-(2*largest + next four magnitudes)² = 37,438,665,997,328.

Thus every cross-layer normalized product is at most 1/2. Dense full support also prevents normalized coincidence with any sparse vector. Total: 5,127*64 + 3,424*128 + 2,908 + 64 = 328,128 + 441,244 = 769,372.

## Controls and benchmark

The fresh C++ checker correctly rejected both a duplicate representative and an invalid nonidentity self-orbit product (negative_controls.json). It is not the stale copied verify_direct.py.

The GitHub connector freshly retrieved the previous result at pinned commit 329642cd3fb22f6d4fd8a33c6153dcd00a10c4da: lower_bound=769308, dense_points=328064, sparse_points=441244. Source: https://github.com/alexlegeartis/KissingNumbers/blob/329642cd3fb22f6d4fd8a33c6153dcd00a10c4da/verifications/improved/dim39-deformed-orbits/continuation/orbit_verification.json . This candidate improves that repository benchmark by 64.

A fresh direct retrieval of https://arxiv.org/abs/2611.03682 gives Lin v1's 39-dimensional bound 768,916, below both values. A targeted web search found no stronger relevant result, but this is not an exhaustive priority claim. Publication of this candidate is outside this audit.
