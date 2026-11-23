# Complete Γ72–D_m tails with all horizontal contact profiles

| Dimension | Exact lower bound | Gain over preceding APN construction |
|---:|---:|---:|
| 81 | 6,234,614,752 | 2,426,496 |
| 82 | 6,238,339,364 | 3,195,264 |
| 83 | 6,242,464,108 | 3,887,968 |
| 84 | 6,247,056,392 | 4,516,160 |
| 85 | 6,251,695,496 | 5,623,296 |

The preceding [APN construction](../dim81-85-gamma72-apn-multipartite/) only assigns horizontal phases to D_m rays joining differently labelled vertices. We keep all preceding phase assignments and clipping directions, and assign a distinct unused phase of depth six to each remaining antipodal root ray. Orient each conflict against an earlier ray towards the new ray. Thus the new witnesses cover all m(m−1) antipodal root directions in each tail, with no old points removed.

For every directed clip the phase difference is nonzero B-only and its only horizontal squared distance below six is four, from a singleton raw norm-eight block. The triple profile 4+4+4 uses the established root-class membership test. The present verifier also retains the profiles 4+8 and 6+6. For a horizontal raw norm-twelve point x and raw norm-eight difference u, the opposing point x−u has raw norm twelve exactly when x·u=4. The verifier reconstructs each local q6 representative up to sign from two q4 roots and every q8 class from one q8 vector and the q4 roots with dot four. Exact Gram tests check length, class and uniqueness. It counts a sign-independent mask for each local class and each singleton blocker. For a fixed pair of local blocks all four sign choices survive together. A cyclic permutation of the three Leech blocks permutes the three nonzero first-label colours, so one table applies to all three.

Run `OPENBLAS_NUM_THREADS=1 python verify.py`. NumPy and SymPy are required. BLAS dot tests use binary64 only on integers, with checked worst-case integer sum below 2^53, so no rounding tolerance enters the proof. All new rays and recovered contacts are compatible with the norm-four B-only tail bonus on the same phase space. This is a Nebe-based constructive lower bound, not a non-Nebe construction or an optimality claim.
