# 37D: all equal-amplitude dense bodies with arbitrary real scalar head are blocked

Live comparator 544,666 at `5fe68ce190db9760cbf36ee0775a59d775d48730`. Dated external comparison, Antoine Olivier with Qihong Lin's support code, 500,008 (https://doi.org/10.5281/zenodo.23569153). The exact incumbent dense list is that of `verifications/improved/dim37-qr-flat-repack/`, derived by QR compression from the Yang et al. parent (https://arxiv.org/html/2609.35051v1, Appendix B, Tables 11 and 12). The literal data blob `dense_parent38.json` is `3343690a0bd2ebdf7b0c30fd80bee29b0fe8702a`.

**Statement.** Keeping only the incumbent 180,224 dense QR words and the two marked-axis poles already prevents adding a point with equal nonzero magnitudes on *all* 36 ordinary sites, whatever its real marked coordinate. No assumption of a new linear code, fixed head list or phase uniformity is needed. The zero ordinary body case is just the two existing poles. Changing some ordinary magnitudes or jointly moving/deleting the old dense rows remains open. This is not an upper bound on the general kissing number.

## One quotient gives all four signed buckets

Use ordinary entries ±2 by choosing the common scale, and write z for the real marked coordinate. Keep the parent's old squared-norm and physical marked values for grades e=0,1,2,3 respectively (+4,+1,−4,−1). Grade is given by the erased old coordinates 36 and 37, in increasing bit order. The original length-38 parent has 11 represented cosets of a common rank-14 binary kernel H. All 11 representatives have grade zero. The zero-grade kernel H0 is rank 12 and lives on the 36 body coordinates; the grade map is onto. Choose H lifts u and v of grade 1 and 2. Quotient the 36 body coordinates by H0 to get 2^24 residues. The 11 seeds for the +4 grade are the projected original representatives; the other grades translate this entire seed set by respective residues u,v,u+v. The distance in the syndrome Cayley graph with the 36 unit-coordinate moves is exactly the minimum ordinary Hamming distance to the relevant old code, because any path encodes a sign difference of its bit weight and conversely any minimum-weight sign difference gives such a path. Therefore a single multi-source BFS yields four per-grade distances for every possible candidate body sign choice.

For each coset q, write D(q) for its distance to grade zero. The four distances are D(q),D(q+u),D(q+v),D(q+u+v); addition is XOR. The exact checked grade-zero histogram has radius 10 with only two residues at 10 (`12518842` gives `(10,1,9,10)` and `14022956` gives `(10,9,1,10)`, with the first coordinate corresponding to the +4 grade), and **no** four-tuple simultaneously dominates (10,9,6,8) in this order. The minimum across all four grades is at most 7, checked also by independently dropping the two grades and using the 22-bit quotient by projected H. Replacing q by q+v exchanges the two opposite signs for both magnitudes, so it suffices to reason z≥0 without assuming spatial antipodal symmetry of the old section.

## Continuous angular proof

The two poles force |z|≤sqrt(48), since the new norm squared is 144+z². The raw dot with an old grade m, squared norm 144+m², at ordinary Hamming distance d is `144−8d+zm`. It may not exceed `sqrt((144+z²)(144+m²))/2`.

- If 0≤z≤2, some old row has d≤7 and then the raw dot is ≥88−8=80, exceeding the largest possible threshold sqrt(148×160)/2<80.
- If 2<z≤4, the threshold with any old grade is at most 80. For m=+4, d≤9 gives raw dot ≥72+4z>80. For m=+1, d≤8 gives ≥80+z>80. For m=−4, d≤5 gives ≥104−16>80. For m=−1, d≤7 gives ≥88−4>80. Thus the necessary distances are (10,9,6,8), already excluded by the quotient.
- If 4<z≤sqrt(48), some m=+4 row has d≤10, because the individual grade has covering radius 10. At d=10, comparing (64+4z)² with 40(144+z²) gives their difference (z−4)(416−24z)>0; a smaller d is worse. Hence this range is excluded too.

The standalone exact quotient verifier regenerates all data from the literal pinned 11 cosets and 14 generators, tests rank and surjectivity of the grade map, computes the 2^24 distances using exact integer arithmetic, and replays histogram/counts. Node 18+; main memory about 85 MB. Optional worktree argument checks file identity.

Companion exact results in this package:

- [Four-bank exchange with no gain](FOUR_BANK_TRADE.md) ([replay](verify_trade.mjs)).
- [Complete one-center menu on the five natural rails](ONE_CENTER_MENU.md) ([replay](verify_menu.mjs)).
- [Real high-rail obstruction on the central 12-flat](LOCAL_HIGH_RAIL.md) ([replay](verify_rail.mjs)).
- [Replay for the arbitrary real scalar head obstruction above](verify_flat.mjs).
