# Certificate proof of a 769,820-point kissing configuration in dimension 39

Status: both complete primary verification and independent direct-integer audit passed. The independent report is `independent/independent_report.json`. The construction is specified by the frozen certificate, its sparse source files, and exact normalization, not by floating-point optimization output.

## Points and count

1. Dense points: 10,268 explicit integer vectors in 39 coordinates, individually normalized, each acted on by the supplied 32-element group of diagonal sign changes. This gives 328,576 points.
2. Middle points: 3,424 supports of size eight, each carrying all 128 even-parity sign patterns with nonzero entries ±2, divided by √32. This gives 438,272 points.
3. Root points: 2,908 retained vectors with two nonzero entries ±4, divided by √32.
4. Heavy points: 64 specified vectors with one entry ±4 and four entries ±2, divided by √32.

Total: 328,576 + 438,272 + 2,908 + 64 = **769,820**.

## Exact normalization and dense pairs

Every dense integer representative v satisfies L² ≤ v·v ≤ U², where L = 99,999,996 and U = 100,000,004. Its exact unit direction is v/√(v·v).

For every unordered representative pair and every relative sign-group element, except the identity pairing of a representative with itself, the complete integer check establishes

v·(g w) ≤ L²/2 = 4,999,999,600,000,008.

The largest observed product is 4,999,995,789,376,822. There are 1,687,063,204 checked relative products. Consequently every distinct dense pair has normalized inner product at most 1/2. The same checks rule out coincident orbits and nontrivial stabilizers, so the dense point count is exact.

## Dense–sparse pairs

For each actual eight-support, the projected sign group determines the exact maximum signed coordinate sum M over its retained even-parity middle fiber. If an odd projected sign mask exists, M is the absolute-coordinate sum. Otherwise, an odd sign parity requires subtracting twice the smallest absolute coordinate.

All 35,157,632 representative/support checks satisfy M² ≤ 2L². The minimum slack is 11,760,726,551. This proves every dense–middle inner product is at most 1/2.

The universal root inequality is 2(a₁+a₂)² ≤ L², where a₁≥a₂ are the largest absolute coordinates. Its minimum slack is 5,598,758,703,657,294.

The universal heavy-vector inequality is (2a₁+a₂+a₃+a₄+a₅)² ≤ 2L². Its minimum slack is 4,100,012,807,429,007.

These two inequalities cover every retained root and heavy point, including every relevant sign choice.

## Sparse pairs

All unnormalized sparse vectors have squared norm 32. Distinct middle supports intersect in at most four coordinates; within one support, distinct even sign patterns differ in at least two positions. Both facts imply integer inner product at most 16. Root–root and root–middle products obey the same bound directly. The specified heavy–heavy, heavy–root, and heavy–middle products are checked exactly from the sparse certificate.

Thus every distinct unit pair in the complete configuration has inner product at most 1/2. Equivalently, placing unit-ball centers at twice these directions gives 769,820 nonoverlapping unit balls touching a central unit ball.

## Reproduction

Run the included sparse verifier and `verify.py`. The receipt is `orbit_verification.json`. The dense certificate SHA-256 is:

`dbc293625abdbc7f4377ee94e77344967877864ea7717482e34479ce794085df`

All final geometric comparisons use bounded integers. No optimization tolerance, active-constraint screen, or inherited dense-coordinate feasibility is used in the proof.
