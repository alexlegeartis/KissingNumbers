# Exact construction and proof

Let L = 9,999,996. Coordinates are in Euclidean R^39.

## Sparse component

The independently replayed `verify_base.py` constructs:

1. All 128 even-sign vectors with entries ±2 on each of 3,424 eight-coordinate supports, after the four explicit support replacements: 438,272 points.
2. All vectors ±4e_i ±4e_j except the four roots on each of fourteen deleted coordinate pairs: 2,908 points.
3. Four heavy stars, each containing the sixteen vectors ±4e_a + 2 sum(s_b e_b), where the four leaf signs have even parity: 64 points.

Divide each sparse vector by sqrt(32). Every sparse norm is 32. Support intersections, sign parity and the explicit heavy-star checks prove all sparse products at most 16. The sparse count is 441,244. The earlier binary dense component checked by the retained baseline program is not included in the new construction; only its sparse layers are retained.

## Dense component

Let V_1,...,V_2563 be the integer rows in `orbit_certificate.json`. Let G be its set of 128 binary masks. The verifier constructs an independent binary basis and checks that its full XOR span is exactly G. For each row V_i and g in G, include

    x(i,g)_j = (-1)^(g_j) V_i,j / ||V_i||.

The verifier checks L² ≤ ||V_i||² ≤ 10,000,004².

For each row let a_1 ≥ ... ≥ a_39 be its absolute coordinates. The checked inequalities are

    (a_1 + ... + a_8)² ≤ 2 L²,
    2(a_1 + a_2)² ≤ L²,
    (2a_1 + a_2 + a_3 + a_4 + a_5)² ≤ 2 L².

They bound, respectively, its absolute products against every middle vector, every root and every heavy-star vector by half the product of norms. They hold for every sign action, since absolute coordinates are unchanged. Thus all dense/sparse pairs are valid.

For a pair of representatives and a relative group element t, the unnormalized product is

    sum_j (-1)^(t_j) V_i,j V_k,j.

All relative elements occur because G is a group. The only excluded case is i=k and t=0, which compares an indexed point with itself. For every other case the integer verifier proves that the product is at most floor(L²/2). Positive products divided by ||V_i||||V_k|| are therefore at most 1/2. Nonpositive products satisfy the same upper bound automatically.

The binary characters of the coordinate actions turn these sums into a length-128 Walsh transform. Its butterfly consists solely of integer additions and subtractions. Selected full transforms are also compared to direct sums against all group masks, independently of the butterfly computation. The complete number of comparisons is

    (2563 × 2564 / 2) × 128 − 2563 = 420,575,485.

The observed maximum is 49,999,014,167,161, while floor(L²/2) = 49,999,960,000,008.

## Arithmetic and distinctness

Every integer coordinate has absolute value at most 10,000,000. Each norm sum and every partial Walsh sum has absolute value at most 39 × 10,000,000². Every sparse sum before squaring is at most 8 × 10,000,000. The verifier checks bounds placing these operations, their squares where used, and all doubled comparisons below 2^63. NumPy int64 arithmetic is therefore exact, without overflow or tolerance.

Repeated indexed dense points would have product one and contradict a checked distinct-index product bound. Dense/sparse coincidences are excluded by the sparse/dense bound. Sparse distinctness is separately checked in the baseline proof.

Hence the number of distinct unit vectors with pairwise inner products at most 1/2 is

    441244 + 2563 × 128 = 769308.

Placing unit spheres at twice these unit vectors gives pairwise disjoint interiors, each tangent to the central unit sphere. This proves τ(39) ≥ 769,308.
