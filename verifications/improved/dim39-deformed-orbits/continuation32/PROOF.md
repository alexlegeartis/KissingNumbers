# Exact proof of the 769,436-point construction

Coordinates lie in Euclidean R^39. Put L = 99,999,996 and U = 100,000,004.

## Sparse layer

The retained baseline verifier reconstructs all 3,424 distinct eight-coordinate supports, including the four explicit support replacements. Their intersections have size at most four. Every support contributes its 128 even-sign vectors with entries ±2. Retain 2,908 roots ±4e_i ±4e_j after the fourteen specified root-pair deletions, and the same 64 heavy-star vectors with one ±4 coordinate and four ±2 coordinates of even leaf parity. These 441,244 vectors all have squared norm 32. The baseline proof checks all required sparse products are at most 16 and that the vectors are distinct. Divide them by sqrt(32). No sparse point is deleted or added in this continuation.

## Dense layer and group

The certificate supplies 10,256 integer rows v_i and a set G of 32 coordinate-sign masks. The verifier checks that G is exactly the span of an independently constructed binary basis. For each i and g, include D_g v_i / ||v_i||, where D_g is diagonal with entry (-1)^(g_j). The exact squared-norm checks give L² ≤ ||v_i||² ≤ U².

For every i ≤ k and every t in G, except i=k and t=0, the verifier checks

    v_i^T D_t v_k ≤ floor(L²/2).

The number of products is (10256 × 10257 / 2) × 32 − 10256 = 1,683,122,416. The maximum is 4,999,992,742,814,284, below 4,999,999,600,000,008. Arbitrary orbit pairs reduce to these products because D_g D_h = D_(g XOR h). For positive products the denominator ||v_i||||v_k|| is at least L²; nonpositive products already satisfy the angle bound. Thus every distinct indexed dense pair has normalized product at most 1/2. The checked products also rule out repeated points within or across orbits.

## Actual middle supports and sign parity

For an ordered support S, project every mask in G onto its eight coordinates to obtain P_S. Let E be all 128 even-parity eight-bit masks. A dense-to-middle dot, including the factor 2 in the middle coordinates, is twice a signed sum whose mask is in

    T_S = {e XOR p : e in E, p in P_S}.

If P_S contains an odd mask, T_S is all 256 masks and the exact maximum signed sum is the sum of the eight absolute coordinates. Otherwise T_S is E: the maximum is the absolute sum, minus twice the smallest absolute coordinate exactly when the preferred signs have odd parity. A zero coordinate makes that correction zero, as required.

The verifier computes this maximum M for every row and all 3,424 actual supports. It checks M² ≤ 2L² using integers. This gives 2M/(sqrt(32)||v_i||) ≤ 1/2 for every retained middle point and every group action. There are 35,116,544 support/representative checks. The largest M is 141,421,267, with minimum squared slack 23,640,114,743. Unlike earlier continuations, no bound against nonexistent supports is imposed.

## Root and heavy-star cross products

Let a_1 ≥ ... ≥ a_39 be a row's absolute coordinates. The checked sufficient inequalities

    2(a_1+a_2)² ≤ L²,
    (2a_1+a_2+a_3+a_4+a_5)² ≤ 2L²

bound every root and heavy-star product, including all group sign actions. They retain positive exact slack. These cross bounds also exclude coincidence between sparse and dense points.

## Arithmetic, count and replay

Every coordinate has absolute value at most 100,000,000. Each norm sum and each partial Walsh sum is bounded by 39 × 100,000,000². The largest sparse absolute sum before squaring is at most 8 × 100,000,000. The checker verifies that all these operations and the doubled squared comparisons fit signed 64-bit arithmetic. The primary Walsh transform consists solely of integer additions, subtractions and products, and selected transforms are compared with direct sums. There is no floating tolerance in the proof.

There are 441,244 + 10,256 × 32 = 769,436 distinct unit vectors with all pairwise inner products at most 1/2. Unit spheres centered at twice these vectors therefore form the claimed kissing configuration. The independent audit separately reconstructs the geometry and checks all products; its report explains its different exact arithmetic route.
