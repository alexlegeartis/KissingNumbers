# Exact saturation of four pure 15D vector families

This is a restricted-alphabet statement, not an upper bound on the 15D kissing number. Use the metric diag(1^14,2), squared norm 16 and pairwise inner products at most 8. Consider only these four families:

1. E: all signed four-support vectors (±2)^4 in the first 14 coordinates, with last coordinate 0.
2. H: all vectors (±1^14,±1).
3. L: all vectors with two signed ±2 entries among the first 14 coordinates and last coordinate ±2.
4. A: all 28 coordinate axes ±4e_i for i≤14.

Every compatible configuration in their union has at most 2,564 points, and the verified source attains this. In particular, no exchange that stays entirely inside these families can beat the source, even if individual signs are selected and arbitrary source points are changed.

For E, fix a global sign pattern on the first 14 coordinates. The supports of matching vectors cannot share a triple, so at most C(14,3)/4=91 supports occur. Every E vector matches exactly 2^10 of the 2^14 global sign patterns, giving |E|≤16·91=1,456.

For H with either fixed last-coordinate sign, two heads at Hamming distance d have inner product 16−2d, so the binary code has minimum distance at least 4. The classical bound A(14,4)=512 gives at most 512 heads per sign. See https://aeb.win.tue.nl/codes/binary-1.html.

For L with fixed last-coordinate sign, the tail alone contributes 8 to a mutual inner product, so their signed two-support heads must have nonpositive pairwise products. Fix a signed head coordinate. All vectors containing it must use the same other coordinate, with at most its two possible signs; otherwise their head product is positive. There are 28 signed coordinates and each vector uses two, so at most 28 L vectors occur per tail sign. Finally |A|≤28.

Adding these independent family ceilings gives 1,456+2·512+2·28+28=2,564. The source contains exactly 1,456 E vectors, 512 heads of each sign, 28 L vectors of each sign, and all 28 axes; the package verifier checks these counts and every source pair.

The search pool contains 29,778 points in these pure families and 17,670 outside them. A strict improvement must select at least one outside point. `search.py` prioritizes those 17,670 points as seeds and adds valid per-family source-saturation cuts while still allowing individual pure points and arbitrary source deletions. The 17,670 points are not claimed to exhaust all non-pure possibilities.
