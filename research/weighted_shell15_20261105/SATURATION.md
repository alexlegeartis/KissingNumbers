# Exact saturation of four 15D candidate families

This is a restricted-family statement, not an upper bound on the 15D kissing number. Use the metric diag(1^14,2), squared norm 16 and pairwise inner products at most 8. Consider only:

1. E: signed four-support vectors (±2)^4 in the first 14 coordinates, with last coordinate 0.
2. H: vectors (±1^14,±1).
3. L: **arbitrary real vectors** of squared norm 16 with last coordinate +2 or −2.
4. A: the 28 coordinate axes ±4e_i for i≤14.

Every compatible configuration in their union has at most 2,564 points, and the verified source attains this. Therefore an improvement using this union plus further candidate shapes must select at least one of those other shapes. The statement permits individual points and arbitrary changes in every listed family; it does not require intact sign orbits.

For E, fix a global sign pattern on the first 14 coordinates. The supports of matching vectors cannot share a triple, so at most C(14,3)/4=91 supports occur. Every E vector matches exactly 2^10 of the 2^14 global sign patterns, giving |E|≤16·91=1,456.

For H with either fixed last-coordinate sign, two heads at Hamming distance d have inner product 16−2d, so the binary code has minimum distance at least 4. The classical bound A(14,4)=512 gives at most 512 heads per sign. See https://aeb.win.tue.nl/codes/binary-1.html.

For L with fixed last-coordinate sign, the tail alone contributes 8 to a mutual inner product. Thus their 14-dimensional heads have squared norm 8 and pairwise nonpositive inner products. There can be at most 28: after normalizing m heads, their Gram matrix G has rank at most 14 and off-diagonal entries in [−1,0], so tr(G²)≤2m; Cauchy's rank inequality gives m²≤14·tr(G²)≤28m. The bound holds for arbitrary real heads, not only signed two-support heads. Finally |A|≤28.

Adding the family ceilings gives 1,456+2·512+2·28+28=2,564. The source contains exactly 1,456 E vectors, 512 heads of each sign, 28 L vectors of each sign, and all 28 axes; the verifier checks these counts and every source pair.

In the 52,373-point search pool, 47,253 candidates belong to the above saturated families and 5,120 are outside them. A strict improvement must select at least one of those 5,120 within this pool. `search.py` seeds from them and adds valid per-family source-saturation cuts while still allowing arbitrary joint selections from the saturated families and arbitrary source deletions. The 5,120 candidates are not claimed to exhaust all vectors outside these families. In particular, the 17,472 denominator-three vectors are new shapes but all lie in L, so they may participate in coupled replacements but cannot give a gain without an outside candidate.
