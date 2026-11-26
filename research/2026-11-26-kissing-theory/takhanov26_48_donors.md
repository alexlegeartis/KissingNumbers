# Four chambers, 48 exact donors, four paid old rows (26D)

This is a constructive **reduction**, not a new bound. It replaces the parity quartets with simultaneous rotated-octahedral modules in **all four chambers** of [Takhanov v0.2](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2). The exact source comparison is 201,212.

Use the ordered six-coordinate charts

| centre | three regular axes | original tail | alternate tail |
|---|---|---|---|
| 0 | 1,2,3 | 9 | 8 |
| 0 | 3,5,7 | 8 | 9 |
| 4 | 1,2,6 | 8 | 9 |
| 4 | 5,6,7 | 9 | 8 |

Put 12 rows in each chart: the centre is ±531, independently the regular triple runs through ±253 times the three mutually orthogonal triples (1,2,2), (2,1,−2), (2,−2,1), and the two tail numerators are fixed +519,+21. Divide by 1062. The squared numerator norm is 1,127,844=1062². All rows obey the signed axis caps and satisfy the [continuous six-coordinate exterior gate](takhanov_exchange_26.md) because the largest sum of regular head magnitudes plus twice the alternate tail is 1838<2124. With equal centre signs in the same chamber and different orthogonal triples the dot numerator is 551,763<563,922; with opposite centre signs the greatest numerator, attained when triples agree, is 563,922 exactly. Among different chambers, aligned two head positions contribute at most 531²+506² (or 2·506² if there is no common centre). With swapped tails add 2·519·21. With unswapped tails the two opposite chambers have no common head. These bounds are all within 563,922. The exact replay checks all 1,128 donor pairs.

Their only paid old skeleton rows are the four rows of [0,4,8,9] with both tails positive, regardless of signs at 0 and 4. This follows from the integer strict threshold 1062 in the mixed supporting dot product, or by exact replay of all 224 skeleton rows. The 200,864 exterior points, 20 signed axes and all 220 remaining skeleton rows are retained.

One may now rebuild each original five-coordinate chamber from scratch. For every rebuilt row, impose unit norm, coordinate caps 1/2 and the local cutoff
`531 |u| + 519 v + 253 max(|r1+2r2+2r3|, |2r1+r2−2r3|, |2r1−2r2+r3|) <= 531`.
Here u is the centre coordinate, v the original tail; the row has no alternate-tail coordinate. Different rebuilt chambers and the old tetrads meet a rebuilt support on at most two coordinates, so the caps suffice. A rebuilt row and a donor in another chamber share coordinates whose donor numerator absolute sum is at most 1058<1062; hence they also are safe, regardless of the rebuilt row's angles. Within each rebuilt chamber enforce pairwise dot≤1/2 as usual.

After all 48 donors are added the fixed population is **201,152**. Any total **61** rebuilt rows gives 201,213, beating Takhanov. In particular a single 16-row chamber code copied into all four chambers gives **201,216**. The 16-row requirement is ruled out for these donors in the analytic obstruction below. This gate is materially different from the [24-row balanced-quartet gate](takhanov26_balanced_quartet.md) and from the earlier [21-row single-octahedron gate](four_octahedra_21.md). The original 26-row template and mere signed rotations are not a proof of this new gate; the task is to produce actual new chamber geometry or tune the donor architecture jointly. [Exact finite checks](verify_48_octahedral.py) use Takhanov's pinned chart file.

## Analytic obstruction: at most 15 chamber rows in each of two exact versions

The conditional 16-row target above is **now excluded for the specified rational donors**. The proof also excludes the mirrored 10089-norm version: in each chamber take the six regular vectors \(\pm24(1,2,2),\pm24(2,1,-2),\pm24(2,-2,1)\), centre +50 and tails +49,+2, and take their reflection in the first regular coordinate with centre −50 and the same tails. All divided by \(\sqrt{10089}\). Within chambers the orthogonal groups have maximum numerator 4905; across the two hot signs use \(Q^TDQ=9I-2(1,2,2)(1,2,2)^T\), maximum absolute entry 8, and the negative core −95, giving 4513. Across chambers the previous cap estimate gives at most 5000. Exactly the same four old skeleton rows are paid. These checks are repeated with integer arithmetic in the replay.

For each version let A,B,k denote the normalized positive amplitudes on the centre, old tail and tripod factor. It is necessary that a candidate row obey
`A |u| + B v + k M(r') <= 1/2`,
choosing the tighter centre sign; r' is r or its one-coordinate reflection for the corresponding hot sign. The orthogonal triples have length 3, so \(M(r')\ge\sqrt3\|r\|\). If v≥0 then \(|u|+v\ge1/2\) and \(\|r\|\ge1/\sqrt2\) by norm and caps, giving score at least \(B/2+k\sqrt{3/2}>1/2\) in both versions. Thus v<0; let t=−v. The caps give each of the other four magnitudes ≥\(\sqrt{1/4-t^2}\), and \(\|r\|\ge\sqrt{3/4-t^2}\). The resulting lower bound \(A\sqrt{1/4-t^2}-Bt+k\sqrt{3(3/4-t^2)}\) exceeds 1/2 at t=1/6 and decreases with t. In each sign layer of u there are six regular sign patterns aligning (up to simultaneous negation) with a tripod form; for those, \(M(r')\ge5\sqrt{1/4-t^2}\). The resulting \((A+5k)\sqrt{1/4-t^2}-Bt\) exceeds 1/2 at t=3/10. Therefore all rows have t>1/6 and these twelve high-type sign patterns require t>3/10.

At most one point fits each full sign orthant (aligned capped magnitudes give dot≥3/4 as in the [2513 proof](takhanov26_2513_23.md)). Since v is always negative, any code of size 16 uses all sign vertices of Q₄ on (u,r1,r2,r3). Every one of those four magnitudes is therefore at least \(c=1/(2\sqrt2)\) by the directed adjacent-orthant rule. Set \(d_i=1/2-|x_i|\) for the four sign coordinates; norm gives \(t^2=\sum_i(d_i-d_i^2)\). On an edge of Q₄ flipping sign i, the squared separation requires
`2(d_i+d'_i) − (d_i+d'_i)^2 <= Σ_(j≠i)(d_j−d'_j)^2 + (t−t')^2`.
Sum on all 32 edges. The left side is ≥\(2\sum t^2\). The four-sign-coordinate part on the right is ≤\(3\sum d_i^2<0.516\sum t^2\) because \(d_i/(1-d_i)\le(1/2-c)/(1/2+c)<0.172\).

For a high–high edge, \((t-t')^2/(t^2+t'^2)\le2/17\); the ratio is at most 2/5 on every edge, by the respective intervals [3/10,1/2] and [1/6,1/2]. The two low sign words in each regular cube are antipodal, so every high vertex has exactly one low regular-direction neighbor and at most one more through the u-sign edge. If L and H sum t² over the four low and twelve high rows respectively, the tail part is at most \(1.6L+(4/5+4/17)H\). The full right side is therefore ≤\(2.116L+1.552H<2(L+H)\), as L≤1 and H>12(3/10)²=1.08. Contradiction.

Each of the four chambers has at most 15 rows, even when fitted independently. Accordingly both specified 48-donor architectures have total ≤201,212 and **cannot improve** Takhanov. The one-tripod 24-donor architecture still has the live 21-row condition and is not covered. The result does not assert that 15 is attained or exclude different donor shapes.
