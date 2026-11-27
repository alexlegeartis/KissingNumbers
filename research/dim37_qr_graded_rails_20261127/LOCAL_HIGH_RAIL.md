# 37D: a continuous upper rail gap on the central 12-flat

Dimension 37's incumbent is **544,666** at main `5fe68ce190db9760cbf36ee0775a59d775d48730`. Its named historical external comparison is Antoine Olivier's 500,008, using Qihong Lin's code (https://doi.org/10.5281/zenodo.23569153). This note is an obstruction for an extension of fixed layers only, not an upper bound on the general kissing number.

Consider an arbitrary real local vector (x,t), with x supported in the 12 ordinary coordinates U=12,…,23, t on the shared 37th axis, and squared norm 32. Keep the existing ordinary roots (including the roots inside U), both unit marked poles, the QR dense layer, and the 168 existing radical six supports, each with all 32 even body signs and both physical tail signs. The 12D octad refill may be discarded entirely. Then necessarily
```
|t| <= 2 − sqrt(36 − 16 sqrt(5)) < 1.53.
```
In particular no local two-sided shared-body split can make two kissing contacts: it would require t²≥8, whereas this interval forces t²<2.341.

## Exact envelope and proof

Write y=|t|, a_i sorted decreasing ordinary absolute magnitudes, m=a_1 and L=sum a_i. The poles constrain y≤sqrt(8). Ordinary roots give a_1+a_2≤4 and m≤4. The dense layer with marked entry +4 or −4 projects **every** 12-bit sign pattern on U separately in each of the two signed strata, as the companion verifier checks. Its body entries are ±2 and its squared norm 160, so
```
2L + 4y <= 16 sqrt(5).
```
The radical six-support body words have entries ±2 and form the complete even sign code; because their support meets U in at most four places, any local restriction of these signs is possible. Thus for each of their outer footprint masks A,
```
sum_{i in A} a_i + sqrt(2) y <= 8,
```
with a_i here labelled by actual local coordinate, not sorted. Every 5-subset of U contains a triple covered by some such footprint (minimum one, verified), which is used below.

For any nonnegative body with a_1+a_2≤4, if m≤2 then sum a_i²≤2L. If m>2 then all other positions have magnitude at most c=4−m and sum a_i²≤m²+c(L−m). If L<8, this upper expression has maximum at m=2 or m=4, so sum a_i²≤16, impossible since 32−y²≥24. Thus L≥8. In the m>2 case, m²+c(L−m) = 2L−(m−2)(L−2m)≤2L. It follows in all cases that
```
32−y² <= 2L <= 16 sqrt(5)−4y,
```
which excludes the open middle interval between 2−sqrt(36−16 sqrt(5)) and 2+sqrt(36−16 sqrt(5)).

For the remaining high branch we have y>12/5 and y≤sqrt(8). Let D=2L−(32−y²) be the slack; then
```
0 <= D <= y²−4y+16 sqrt(5)−32 <= 16 sqrt(5)−8 sqrt(2)−24 < 0.465.
```
If m≤2, each a_i(2−a_i)≤D. If m>2, since L≥(32−y²)/2≥12 and m≤4, the inequality (m−2)(L−2m)≤D forces m<2.117 and c>1.883. Also the individual deficit a_i(c−a_i)≤D for every i≥2, because their sum is m²+c(L−m)−(32−y²)≤D. In either case all coordinates except possibly the largest are either <0.3 or >1.58: the concave quadratic on [0.3,1.58] exceeds 0.465 at both endpoints already for c≥1.883. At least five coordinates must exceed 1.58, since otherwise the squared norm in x is at most (2.117)²+3·4+8·0.3²<24. Some retained radical footprint meets these five in a triple, of sum >4.74, but requires its total sum at most 8−sqrt(2)y <8−(7/5)(12/5)=4.64. Contradiction.

The result does not exclude weak local rails in the surviving interval, rotating the exterior or removing ordinary roots/radical supports. It does show that assigning a real scalar head to a common local polar body cannot achieve the required antipodal twofold replication under the actual old scaffold.

## Replay

The short companion file independently expands all 180,224 dense parent words, records 4,096 local signs in each marked value, and tests all 792 local five-sets against the literal 168 radical supports. Its optional worktree argument verifies pinned Git blobs. No numerical optimization is a premise.