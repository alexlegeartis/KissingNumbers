# Free-head completion in dimensions 26 and 27

Date: 26 October 2026.

The frozen candidates below extend the published shared-midpoint constructions. Construction-side integer checks passed; the separate independent-verification directory records the authoritative replay status. No novelty claim is made for the elementary construction principle.

## General midpoint lemma

Work at squared norm 32, with allowed inner product at most 16. Let u and v be two original Leech minimal vectors with u dot v = 0. Their midpoint h = (u+v)/2 has squared norm 16 and clears every original Leech shell point, without deleting either endpoint:

- h dot u = h dot v = 16;
- for every other shell point z, both endpoint products are at most 16, so their average is at most 16.

If t is a tail vector of squared norm 16, then (h,t) has squared norm 32. It remains only to check its products with the earlier inserted points, other new points, and pure tail points.

## Existing coordinate split

The baseline certificates use octad coordinates I = {0,2,4,5,6,10,11,23}; J is the complementary set of 16 coordinates. An original residual fibre consists of r plus 16 normal vectors, with r supported on J and squared norm 16. These normal vectors form an 8-dimensional cross-polytope, also of squared norm 16. Opposite fibre points have inner product zero, so each r satisfies the midpoint lemma.

The 4,320 residual fibres group into 135 normal frames, each with 32 residual vectors forming a 16-dimensional cross-polytope. Distinct frame pairs have maximum residual inner product 4 or 8. These facts were regenerated and checked exactly from the Golay description of the Leech shell, rather than assumed from the old search.

Every old inserted point has all normal-octad coordinates and all tail coordinates bounded in absolute value by 2. This bound does NOT hold for all residual coordinates; the failed global-coordinate shortcut must not be used.

## Dimension 27: 288 added points

First choose one occupied residual frame for each tail axis such that its old midpoint group uses the perpendicular tail plane. The frozen candidate uses frame 3 for x, frame 0 for y, and frame 102 for z. For each of its 32 residual vectors r, add (r,+4 e_a) and (r,-4 e_a). This gives 192 points.

For a different old frame, the head product is at most 8 and the tail product at most 8. In the same old frame, the chosen tail is perpendicular to the old tail plane; the head product is at most 16. Thus all earlier insertions are cleared. For one fixed signed tail, the 32 residual heads have pairwise product at most zero. Different tail axes are orthogonal, and opposite tails have negative product. Hence the 192 added points are mutually compatible.

Next, for every i in I and each tail coordinate a, add all four sign choices

    (±4 e_i, ±4 e_a).

There are 8*3*4 = 96 such points. Their head parts clear the full Leech shell because shell coordinates have absolute value at most 4. Their products with earlier insertions are at most 4*2+4*2 = 16. They are mutually compatible because they are distinct norm-32 roots of a scaled D lattice. Their head parts are orthogonal to the 192 residual heads, so their cross products are at most the tail product 16.

All 288 additions clear the retained 12 pure tail cuboctahedron points (±4,±4,0) and permutations. The resulting count is

    196560 - 8064 + 16128 + 288 + 12 = 204924.

Frozen extension: extension27_204924.npz.
SHA-256: 75729558f81bebedcaf50394aaccde6a86ffa62d27966f2a8c41fc41b19d166e.
It contains only the 288 additions and the baseline certificate hash. The older 192-point extension remains unchanged as a superseded checkpoint.

## Dimension 26: 192 additions and a pure-tail exchange

Choose a residual frame not used by the 34 old midpoint groups; the frozen construction uses frame 0. Its 32 residual heads, together with the 16 normal coordinate heads ±4 e_i for i in I, form 48 points with pairwise nonpositive products in 24 dimensions. For each of these heads, use four tails ±4 e_x and ±4 e_y. This gives 192 points.

The same midpoint, coordinate, and cross-frame arguments clear all old non-pure points. Replace the old six pure hexagon points by the four square points (±4,±4). These four points clear both the old insertions and the 192 additions. Thus the gain is 192-6+4 = 190, not 192, and the count is

    196560 - 4352 + 8704 + 192 + 4 = 201108.

Frozen extension: extension26_201108.npz.
SHA-256: 8499ec82f0c9633b0f2439fb14a298cbe26643ceb4a6e0c6c0ee9cfa73319309.
The file contains 192 additions, four replacement pure points, the selected frame, and the baseline hash. The old six pure points are not retained.

## Saturation of this fixed head-norm class

In 27D, retaining all 12 pure D3 tail points forces |t_a|+|t_b| <= 4 for each coordinate pair. If a new head has squared norm 16, its tail also has squared norm 16. Write a>=b>=c>=0 for the absolute tail coordinates. Then a+b<=4 and

    a^2+b^2+c^2 <= (4-b)^2+2b^2 <= 16,

with equality only at b=c=0 and a=4. Thus its tail is one of the six signed coordinate axes. Within a single signed-tail class, distinct heads must have nonpositive inner products. The orthoplex bound in 24 dimensions permits at most 48 such heads, hence at most 6*48 = 288 additions. Our candidate attains this number.

With the four pure square points in 26D, the analogous inequality |t_x|+|t_y|<=4 forces four signed-axis tails for head norm 16, so the corresponding upper limit is 4*48=192, also attained.

This proves saturation only for additional points with head squared norm 16 and the stated retained pure-tail set. It does NOT prove global optimality or rule out changing head norms, sacrificing earlier insertions or pure points, or rebuilding the construction.

## Inputs and provenance

Baseline 27D certificate SHA-256:
592ba5d22deebeb43f15ab3600d8ab75335c3b0591ab5ebbd7a098a0fa1e54ad.

Baseline 26D certificate SHA-256:
6a1fcba8c37657cd6111bb2954c354e5a2fc3668ade9c6551c5c13db3156ca43.

Both were fetched from the current default branch and byte-checked against the repository's recorded certificate hashes:
https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim26-27-midpoint-replacements/source

The existing shell generator was read for provenance. Our discovery scripts independently implement its elementary Golay description; independent verification uses separately written checks. The NPZ files contain only numeric arrays and non-secret metadata, loaded with allow_pickle=False.

## Reusable completion for future group selections

For either dimension, a simpler uniform completion is available whenever at least one of the 135 frames is unused. Choose all 32 residual vectors from that unused frame, and all 16 signed normal coordinate vectors ±4 e_i. These 48 heads form a 24-dimensional cross-polytope. Pair them with every signed tail axis: four directions in 26D, six in 27D. Every old group is in a different frame, so all residual-head/old-insertion products are at most 8+8=16. The normal-coordinate argument is unchanged.

Thus every compatible two-axis selection on fewer than 135 distinct frames receives the same free bonus: 192 additions with four pure square points in 26D, or 288 additions with twelve pure cuboctahedron points in 27D. Optimizing the old group selection therefore does not require preserving the particular three frames used by the frozen 27D candidate. A future candidate may select a fresh unused frame and reverify its coordinates.

Each distinct frame has distinct residual fibres, so its deleted equator points are disjoint from those of every other frame. A group deletes 4*32=128 points and inserts 8*32=256 distinct points, for net gain 128. Our complete two-axis label pool contains 135*28*2*3=22,680 labels, with at most one selected label per frame. The search enforces that rule explicitly, rather than count overlapping groups as disjoint. This is a search space, not an optimality theorem.
