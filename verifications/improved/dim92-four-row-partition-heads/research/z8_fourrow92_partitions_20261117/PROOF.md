# 92D four-row three-alphabet construction

**Producer verification PASS; independent audit pending. No publication authorized by this package.**

Dense count29,796,335,616; with the specified disjoint RS/root layers the total is29,797,911,144. This is a complete replacement construction, not a union with the old92D dense or sparse layers.

## Inherited parent

Use the free rank12 Z8 length24 code with homogeneous minimum24 from the frozen `dim88-96-z8-golay` source package. Its generator SHA256 is903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73. The parent minimum is an explicit inherited premise, independently replayed previously by cutoff_audit for the three-basis93–96 package. This verifier does not re-enumerate that parent. The source is the Duursma–Greferath–Litsyn–Schmidt Z8 Golay construction; no exhaustive priority assertion is made.

The first four generator columns form an upper-triangular invertible4×4 block with odd diagonal and are zero in the remaining eight generator rows. Thus prescribing any four-symbol prefix gives exactly8^8 codewords.

## Three alphabets and heads

For symbol a, A(a) is the4-sign affine Boolean truth table with constant bit a_2 and linear bits a_0,a_1. B(a) negates its first physical coordinate. C(a)=2(−1)^a_2 e_(a mod4). All three have identical internal Gram matrices:4 on equal symbols,−4 on symbols differing by4, and0 otherwise. Every cross-alphabet inner product is±2. These finite statements are directly checked.

Delete the first four blocks of each corresponding parent-code copy, retaining20 blocks inR80. For each body color choose a partition of the12 head axes into four3-element groups. A selected4-support is a transversal of those groups. Assign its chosen axis in group j to the corresponding one of the first three cross-polytope symbol axes in prefix row j. Include all16 choices of signs, with head2 times the signed support indicator. Each head has norm16; its within-color Gram equals exactly the Gram of the four removed alphabet blocks. Hence within a color the parent squared-distance minimum96 is preserved.

The selected support counts by color are81,24,6. Every cross-color pair of supports intersects in at most2 axes, verified exactly. Therefore head cross-inner-product is at most8. The retained20 body blocks have cross-inner-product at most40, so the full cross-inner-product is at most48. All vectors have norm96 and pairwise squared distance at least96.

There are111·16=1776 heads, each with8^8 codewords. Prefixes are distinct within each color, and all signed heads are distinct across colors, so the dense count is1776·8^8=29,796,335,616. `head_certificate.json` records every head, prefix and color; the verifier checks all head Gram blocks directly. The earlier105-support two-color candidate is retained in `candidate.json`, but `candidate_three.json` is the stronger frozen result.

## Safe RS6 and root augmentation

Every dense physical coordinate has absolute value≤2. Partition the92 coordinates into one group of size12 and five groups of size16. Over GF16 (modulus x^4+x+1), evaluate every degree≤3 polynomial at six distinct field points0,...,5, restricting its value at0 to the12 labels of the first group. This gives12·16^3 distinct transversal6-supports. Any two supports intersect in at most3 positions, because a nonzero cubic has at most3 roots; all4-column evaluation minors are checked. On every support use amplitude4 and all32 even-parity sign patterns. Same-support sign distance≥2; distinct supports have inner product≤48. Dense cross-inner-product is≤6·4·2=48. This contributes1,572,864 vectors of norm96.

Add √48(±e_i±e_j) for every pair of axes within one RS coordinate group, all four signs. There are4[binom(12,2)+5binom(16,2)]=2664 such roots. They form a subset of scaled D92 roots, hence have norm96 and mutual inner product≤48. An RS support meets any such pair in at most one axis, giving cross-inner-product≤4√48<48. Dense cross-inner-product is likewise≤4√48<48. These layers are disjoint: roots have2 nonzero coordinates, RS vectors6, whereas dense vectors have at least24 (20 body blocks plus4 head coordinates).

Total29,796,335,616+1,572,864+2664=29,797,911,144.

## Replay

Run `python verify.py` with numpy. It checks source binding, prefix uniformity, all alphabets, all1776 signed heads and pair inequalities, GF16 arithmetic/evaluation minors, signs and exact augmentation counts/bounds. The selected-support discovery is not a proof dependency. No old support22 layer is included because the third alphabet can have coordinates±2.
