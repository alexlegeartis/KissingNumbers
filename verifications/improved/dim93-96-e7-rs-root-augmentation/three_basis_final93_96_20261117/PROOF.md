# Strongest three-basis constructions with safe RS16 and root layers,93–96D

The included `base/PROOF.md` proves the three-alphabet Z8 construction and compressed93/94D configurations. Its independent audit includes a fresh replay of the inherited parent minimum24. The included `e7/PROOF.md` supplies the explicit100-head95D improvement; independent audit of that new head table is requested separately. No old sparse22 or support24 layers are retained.

## Orthogonal coordinates and eligible cells

In93/94D use88 bodycoordinates and the displayed five/six headcoordinates. In96D use the96 physical rowcoordinates. Every dense coordinate has absolute value at most2.

For95D use88 bodycoordinates, followed by the orthonormal E7 head coordinates sqrt(2)*x0,x1,x2,x3,x5,x6,x7. The first headcoordinate (ambient index88) may have absolute value2sqrt(2); all others are at most2. EXCLUDE index88 from the Reed–Solomon support cells. Thus the eligible cell counts are93,94,94,96 for dimensions93,94,95,96 respectively.

Partition eligible cells, in increasing order, into one initial group of size r and five groups of size16, where r=13,14,14,16 respectively. Exact groups are recorded in augmentation.json.

## Six-support Reed–Solomon layer

OverGF16=F2[x]/(x^4+x+1), evaluate all cubic polynomials at the six distinct field elements0,...,5. Restrict the constant coefficient to0,...,r-1; the other three coefficients are arbitrary. Each evaluated word chooses one coordinate in each of the six groups. Distinct cubics agree at at most three positions, so distinct supports intersect at most three times. The verifier checks field arithmetic, all15 four-column ranks of the evaluation generator, and enumerates every support to prove its count and injectivity.

Each support receives all32 even-parity sign patterns with amplitude4. Every resulting point has norm squared96. Distinct supports have inner product at most3*16=48. Distinct signs on the same support have Hamming distance at least2 and squared distance at least128. Every dense coordinate used by this layer has absolute value at most2, hence dense/RS inner products are at most6*4*2=48. The layer contributes32*r*16³ points.

## Two-support root layer

For every unordered pair of coordinates within any one RSgroup, add all four sign choices sqrt(48)*(+/-e_i+/-e_j). In95D also use every pair consisting of the excluded coordinate88 and one eligible coordinate. Each root has squared norm96, and distinct roots in this D_n root system have inner product at most48.

Every RS support meets each root support at most once: within-group pairs use one group, and index88 is absent from RS supports. Thus absolute RS/root inner products are at most4sqrt(48)<48.

For93/94/96D, root/dense inner products are at most4sqrt(48)<48. For95D, pairs involving the special coordinate are bounded by(2sqrt(2)+2)sqrt(48)<5sqrt(48)<48. All inequalities are verified exactly after positive squaring (8<9 and25*48<48² suffice for the exceptional case).

The root counts are2r(r-1)+5*2*16*15, with an extra4*94 in95D. There are no coordinate poles in these totals.

## Counts and disjointness

Dense points have at least22 nonzero bodycells (24 in the full96D coordinate-alphabet family). RS points have exactly six nonzero ambient coordinates; roots exactly two. Thus the three layers are disjoint. The dense injection is proved in the included packages; support/sign and root-support injection are explicit here.

|Dimension|Dense|RS16|Roots|Total|
|---|---:|---:|---:|---:|
|93|42,949,672,960|1,703,936|2,712|**42,951,379,608**|
|94|64,424,509,440|1,835,008|2,764|**64,426,347,212**|
|95|107,374,182,400|1,835,008|3,140|**107,376,020,548**|
|96|206,158,430,208|2,097,152|2,880|**206,160,530,240**|

## Replay and scope

Run `python base/verify.py`, `python e7/verify.py`, and `python verify_augmentation.py`. The latter checks the actual listed head-coordinate bounds, GF16 arithmetic, all support counts and code rank conditions, group/root disjointness and exact geometric inequalities. Source parent minimum24 is an inherited audited premise, not silently inferred from these small scripts.

No publication or comprehensive worldwide-priority claim is made by this package.92D remains the separately verified two-basis construction because its sparse22 augmentation is not compatible by the arguments used here.
