> Archived independent pre-integration audit. Development paths below refer to the research workspace. All seven claims are now pinned in data/sections.json; packaging gaps listed below have been resolved by the production verifier.

# Independent audit of the seven Z8 row-section claims

## Scope and provenance

The five zero-sum candidates are in `../z8_sections_20261110/verified_sections.json`; the two equal-sum candidates are in `verified_equal.json`. They use the published Duursma–Greferath–Litsyn–Schmidt binary parent (2001), not a newly discovered parent code. This audit assesses correctness and improvement over the repository, not worldwide priority or optimality.

Current remote `RESULTS.md` was fetched on 2026-11-10, Git blob `1912c0d9aef1fafb6f9db2ee8d9a93cc7f6a65ab`. Every claimed baseline matches it. Its 70D entry remains 1,249,778,250 through the original cross-section package; a pending change is not reflected in this fetched file. Source URLs and blobs are in `source_metadata.json`.

## Fresh parent verification

Both full source proof implementations were run afresh on the candidate's actual data. Both PASS, and their Z8, support, and dimension fields agree exactly. The local canonical verifier is byte-identical to the freshly fetched repository version. The second verifier was fetched separately and imports none of the first verifier's proof code. Its input root was supplied by the runner without changing its mathematics.

Each proof checks all 16,777,216 Z4 words, 759 octads against 5,489 low-Lee-weight errors each, and 338,559,231 distinct support pairs. The proofs validate the free rank-12 Z8 generator, Gray-map distance, residue weights and octad-coset exclusion. The union with the first-bit-flipped coset has binary distance at least 24: cross-coset rows have opposite parity, hence contribute at least one differing bit in each of 24 rows.

The support replay validates unique weight-24 supports, intersection at most 12, disjoint odd trace obstructions against every GF(4) transversal, the GF(4) [24,9,12] code, and binary [24,14,6] signs. The canonical implementation's binary32 support products are exact sums of 0/1 entries; the independent proof checks the same pairs using packed integers. Neither proof enumerates the 2^37 binary words.

## Actual section transfer

Write S_t for the sum of the four physical coordinates in row t.

For dimensions 88–92, impose S_t=0 for the listed rows T. These equations have a displayed identity minor on one cell from each row, so rank is |T| and the ambient section dimension is 96-|T|.

For dimensions 93–94, impose S_t=S_(T[0]) for every other listed row. An identity minor on the nonreference rows proves rank |T|-1 and dimension 97-|T|. The two sets are [0,1,4,8] and [0,1,8], respectively. Their Z8 projections are surjective, since they are subsets of the first twelve columns of the triangular unit block. The independent checker also enumerates their complete Z8 projection images and checks both flipped and unflipped row-sum multiplicities directly.

This is a linear section followed by an orthogonal identification with R^d, not coordinate deletion or a distance-contracting projection. Inner products, norms and distinctness are unchanged. Zero-section configurations span the section because their A3 roots span each selected row's sum-zero subspace and their poles span all outside cells. Equal-section configurations also contain the all-+1 dense word, adding the one missing common-sum direction. Their span is therefore exactly the claimed dimension.

Unflipped dense words obey zero-sum conditions precisely when selected Z8 symbols are nonzero modulo four. Uniform homomorphism fibers give the exact count 8^12 times the fraction of valid projected Z4 images. In 88D the image has 32,768 elements and 3,281 valid images, so the non-surjective projection is explicitly accounted for. Uniform independent-row counts cannot be used blindly there.

The 26,022 sparse supports retain the original increasing physical-coordinate sign assignment. No permutation of the sign code is assumed. `check_counts.py` independently expands all binary sign words by matrix multiplication, evaluates row sums with signed-coordinate prefix sums, and aggregates only identical 24-row occupancy profiles. It imports no source proof functions and does not use the search histogram, zeta transform, or source packed bit-intersection count.

For equal sections, the transversal signs occupy the row positions themselves. Direct enumeration of all 16,384 sign words gives 2,048 valid patterns in 93D and 4,096 in 94D. Multiplication by 4^9 transversal supports gives the claimed 536,870,912 and 1,073,741,824. Zero sections have no transversal points.

The six-support layer is placed on 6q distinct outside coordinates, in increasing order. This is an explicit injective coordinate relabeling of its RS construction. There are enough outside coordinates in every case: (dimension,q,available)=(88,9,64),(89,11,68),(90,11,72),(91,11,76),(92,13,80),(93,13,80),(94,13,84). The pure layer has 12|T| A3 roots and 2(96-4|T|) poles, hence 192+4|T| points.

## Exhaustive pair-category argument

All points have squared norm 96, and the required product bound is 48.

- Dense/dense: binary distance >=24 gives product <=96-2*24=48.
- Dense/weight-24: <=24*2=48. Dense/weight-six: <=6*4=24.
- Two weight-24 points on one support: sign distance >=6 gives <=4*(24-12)=48. On distinct supports: intersection <=12 gives <=4*12=48. This includes additional/transversal intersections by the trace obstruction.
- Two weight-six points on different supports: RS distance >=3 gives intersection <=3, hence <=16*3=48. Same support: odd-parity signs differ in at least two places, giving <=16*(6-4)=32.
- Weight-24/weight-six: <=6*2*4=48, independent of the relocated support.
- A3 root/dense: <=2*sqrt(48); A3 root/weight-24: <=4*sqrt(48), both <48. Root/weight-six is zero because the latter avoids all selected rows.
- Outside pole/nonpure: <=4*sqrt(96)<48. Poles and roots have disjoint support. Different-row roots are orthogonal; same-row A3 roots have product <=48; distinct poles have product <=0.

Layers with support sizes 96,24,6,2,1 are disjoint. The two weight-24 families are disjoint because transversals occupy all 24 rows, whereas additional supports occupy at most 14. Generator freeness, distinct supports, and distinct sign words establish uniqueness within each family. Isometric section identification preserves it.

## Remaining packaging gaps

No mathematical defect was found in the section mechanism or inherited parent assumptions. The candidate README documents only the five zero-sum gains; the two equal-sum candidates need their own explicit conditions, counts and proof in the final package. The original candidate verifiers copy unchanged terms from baseline JSON rather than checking all seven full constructions end-to-end; this audit separately checks these terms and geometry. The audit supplements do not yet constitute a registered, manifest-pinned production verifier with new candidate-specific corruption controls. No remote files were changed or published. No optimality or worldwide novelty claim is supported by this audit.

## Final count outcome

All seven independent count checks PASS. Together with both parent replays and the section/layer checks, all seven claimed gains are substantiated against the freshly fetched repository baseline.

**Historical comparison.** The audited totals and gains are retained as recorded. The added current registry column follows `nonlinear/data/verified.json`.

|Dimension|Historical audited total|Gain over repository baseline at this audit|Current registry|
|---|---:|---:|---:|
|88|6,922,487,106|1,235,522|6,967,552,256|
|89|9,229,785,598|431,462|9,270,596,604|
|90|12,307,197,352|527,024|12,341,351,672|
|91|16,410,365,876|766,760|16,435,507,572|
|92|21,881,382,480|892,192|21,896,439,920|
|93|31,041,848,144|889,504|31,056,799,856|
|94|47,698,308,180|800,104|47,703,347,340|
