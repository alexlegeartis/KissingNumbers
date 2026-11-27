# 60D: affine-cube compression into a D4 tag constellation

**Exact construction: 82,821,160 kissing points in 60 dimensions.** This is **13,116,536** above the opening project bound of 69,704,624 (opening main ab408918cfddac0bc96d6c3585c608e5bb139d41), or about **18.82%**. [Kamenetsky, Table 3](https://arxiv.org/html/2611.01906v1#S6) gives the separately attributed published external comparison 57,573,374. Neither is an upper bound or a worldwide priority claim. This is an entirely new construction; do not add its layers to the old 60D stock.

## The folding principle

Use the [published 80-cubic-coset parent](https://github.com/alexlegeartis/KissingNumbers/tree/ab408918cfddac0bc96d6c3585c608e5bb139d41/verifications/improved/dim62-67-rm2-cubic-cosets). Its length-64 sign words have Hamming distance at least 16, with each homogeneous cubic translate containing every quadratic/affine Boolean polynomial. The minimum-distance verification and literal 80 cubic representatives are inherited; stored witness text SHA256 5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f. Use coordinate labels in F2^6, numeric bit order.

Mark the eight sites in U=span(3,20,41), ordered as the four matched F sites followed by their four P partners. Match each f_i of F=(0,20,41,61) with p_i of P=(3,23,42,62), in the stated order. To these marks assign signs along the four tail axes: π(f_i)=−e_i and π(p_i)=+e_i. Keep the other 56 coordinates unchanged. Substituting a four-dimensional tag for the eight marked coordinates gives 60 dimensions.

For a marked binary label S (negative signs at sites in S), use these tags:
* for every two-site label {a,b} not a matched pair, tag h=−2(π(a)+π(b)); there are 24, the scaled D4 roots, squared tag norm 8;
* for every singleton {a}, tag h=−π(a); there are eight, squared tag norm 1;
* for the empty label, tag h=+e_1, also squared tag norm 1 (using the first P representative).

The binary evaluation of RM(2,6) on U is the entire even subspace of F2^8: the 3-flat restriction has all seven degree≤2 coefficients, and every such word has even weight. Of the 80 actual cubics, 64 have even parity on U and 16 odd. Thus every even label has 64×2^15 parent words and every odd label 16×2^15, giving **56,623,104 dense points**. Their squared norms are 64 or 57. For distinct parent words with labels S,T the ordinary dot is at most 24+2|S△T|, directly from parent Hamming distance 16 and subtracting the literal eight signs. The checker verifies
  24+2|S△T|+h_S·h_T ≤ √((56+||h_S||²)(56+||h_T||²))/2
for every pair of labels, also S=T. No pointwise conditional Hamming improvement is assumed. The combinatorial reason the 24 D4 heads work: disjoint two-site labels give nonpositive head dot (any repeated physical axis comes with opposite sign), one common site gives head dot ≤4, and the same label gives dot 8. The singleton and empty tag geometry uses the lower norm, not a preliminary floating-point fit.


## Rephase the whole previously verified sparse family through the same coordinates

Take the **13,013** distinct weight-16 supports of the prior 60D signed construction avoiding F. These comprise 12,005 parent weight-16 words plus 672 masks in the published first-stage support guards and 336 in the published guard transfer. The replay reconstructs the first 12,005 by restricting all 80 cubic cosets to vanish on F and enumerating the 2^18 solutions per coset. It reads the other 1,008 from the two pinned source fixtures:
* [extra60.json.gz.b64, first-stage source](https://github.com/alexlegeartis/KissingNumbers/blob/ab408918cfddac0bc96d6c3585c608e5bb139d41/verifications/improved/dim60-61-rm80-puncture/extra60.json.gz.b64), SHA256 8420af6106dcc6a6a2bda865593813e14e36f76dc547bac2accfcce320d05ac2.
* [transfer60.txt.gz.b64, guard-transfer source](https://github.com/alexlegeartis/KissingNumbers/blob/ab408918cfddac0bc96d6c3585c608e5bb139d41/verifications/improved/dim60-61-rm80-guard-transfer/transfer60.txt.gz.b64), SHA256 7fb3b5ed5b309bc87f05df08985e45e0f9d14c240835f4404b6912b748dcf65e.

On each old support the sign bank has overall odd parity. The 3,112 literal pair guards set opposite parities on their shared eight sites. Independently replay all pairs of support masks: intersections above eight are exactly those guarded pairs, and their sizes are nine or ten. Thus any such pair must have an ordinary sign disagreement in the common physical chart and its raw dot decreases by at least two, returning to the threshold eight. The joint affine checks on each support are consistent. Every nonzero homogeneous row other than the full support has weight eight; hence when independent rows modulo the all-one row are restricted to sites, the induced binary column classes have equal sizes (Walsh inversion). Split each class into distinct extra binary labels to complete all columns (1,a), a in F₂⁴. Use zero syndrome on the new independent rows, preserving the old prescribed signs on the original rows. This gives a literal affine [16,11,4] sign bank with 2,048 words per support; the replay completes these charts in sorted physical coordinate order. It does not trust a prior solver or assume full pair compatibility without checking it.

Embed each surviving sparse vector through the same isometry: P is sent to its four distinct +e_i and the 56 ordinary coordinates are kept; supports never use F. Hence no prior sparse–sparse distance changes. Binned by c=|S∩P|, the exact support counts are:

| Number c of included P marks | Mother supports | Additional guarded supports | Total |
|---:|---:|---:|---:|
| 0 | 5,261 | 54 | 5,315 |
| 1 | 1,924 | 811 | 2,735 |
| 2 | 4,456 | 63 | 4,519 |
| 3 | 244 | 80 | 324 |
| 4 | 120 | 0 | 120 |

Keep the entire odd sign bank for c≠3. For c=3 fix the three marked signs positive; the completed dual has minimum weight eight, so its affine primal coset projects onto all patterns on those three positions and exactly 256 signs remain. The shorter sign banks retain all old sparse/sparse checks. The new sparse population is **26,070,016**.

## Dense–sparse interface

The exact replay verifies that every one of the 13,013 sparse masks is orthogonal modulo two to every dense generator: the 22 quadratic/affine rows and the 80 cubic truth tables. Thus a sparse sign pattern of odd parity must disagree with the original parent's signs somewhere on the support.

For a tag h and a subset J of the available marked sites P, let c=|J| and a_j=+e_j·h. If c≠3 the unguarded mixed upper bound is (16−c)+∑_{j∈J}|a_j|. Whenever this exceeds the threshold for the actual dense norm (64 or 57) and sparse norm 16, every marked coefficient present is nonzero and its sign is the original dense sign there. The ordinary coefficients already have the original dense signs and magnitude one. The forced parity disagreement thus decreases this upper bound by at least two. For c=0 the condition on marked coefficients is vacuous. For c=3, since the three marked sparse signs were fixed positive, use the stronger bound (16−c)+∑ a_j. If it is still above threshold then all three marked original dense signs are positive (verify this for every literal label), so the forced disagreement occurs on an ordinary site, again subtracting two. These finite checks include all 33×16 tag/subset pairs and use exact squared comparison; no floating-point criterion.

Add the 7,995 affine 2-flat supports avoiding F, with all sixteen magnitude-one signings, mapped by the same sparse isometry. These **127,920 points** have squared norm four and distinct supports overlap at most two. Their same-support distinct signs have raw dot at most two; mixed with a sparse support the raw dot is at most four. Against dense vectors the absolute-value bound (4−c)+∑ |a_j| is at most six and verified. Add the **120 coordinate poles**, since the largest normalized component of any other family is at most one half. These cross-angle checks also eliminate coincident normalized points. The original project bound is not added.

| Layer | Population |
|---|---:|
| 24 D4 double-mark labels | 50,331,648 |
| 8 singleton labels | 4,194,304 |
| even empty label | 2,097,152 |
| complete and shortened odd sparse banks | 26,070,016 |
| all-signed affine four-sets | 127,920 |
| poles | 120 |
| **Total** | **82,821,160** |

## Exact replay

Run the companion as Node.js with the local repository root, for example:
    node verifications/improved/dim60-d4-affine-cube-fold/verify.mjs .
It requires only Node's standard library and the three pinned existing files, no native compiler, external solver or tolerance. The Hamming-16 distance of the 80-coset input remains the published and hash-pinned premise; it can be independently reproduced with its own package. The new checker independently regenerates the sparse mother stock, the added guarded stock and all guard phase completions, all 84,662,578 support pairs, the elementary angle inequalities and every count. The numeric hash pins compare exact compressed-text bytes before parsing. The same paper's sparse blocks are not counted on top of the new result. An intermediate proof that used only mother supports obtained 80,747,560 and is strictly weaker; it is not needed for replay.
