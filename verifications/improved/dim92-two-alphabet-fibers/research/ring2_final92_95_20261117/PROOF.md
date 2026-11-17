# Final two-coset construction with sparse22, RS17 and poles,92–95D

The dense proof and its inherited Z8 minimum24 premise are given in DENSE_PROOF_92_94.md and SPARSE22_PROOF.md (the latter includes the95D dense extension). The sparse22 layer is proved there and exactly verified by verify_sparse.py. This document adds two mutually compatible layers to those new constructions; no legacy layer is silently retained.

## Reed–Solomon six-support layer

Let d=dimension-88, so4<=d<=7. OverF17 evaluate every polynomial p(x)=a0+a1*x+a2*x²+a3*x³ at x=0,1,...,5, restricting a0 to0,...,d-1 and allowing all17³ choices of a1,a2,a3. Use p(0) as a head coordinate (physicalindex88+p(0)); put p(j),j=1,...,5, into the jth disjoint17-cell bodygroup (physicalindex17*(j-1)+p(j)). All supports therefore have one head cell and five body cells, within the88 available body coordinates.

Distinct cubics agree at at most three evaluation positions, hence distinct supports intersect in at most three cells. The verifier also checks all15 four-column minors of the evaluation generator are nonzero modulo17, and enumerates every support for cardinality/injectivity. Every support receives all32 even-parity sign patterns, with amplitude4. These vectors have squared norm96. Distinct supports have inner product<=3*16=48. Same-support signs have Hamming distance>=2, giving squared distance>=4*16*2=128.

Against dense points, the five body entries contribute at most20 and the one head entry at most8, totaling28<48. Against sparse22 points, the head entry contributes0 and only five body entries can meet the support, so absolute inner product<=20*sqrt(48/11)<48. Thus the RS layer is compatible with both preceding layers.

Its size is32*d*17³=d*157216.

## Coordinate poles

Add +/-sqrt(96)*e_i for every ambient coordinate,2*(88+d) points. Their own pairwise inner products are0 or-96. Their largest possible absolute cross inner product with any RS point is4*sqrt(96)<48; with sparse22 it is sqrt(96*48/11)<48; with dense points it is at most2*sqrt(96)<48. All comparisons are checked exactly by squaring positive quantities.

All four layers are disjoint: dense vectors have88 nonzero body cells and two nonzero head cells; sparse22 vectors have22 nonzero body cells and zero head cells; RS vectors have five body cells and one head cell; poles have a single nonzero coordinate.

## Final counts

| Dimension | Dense | Sparse22 | RS17 | Poles | Total |
|---|---:|---:|---:|---:|---:|
|92|25,769,803,776|268,435,456|628,864|184|**26,038,868,280**|
|93|38,654,705,664|268,435,456|786,080|186|**38,923,927,386**|
|94|55,834,574,848|268,435,456|943,296|188|**56,103,953,788**|
|95|77,309,411,328|268,435,456|1,100,512|190|**77,578,947,486**|

Run verify_dense.py, verify_sparse.py and verify_rs17.py. All new finite sources, ranks, distances, heads and geometric comparisons are explicit. The Z8 parent minimum24 is the only inherited distance premise. No publication is performed. Independent audit is requested separately.
