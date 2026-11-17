# Independent audit: direct91D/92D latitude layers

PASS under the already established Pless [48,24,15] minimum-distance premise and the previously audited21,943,640,934-point90D equatorial union. These premises are inherited, not independently exhaustively re-enumerated in this audit. The generator is checked identical to the prior hash-pinned sparse15 certificate.

The independent checker verifies all six31680-word weight14 body arrays: ternary alphabet, uniqueness, systematic membership after indicated three-symbol prefixes, weight and exact counts. It verifies the2609-word shortened phase14 code by a separately coded packed-symbol XOR/popcount traversal of all3,402,136 pairs, minimum5. It also verifies four6930-word weight13 arrays and all393,828 pairs of the888-word phase13 code. No completeness of the outer shells is needed for these lower bounds.

Signed A2 roots have squared norm3, distinct same-sign roots distance9, opposite-sign roots distance at least3, and root/triangle absolute inner product at most3/2. These are checked with exact integer dot products in coordinates(a/2,b*sqrt3/2). Disjoint positive/negative alphabets and distinct phases make the construction injective.

For a full parent prefix p of weightw, the body has norm squared45-3w and its head has norm squared3w. All points therefore have norm squared45. For different outer codewords with prefixes p,q, body Hamming distance is at least15-dH(p,q), giving body squared distance at least45-3dH(p,q). The signed grid heads sqrt3(f(p1),f(p2)), f(0)=0,f(1)=1,f(2)=-1, have squared distance at least3dH(p,q), checked on every ordered pair of the nine prefixes. This verifies all cross-layer pairs, including equator prefix00, four axis layers and four corner layers. For the same outer word, phase distance5 gives body distance at least45 directly. In91D use only opposite scalar heads.

Against dense triangle equator points, a weight14 body has absolute inner product at most21, and weight13 at most19.5, both below22.5. Against the old sparse15 equator, the common Pless parent plus the prefix-distance argument above gives total distance at least45. No additional unverified old sparse layers are unioned.

The two optional91D poles have norm squared45 and zero body. Their cross inner product is bounded bysqrt135<22.5. The six optional92D orthogonal hexagon points also have norm squared45; against corner heads the bound is sqrt270<22.5, and their own maximal inner product is22.5. All are disjoint from nonzero-body layers. New nonzero heads distinguish latitudes from the equator, so no point or antipode is accidentally counted twice. The count is not doubled again for antipodal completion.

Counts:
- Each axis layer31680*2609=82,653,120.
-91D:21,943,640,934+2*82,653,120+2 = **22,108,947,176**.
-92D square version:21,943,640,934+4*82,653,120+6 = **22,274,253,420**.
- Four corner layers add4*6930*888=24,615,360.
-92D full grid version: **22,298,868,780**.

This audit establishes the stated construction from inherited validated premises and exact finite data. It does not assert global novelty, a new exhaustive parent minimum proof or permission to publish. Run python -B check.py; NumPy is the only nonstandard dependency.
