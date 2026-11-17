# A new15-block A2-root sparse layer compatible with triangle-body codes

This is an unsubmitted research construction. It REPLACES prior sparse layers unless their cross compatibility is separately proved. No external edits or publication are authorized by this package.

Let the dense configuration have squared norm45. Suppose m of its orthogonal two-dimensional blocks are unit triangle symbols
 t0=(1,0), t1=(-1/2,sqrt(3)/2), t2=(-1/2,-sqrt(3)/2),
possibly with additional head coordinates, and possibly globally negated. Sparse vectors below have zero in all additional head coordinates.

Define three A2-root directions
 h0=t0−t2=(3/2,sqrt(3)/2),
 h1=t1−t0=(-3/2,sqrt(3)/2),
 h2=t2−t1=(0,−sqrt(3)).
They have squared norm3 and pairwise inner product−3/2. Every inner product between any ±hi and any triangle symbol tj belongs to{−3/2,0,3/2}.

Choose an outer ternary linear code C of lengthm and minimum Hamming distance15. Let S be any explicitly verified set of its weight15 words. Choose a ternary phase code K of length15 and minimum distance5. For each c∈S, list its15 nonzero coordinates in increasing order. For each k∈K, put h_(k_j) at the jth occupied block when c has symbol1, and −h_(k_j) when c has symbol2. Put zero at all other blocks. This defines |S|*|K| DISTINCT sparse vectors, all of squared norm45.

## Sparse–sparse distances

For different outer words c,c', every coordinate where their ternary symbols differ contributes squared Euclidean distance at least3, regardless of phases:
- zero versus a root contributes3;
- positive-root versus negative-root contributes either3 or12.
Thus different outer words give squared distance at least3*d_H(c,c')≥45.

For the same outer word and different phase words, a differing phase gives squared distance9, so total squared distance is at least9*d_H(k,k')≥45.

The norm45 sparse points therefore have mutual inner products at most22.5. If S is closed under negation, the layer is antipodal by construction, since negating c flips every root without changing k.

## Sparse–dense compatibility

Only15 blocks contribute. Each sparse–dense block inner product has absolute value at most3/2, so the total absolute inner product is at most15*(3/2)=22.5. Additional dense heads contribute zero. This proves compatibility with the ENTIRE dense configuration and its antipodes, independent of its code parent or prefix-head details, provided the m body blocks really are unscaled unit triangles and total norm is45.

Sparse vectors cannot coincide with dense vectors when m>15: sparse vectors have zero body blocks whereas every dense triangle body block has norm1. Thus union cardinalities add exactly.

## Explicit nonlinear phase code and verification

The phase code K is the7812-word ternary length15 code released by Antti Laaksonen and Patric R. J. Östergård with their2017 ICMCTA paper, New lower bounds on error-correcting ternary, quaternary and quinary codes. The authors' page is https://www.cs.helsinki.fi/u/ahslaaks/codes/ and the data archive is https://www.cs.helsinki.fi/u/ahslaaks/codes/codes2a.zip . The original member3_15_5_7812.txt is copied unchanged as phase7812.txt. Its SHA256 is d5ba3940239da9750271082790a2ac3f62aa52a14c4ee66727f9a66121f593e2.

All7812 words are distinct, ternary, and of length15. All30,509,766 unordered pairs were checked exactly; minimum distance is5. No linearity or phase-code antipodal closure is required by the theorem. The outer set is antipodal, and negating its sign word while retaining the same phase word negates the full sparse vector.

The explicit outer files and systematic Pless generator are included. All selected words are distinct weight15 code members and satisfy their specified shortening zeros. Completeness of the parent min15 shell is NOT needed for these lower counts. The parent code's minimum15 is retained as the established published/full-enumeration premise.

| Body blocks | Shortened coordinates, zero-based | Outer words | Sparse points |
|---:|---|---:|---:|
|45|0,1,2|130944|1022934528|
|44|0,1,2,3|87296|681956352|
|43|0,1,2,3,4|57536|449471232|
|42|0,1,2,3,4,5|37432|292418784|
|42|12,14,25,41,42,44|37468|292700016|
|41|0,1,2,3,4,5,6|24036|187769232|

A90-dimensional dense triangle layer of2*3^21=20,920,706,406 points plus the45-block sparse layer yields21,943,640,934 points. The resulting90-dimensional configuration embeds in91 and92 dimensions; respectively two or six norm45 points may be added in the unused orthogonal complement, giving21,943,640,936 and21,943,640,940. The six complement points are a regular hexagon; their pairwise inner products are at most22.5 and their cross inner products with the entire90-dimensional construction are zero.

This is a research result pending independent audit and owner direction. It must not be silently combined with old RS9/SQS sparse counts, whose cross compatibility is not established. No public repository edits were made.
