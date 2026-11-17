# Two-coset root-head compression in dimensions92,93,94

## Scope and inherited source

The source is the previously verified free rank12 Z8 code C of length24 in `dim88-96-z8-golay`, with homogeneous minimum24 and its displayed generator. Its Gray image has binary minimum24. This minimum premise is inherited from that frozen source proof; the present small verifier does not repeat the large Z4 enumeration. The public parent construction is attributed to Duursma,Greferath,Litsyn,Schmidt (2001), DOI10.1109/18.923742. No new parent code is claimed.

The copied generator JSON SHA256 is903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73. Its first12 columns are upper triangular with odd diagonal; hence C is free and its first two ring coordinates are uniform. The upper-left2x2 block is[[7,1],[0,7]], invertible modulo8, and all later generator rows vanish in these two positions. Thus every fixed pair(a,b) occurs in exactly8^10 distinct codewords.

For ring symbol a, the four Gray bits indexed by(x,y) in F2² are

    (a>>2 mod2) XOR ((a mod2)*x) XOR ((a>>1 mod2)*y).

Let s(a) be their sign vector. Each has norm squared4. The eight vectors form a four-dimensional crosspolytope: a mod4 selects its axis, while floor(a/4) selects its sign. All Gray rows have even parity. Flipping the first physical bit of every row produces a second coset of odd row parity and preserves all within-coset distances and inner products.

## Construction

Delete the first two four-bit rows, retaining22 sign rows in R88. Replace the deleted rows by the head

    2*((-1)^floor(a/4)*e_A[a mod4] + (-1)^floor(b/4)*e_B[b mod4]).

Only listed indices in A and B are used. In each coset the two axis lists A,B are disjoint. Consequently the head's Gram matrix equals the original two-row Gram matrix exactly for every selected pair of prefixes in that coset. Heads have norm squared8; bodies have norm squared88; all resulting points have norm squared96.

Use the following head maps and physical support edges. Each edge permits all four independent sign choices.

Dimension92 (head R4):
- Coset0: A=(0,1),B=(2,3), all four edges between these lists,16 heads.
- Coset1: A=(0,2),B=(1,3), only edges01 and23,8 heads.

Dimension93 (head R5):
- Coset0: A=(0,1),B=(2,3,4), all six edges,24 heads.
- Coset1: A=(0,2),B=(1,3,4), only edges01,23,24,12 heads.

Dimension94 (head R6):
- Coset0: A=(0,1,2),B=(3,4,5), all nine edges,36 heads.
- Coset1: A=(0,3),B=(1,2,4,5), only edges01,02,34,35,16 heads.

`certificate.json` lists every ringprefix, coset index and head vector explicitly. A ring axis index is its position in A or B, not necessarily the physical head axis number.

## Separation

Within a coset, replacing the prefix preserves its exact Gram matrix and norm, while the body is unchanged. Thus every squared distance equals its original full96-bit sign-code distance, at least4*24=96.

Across cosets, each retained row compares an even-parity binary word with an odd-parity word. They differ in at least one bit, giving body squared distance at least4*22=88. The two head support-edge lists are disjoint. Heads in the two cosets are therefore distinct vectors of the form2*(+/-e_i+/-e_j), and any two distinct such roots have squared distance at least8. Total cross distance is at least96. This argument requires no additional unknown cross-distance property of C.

Since all vectors have squared norm96 and all distinct squared distances are at least96, normalization gives a kissing configuration. For fixed prefix, deleting it is injective on full codewords. Different prefixes within a coset have different heads; the cosets also have different heads. Hence every counted point is distinct.

## Counts

-92D:24*8^10 = **25,769,803,776**.
-93D:36*8^10 = **38,654,705,664**.
-94D:52*8^10 = **55,834,574,848**.

These are new configurations, NOT unions with previously published dense/sparse layers. No compatibility with those other layers is asserted. No95D improvement is claimed. This package performs no publication. Independent audit is requested separately.

## Replay

Run `python verify.py`. It checks generator uniformity, all64 scalar Gray comparisons implicitly through Gram/parity arrays, every selected same/cross-coset head pair, head uniqueness, ambient dimensions and exact counts. The source homogeneous minimum24 remains the explicit inherited premise.
