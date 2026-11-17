# Sparse22 augmentation of the new two-row root-head configurations,92–95D

This package extends the source-derived dense construction proved in DENSE_PROOF_92_94.md, without changing that frozen package. Its only inherited distance premise is the previously audited Z8 parent Gray minimum24. All new sparse code distances are checked here by complete exact finite enumeration.

## Dense95D extension

Retain the same22 Grayrows(bodydimension88) and use seven head coordinates. Coset0 uses disjoint row-axis lists A=(0,1,2),B=(3,4,5,6), retaining all12 support edges and48 signed heads. Coset1 uses A=(0,3,4),B=(1,2,5,6), retaining edges01,02,35,36,45,46 and24 signed heads. Each head is2*(+/-e_i+/-e_j), norm squared8. The same-coset Gram preservation and cross-coset body88 plus head8 separation proof is unchanged. Thus72*8^10=77,309,411,328 dense points lie in95D. `verify_dense.py` checks all maps, head pairs, uniformity and counts in all four dimensions.

## Explicit outer support code

The archived source supplies a GF4 generator of shape9x24, with GF4 represented by F2[x]/(x²+x+1). Shorten coordinate0 to zero and delete coordinates0,1. Row elimination constructs a8x22 generator. Complete enumeration proves65536 distinct words and minimum nonzero weight11. Each word chooses one cell in each of the22 retained four-cell bodyrows. Thus these65536 distinct22-element supports intersect in at most22-11=11 cells.

## Explicit sign code

Shorten the source binary14x24 sign generator at coordinates0,1 and delete both. Row elimination produces a12x22 generator. Complete enumeration proves4096 distinct words with minimum nonzero weight6. Use its signs on the22 cells, ordered by bodyrow. Every support carries all4096 signs.

Put amplitude a=sqrt(48/11) on these cells and zero elsewhere, including all head coordinates. Every point has squared norm22*a²=96. The sparse layer has65536*4096=268,435,456 distinct points.

- Different supports: inner product≤11*a²=48.
- Same support: binary distance≥6 gives squared distance≥4*a²*6>96.
- Against any new dense point: its retained body entries are+/-1, so absolute inner product≤22*a=sqrt(2112)<48. The head contributes0.

Hence this sparse layer is compatible with every new dense92–95D configuration. Sparse and dense points are disjoint: sparse bodies have exactly22 nonzero cells, dense bodies all88. Within the sparse layer the support and signs recover the construction data.

## Exact totals

92D:25,769,803,776+268,435,456 = **26,038,239,232**.
93D:38,654,705,664+268,435,456 = **38,923,141,120**.
94D:55,834,574,848+268,435,456 = **56,103,010,304**.
95D:77,309,411,328+268,435,456 = **77,577,846,784**.

No prior dense/sparse layers are unioned, and no compatibility with them is assumed. No coordinate poles or other optional additions are included in these totals. No publication is performed.

## Replay

Run `python verify_dense.py` and `python verify_sparse.py`. Both produce exact receipts. Sparse generator elimination, cardinalities, full nonzero-weight histograms and all geometric comparisons use finite-field/integer/rational arithmetic. The square-root inequality is verified after squaring positive quantities. The copied source_supports.npz is pinned by SHA256 in sparse_certificate.json. Only its two generator matrices are used; its other legacy arrays are not mathematical inputs to the new sparse proof.
