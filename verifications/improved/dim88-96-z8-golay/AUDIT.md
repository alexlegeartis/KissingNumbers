> Historical first-pass report: its 88D candidate is superseded by FINAL_88_AUDIT.md. report.json and audit_final.log contain the final replay; audit.log retains the earlier run. Other dimensions below remain unchanged.

# Independent adapted-sign audit

All five claims PASS. The audit is `audit.py`, with machine-readable receipt `report.json` and log `audit.log`.

|Dimension|Total|Gain over current published repository|
|---|---:|---:|
|88|6,937,411,712|14,924,606|
|89|9,243,338,748|13,553,150|
|90|12,316,642,276|9,444,924|
|91|16,416,373,024|6,007,148|
|92|21,883,990,128|2,607,648|

Fresh remote RESULTS.md, blob `5ee520eafa6c797f191b03a79bc9d948c541ac74`, already contains the preceding seven row-section improvements. These five candidate baselines match it, not the earlier superseded baseline. The same fresh file now lists 70D at 1,360,986,112.

## Independent reconstruction

The auditor uses scalar XOR to enumerate the binary sign code, verifies 16,384 unique words and minimum distance six, and recounts every stored adapted pattern using scalar bit counts. It imports no source proof code. It validates disjoint pattern groups of the prescribed sizes and all position bounds.

For every retained support in every dimension, it builds an actual bijection from its 24 occupied physical cells to the 24 sign positions. The deterministic convention is:

1. List selected rows with two occupied cells, in increasing physical row order, then selected rows with four occupied cells, likewise ordered.
2. Assign the stored pattern groups in that order, with cells inside each row ascending and sign positions in the supplied group order.
3. Match all remaining occupied physical cells in increasing order to all remaining sign positions in increasing order.

Every constructed mapping is checked to be bijective and to pull the physical selected-row equations back to exactly the pattern groups whose signs were independently counted. Selected rows of odd occupancy cannot sum to zero and are omitted. Selected empty rows impose no restriction.

The numbers of explicitly reconstructed support maps for dimensions 88 through 92 are 24,349, 24,479, 24,557, 24,747, 25,036. Their pattern histograms exactly match the candidate receipts. Their weighted sparse totals are respectively 56,445,824; 69,927,168; 85,583,084; 108,450,348; 139,804,032.

The projected Z4 images are independently enumerated, without importing a baseline count. Their sizes for 88 through 92 are 32,768; 16,384; 4,096; 1,024; 256, and their all-nonzero counts are 3,281; 2,187; 729; 243; 81. The 88D projection is correctly treated as non-surjective. Exact uniform fibers of the free rank-12 Z8 code give the claimed dense counts.

## Compatibility, dimension and injectivity

Different supports do not need a common sign permutation. They intersect in at most 12 physical cells, so any two signed weight-24 vectors on different supports have inner product at most 4*12=48. Within one support a coordinate bijection preserves Hamming distance at least six, hence the same bound 4*(24-2*6)=48. The new sparse code may differ from the old code on each support; this does not affect either argument.

The underlying support file is byte-identical to the file covered by both fresh full parent proofs in `../z8_sections_audit_20261110/`. Its SHA-256 is `d0def7d6587ca3aa1b39c9f3e7ccb844310285a7f1794dc8e9492f0fff9fc65d`. Those two proofs already rechecked all 338,559,231 support pairs and the parent Z8 distance. This follow-up does not claim another complete support-pair replay.

All points retain squared norm 96. Dense/sparse products are at most 24*2=48. Sparse/six-support products are at most 6*2*4=48 regardless of sign permutations. An A3 root has product at most 4*sqrt(48)<48 with a sparse point, and an outside pole at most 2*sqrt(96)<48. Other pair categories are unchanged from the preceding audit: six-support vectors remain outside selected rows, and the dense and pure layers use the exact same parent construction rules. No transversals survive the zero-row equations.

Each selected row contributes one independent zero-sum equation, as certified by the identity minor on one physical cell per row. Hence the section dimension is exactly 96-|T|. Roots span the three-dimensional zero-sum subspace in each selected row, and poles span all outside coordinates, establishing actual span of that dimension. Orthogonal identification with R^d preserves distances and distinctness. Different support sizes separate dense, sparse, six-support, root, and pole layers; unique supports and within-support bijections prevent duplicates.

## Gaps and limits

There is no substantive mathematical defect in the five proposed configurations. The source verifier counts abstract groups but does not instantiate their physical-coordinate assignment; the deterministic convention above supplies that missing explicit recipe and checks every support. The source verifier also does not itself establish each dimension's constraint rank or replay inter-support geometry; those must remain linked explicitly to the parent proof and the section argument. Its baseline is loaded from local files rather than checked against the remote repository; this audit performed that fresh comparison.

These remain research artifacts until integrated with pinned manifests, a full end-to-end verifier and candidate-specific rejection controls. Pattern search is heuristic, so neither optimality of each pattern nor optimality of each row selection is established. No worldwide-priority conclusion follows. No remote files were changed by this audit.
