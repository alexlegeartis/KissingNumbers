# Two complete trine fibres in 36 dimensions

The preceding [Coxeter–Todd certificate](../dim36-ct-triangle-refill/) gives 502,902 points. Over the two literal new triples of Coxeter–Todd phase cells, use all 27 independent trine choices per triple. The two sets of 27 new points have the same 27 blockers in the old root stock. Remove those blockers and retain both sets, giving **502,929**.

The exact integer verifier checks the literal phase data and source hashes, every pair of new points, all 5,356,692 products with the old root stock, and 77,490 checks for the other inherited base families. Run `python3 verify.py` with NumPy and the adjacent preceding package. The replaced stock, unrelated inherited families, and separating inner-product threshold are not treated as numerical heuristics.

We compare with the already verified repo bound, and retain the table's named older [Yang et al. comparison](https://arxiv.org/abs/2609.35051v1). This is a lower bound, not an optimality claim.
