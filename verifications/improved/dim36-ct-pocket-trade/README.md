# Coxeter–Todd pocket trades in 36 dimensions

The [preceding two-trine construction](../dim36-ct-double-trine/) gives 502,929 points. Here 19 local trades use 38 root phase cells across three 12-dimensional Coxeter–Todd blocks. In each pocket remove nine complete old trines (27 points). Insert the two compatible (3\times3\times3) root phase cells (54 points), except that two new cells are partially filled with 24 and 21 points. The removed trines are the shared blockers of those new cells. The 36 complete cells and two partial cells give 1,017 new vectors; 513 old vectors are removed. Refill 115 points into the previously empty block of the vacated two-block fibres. These use 47 different vacated fibres and preserve every accepted root phase. The resulting configuration has **503,548** vectors.

No old root triple is removed. Every new point has squared norm 12. The verifier reconstructs the preceding stock, checks all new-to-old integer Gram products and recomputes the removal set (100,911,825 old-root and 408,260,412 old two-block products). It also checks 57,894,795 exact products involving the restored points, all new pairs and the exact count. The inherited old-to-old checks are supplied by the linked predecessor, whose input hashes are pinned here. No numerical tolerance or search step is involved.

Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` with NumPy. The historical benchmark is [Yang et al.](https://arxiv.org/abs/2609.35051v1), 484,760; this is a lower bound, not an optimality claim.

For an independent literal expansion and coordinate Gram replay, run `OPENBLAS_NUM_THREADS=1 python3 audit.py`.
