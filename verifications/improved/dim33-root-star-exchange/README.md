# 370,556 kissing points in dimension 33

**Certified construction, not a reclaimed record.** The witness improves [the previous repository configuration](../dim33-joint-signed-refill/) from 370,528 to 370,556. [P. Carder reports 370,568](https://doi.org/10.5281/zenodo.24577265); his coordinate artifact was not independently replayed here.

| Layer | Points |
| --- | ---: |
| Dense binary signs | 131,072 |
| Complete signed octads | 235,136 |
| Individual signed octads | 2,216 |
| Signed stars | 132 |
| Signed two-coordinate roots | 2,000 |
| **Total** | **370,556** |

All coordinates are integral at common squared norm 32; the separation threshold is inner product 16.

## Rebuilding the two-coordinate roots

Every one of the 2,112 signed roots with two entries of magnitude 4 is compatible with the dense code and with all octads. Among the previous 172 stars (one head of magnitude 4 and four leaves of magnitude 2), build a bipartite graph of star-root conflicts. Each star excludes exactly four roots, one per signed leaf.

A matching of size 152 is a literal obstruction to keeping more than `172+2112−152=2132` points from the two families. Its alternating forest provides an attaining choice: 132 of the old stars exclude only 112 roots. Retain every other one of the 2,112 roots. The dense and octad layers do not move.

The compact certificate is [witness.npz](witness.npz). The independent checker verifies all norms and each pair category exactly. The optional expanded replay also checks all coordinate products directly.

```sh
OPENBLAS_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 python verify.py --all-pairs
```

To rediscover the matching construction, run `build_exact.py` on the previous, SHA-256-pinned [370,528 witness](../dim33-joint-signed-refill/witness.npz); NumPy and SciPy are required for that step.

## Why retaining all the old individual octads cannot reclaim 33D with conventional stars

This is an upper bound only for the subclass in which all complete banks, all dense signs and all 2,216 individual octads remain fixed, and the only other point types are the signed roots above and conventional stars. It is not an upper bound on the 33D kissing number.

Group a compatible collection of stars by signed head. In one group the residual leaf vectors have pairwise nonpositive products, since their head already contributes 16 to the dot product. If the group uses one sign in each of `s` leaf coordinates and both signs in each of `d` others, it has at most `s+2d` stars, the number of signed leaf coordinates in use. To see this, decompose the positive semidefinite leaf Gram matrix into its connected negative-edge blocks. Each singular block has nullity one and a strictly positive null vector (Perron–Frobenius). Such a positive dependence cannot involve a coordinate used with only one sign. The spans of the singular blocks are orthogonal, so there are at most `d` singular blocks, and the entire Gram rank is at most `s+d`.

For each signed head `h`, let `R_h` be the set of roots it excludes and let `d_h` be `|R_h| − (# stars at h)`, which is nonnegative. A root can be excluded from at most its two signed endpoints. Thus, if every compatible root is kept, the gain over using only the full root layer is precisely the number of roots excluded from both endpoints, less `sum_h d_h`.

Exact enumeration of the four-leaf supports and their signs finds 1,486 eligible stars and only 34 potential roots that could be excluded from both endpoints. Outside the four signed heads on coordinates 20 and 30, there are 28. If any such exceptional signed head is used, its local deficit is at least 3, verified by exhaustively checking its at most four available stars; yet it can add at most two reciprocal exclusions. Therefore the entire gain over the full-root baseline of 370,536 is at most 28. The subclass ceiling is **370,564**. To reach 370,569, one must move some of the old individual octads or use a different layer mechanism.

Run `OPENBLAS_NUM_THREADS=1 python reciprocity.py` with NumPy 2 to verify the finite inputs to this bound. The theorem uses all 1,486 eligible stars, not just the 172 stars used in the prior witness.
