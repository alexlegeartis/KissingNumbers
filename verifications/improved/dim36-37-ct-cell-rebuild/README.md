# 36D and 37D: 502,767 and 506,709 by rebuilding Coxeter–Todd cells

This construction raises the [preceding 501,822-point bound](../dim36-37-orbit-exchange/)
by **945** and raises the preceding 505,638-point 37D bound by **1,071**. The 36D code is given by 5,098 outer cosets of a
product of three Coxeter–Todd lattices. The plain-text certificate is
`certificate.json`. The [external Yang et al. comparison](https://arxiv.org/abs/2609.35051v1)
is 484,760; this is a lower-bound construction, not an optimality
claim. The previous proof and manuscript are unchanged.

Use a copy of the Coxeter–Todd lattice in each of three orthogonal
12-dimensional blocks. The verification checks the 12-dimensional Gram
matrix and an index-729 similar sublattice, and enumerates all vectors
of squared norm at most eight. This gives exact smallest-norm distances
and point multiplicities in each coset of the sublattice. The final
vectors all have squared norm 12.

| Tuple of block minima | Outer cosets | Points per coset | Total |
|---|---:|---:|---:|
| (0,0,0), using a sublattice minimum in one block | 1 | 2,268 | 2,268 |
| Permutation of (0,4,8) | 756 | 243 | 183,708 |
| Permutation of (0,6,6) | 672 | 324 | 217,728 |
| (4,4,4) | 3,669 | 27 | 99,063 |
| **Total** | **5,098** | | **502,767** |

The first three rows retain the earlier full 403,704-point base. The
search used the ternary discriminant form to build a 729,000-cell
candidate palette and repeatedly reoptimized whole coordinate fibres.
The construction does not assume that this palette or the output is
optimal. Exact replays check the finished certificate directly,
independently of how it was found.

Two separate verifiers enumerate the short vectors, recompute coset
multiplicities, and check all 12,992,253 distinct outer-coset pairs.
Their minimum squared quotient distance is 12. Within an outer
coset, distinct points differ by a nonzero vector of the product
sublattice, also of squared norm at least 12. Thus every distinct
normalized pair has inner product at most 1/2.

## 37D lift

Use the entire 36D set as the equator. Select 219 outer cells with
permuted block minima (0,4,4), each containing nine norm-eight points.
At heights +2 and −2 these give 3,942 additional norm-12 points:
**506,709** in total. Both replays check all 23,871 latitude cell
pairs on one side at squared quotient distance at least 12, and all
1,116,462 latitude-to-equator pairs at squared quotient distance at
least eight. The height coordinate supplies the other four, and
opposite heights alone supply squared distance 16.

Run with Python 3; the second replay also needs NumPy.

```sh
python verifications/improved/dim36-37-ct-cell-rebuild/verify.py
python verifications/improved/dim36-37-ct-cell-rebuild/audit.py
```

The Gram data trace to the [Nebe–Sloane K12 entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/K12.html),
as in the preceding Coxeter–Todd construction.
