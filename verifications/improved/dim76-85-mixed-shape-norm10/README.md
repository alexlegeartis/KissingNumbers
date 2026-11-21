# Gamma72 mixed-shape norm-ten layers in 76–85 dimensions

This is a kissing-number lower-bound construction, not a bound on packing density or an optimality claim. The shell additions are certified by a finite independence bound, not by listing every member of the independent sets.

## Bounds

| Dimension | Previous repo bound | New lower bound | Gain | Tail shape | Head scale squared |
|---:|---:|---:|---:|:---|:---|
| 76 | 6,218,315,288 | 6,218,371,050 | 55,762 | 2,1,1 | 8/13 |
| 77 | 6,218,490,542 | 6,218,669,616 | 179,074 | 2,2,1,1 | 16/25 |
| 78 | 6,218,688,770 | 6,218,991,624 | 302,854 | 2,2,1,1 | 16/25 |
| 79 | 6,219,067,454 | 6,219,415,042 | 347,588 | 2,1,1,1,1 | 12/19 |
| 80 | 6,219,714,034 | 6,219,960,348 | 246,314 | 2,1,1,1,1 | 12/19 |
| 81 | 6,220,589,740 | 6,220,773,344 | 183,604 | 1,1,1,1 | 4/7 |
| 82 | 6,221,706,936 | 6,222,087,772 | 380,836 | 1,1,1,1 | 4/7 |
| 83 | 6,223,067,174 | 6,223,765,754 | 698,580 | 1,1,1,1 | 4/7 |
| 84 | 6,224,663,766 | 6,225,860,206 | 1,196,440 | 1,1,1,1 | 4/7 |
| 85 | 6,226,484,034 | 6,228,414,172 | 1,930,138 | 1,1,1,1 | 4/7 |

Replay with `python verifications/improved/dim76-85-mixed-shape-norm10/verify.py` (NumPy required; no LP solver). Supply selected dimension numbers as optional arguments. The supplied rational duals, not the numerical optimization used to discover them, certify every common-neighbor bound.

## Geometry

Let (k=n-72), let (S) be the full norm-ten shell of Gamma72, and let (W) contain every signed permutation of the indicated shape, with zeros appended in `k` coordinates. Put (w=|r|^2) for (rin W) and (q=8-10x) for the listed head scale (x). Each candidate is the antipodal pair represented by

    (sqrt(x) v, sqrt(q/w) r),  v in S, r in W.

The graph on pairs has an edge for `abs(x*(v.u)+(q/w)*(r.s)) > 4`. Thus an independent set gives two kissing points per vertex at squared norm eight. All 6,218,175,600 Gamma72 minima may be kept on the equator: their absolute dot product with any unscaled norm-ten vector is at most five by the minimum norm, and (25xle16). The pure tail consists of the `2k(k-1)` scaled Dk roots (2(\pm e_i\pm e_j)), and the maximum possible mixed product is at most twice the sum of the two largest shape coordinates times `sqrt(q/w)`. The verifier checks its square is at most 16. The three layers are disjoint.

The 81–85D choice (x=4/7) makes the graph condition the integer comparison `abs(v.u+r.s)>7`. With distinct head lines, only head products (j=4,5,6) can contribute edges. For (j=6) this latitude needs tail product at least two instead of at least one in the previous (8/13) construction, removing many edges. In 76–80D, larger signed-permutation orbits with unequal amplitudes work better.

## Rigorous counting

The published premises are [Nebe's even unimodular Gamma72 with minimum eight](https://arxiv.org/abs/1008.2862) and [Venkov's 11-design theorem for these shells](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf). The checker recomputes the weight-36 theta coefficients to norm ten and the complete one-dimensional shell populations using exact fractional Lagrange interpolation. It does not reconstruct the lattice itself.

Degree is constant: the shell moments give the same marginal populations about any head, and signed permutations act transitively on each full tail shape. For a neighbor on a distinct head line, orient it to have positive head product (j). Permutations preserving the reference tail and sign changes on its zero coordinates group neighbor tails by their exact orbit keys. For equal-amplitude weight four, this is just overlap and signed dot. For general shapes the key records the multiset of signed amplitudes within each equal-amplitude reference group, and the multiset of absolute amplitudes in the zero group.

For any such pair of vertices, the remaining oriented norm-ten head vectors are indexed by two integral head products. Gram positivity and the lattice minimum discard impossible profiles. A finite rational majorant of the common-neighbor tail count is then verified on every profile. Its expectation is known exactly from the two-coordinate spherical moments through degree ten and the two exact marginals. The exceptional heads on the two distinguished lines are counted directly, dividing by two everywhere to quotient simultaneous negation of head and tail. Same-head-line edges instead use the trivial common-neighbor bound of the full degree. This gives a uniform exact triangle upper bound.

The final graph theorem and block inequalities are as in the preceding [triangle-sensitive certificate](../dim79-86-norm10-triangle-layers/PROOF.md), using the [Boßecker–Rautenbach interpolation theorem](https://doi.org/10.1016/j.disc.2010.05.026). Here roughly 512 rational arithmetic blocks propagate the bound on `g_d=(d+1)f(d,T)-1`, and each block rounds its fixed-point integer value downward. This can only decrease the certified lower bound. A final ceiling is valid because independence numbers are integers.

- `cert_4.json` through `cert_13.json`: frozen rational majorants and graph counts.
- `verify.py`: exact proof replay, apart from displayed decimal summaries not used for proof.
- `claims.json`: lower bounds and comparison with the repository immediately before this change.
