# Verified asymmetric exchanges in dimensions 26, 28, 29 and 30

28 October 2026. These four integer configurations exceed the checked
PackingStar coordinate files. Together with the separate 27D and 31D
packages, every dimension from 26 through 31 exceeds that pinned comparison.

| Dimension | PackingStar | This construction | Gain | Removed | Added | Unpaired points |
|---|---:|---:|---:|---:|---:|---:|
| 26 | 201,134 | **201,136** | 2 | 59 | 61 | 58 |
| 28 | 221,664 | **221,694** | 30 | 130 | 160 | 32 |
| 29 | 231,636 | **231,638** | 2 | 14 | 16 | 22 |
| 30 | 247,956 | **247,960** | 4 | 80 | 84 | 104 |

PackingStar's retained points are credited as the baseline. The comparison
is pinned to commit
[`b67832ad9970895a22ab667bb06d1d43a45c6b18`](https://github.com/CDM1619/PackingStar/tree/b67832ad9970895a22ab667bb06d1d43a45c6b18/25D-31D/new_bounds_20261026).
Our contribution is the explicit asymmetric replacement in each dimension.
These are lower bounds, not claims of optimality.

## Complete independent verification

```bash
python final.py
python verify.py
python verify.py 26 30
```

Python 3.9+ and NumPy suffice. Run without `-O`. The default verifier checks
all four complete configurations and does not assume that the public input,
the optimizer, or a search log is correct. The final replays took roughly
70-122 seconds per dimension on the discovery machine, with some running
concurrently. Reports are supplied under `data/verification*.json`.

Every point has squared norm 32, and every distinct pair has inner product
at most 16. Dividing by sqrt(32) produces unit vectors with pairwise inner
products at most 1/2. The NPZ files use arrays `rational` and `radical` in
the format A+B sqrt(2); B is identically zero in these four certificates.

The verifier checks norms and distinctness, then identifies actual opposite
pairs. Between two such pairs it checks an absolute inner-product bound;
between an opposite pair and an unpaired point it also checks an absolute
bound. Between two unpaired points it checks the signed product only.
Antipodality is never assumed. With P opposite pairs and U unpaired points,
the checks cover exactly

    4*binom(P,2) + P + 2*P*U + binom(U,2) = binom(2*P+U,2)

unordered point pairs, including the P opposite pairs whose product is -32.
Small integer dot products are evaluated exactly by binary64 BLAS. Bounds
on all coefficients and the dimension ensure that each product and partial
sum is exactly representable. The generic checker also handles quadratic
coefficients using sign-sensitive integer comparisons. No tolerance is used.

## Construction and search scope

Let I = {0,2,4,5,6,10,11,23}. Every changed point is supported on I and
the d-24 tail coordinates. The unchanged complement is retained from the
pinned public file. The normal-sector counts change as follows:

| Dimension | Previous sector | New sector |
|---|---:|---:|
| 26 | 334 | 336 |
| 28 | 648 | 678 |
| 29 | 728 | 730 |
| 30 | 1,216 | 1,220 |

The finite candidate dictionaries combine E8 heads with tail vectors of
complementary squared norm. The search permits each point and its negative
to be chosen independently, using signed conflicts. The final configurations
therefore escape the antipodal restriction of the preceding searches.
Candidate pools were enlarged and then recomputed around the improved codes;
the latter step yielded further gains in 26D and 28D.

The initial global low-conflict models and subsequent neighborhood solves
do not establish an unrestricted maximum. In particular, the final
recentered low-conflict pools in 28D and 30D were solved to their current
counts, whereas the corresponding 26D and 29D models retain unresolved
gaps. Larger conflicts, different coordinates and changes outside the
normal sector remain open. The coordinate verification is independent of
all such optimization statements.

The final 28D and 30D finite-pool optima also have independent certificates.
The 928-point 28D pool is partitioned into 678 conflict cliques; the
1,642-point 30D pool is partitioned into 1,220 conflict cliques. A compatible
subset can use at most one point from each clique, and the supplied selected
subsets attain those bounds. `verify_pool.py` checks the entire covers and
the matching subsets with integer arithmetic, without an optimizer. The
default replay checks these certificates too. These bounds concern only the
explicit supplied pools, not all possible normal sectors or kissing codes.

The earlier [27D quadratic module](../dim27-quadratic-normal-module/) gives
209,600, and the [31D normal-module exchange](../dim31-normal-module-exchange/)
gives 251,110. Their coordinates and proofs are unchanged.

## Rebuild from the pinned public files

Each `data/exchange<D>.npz` contains the exact removed and added integer
points. For example:

```bash
python rebuild.py 26 D26_201134_integer_coordinates.npz --output rebuilt26.npz
python verify_core.py rebuilt26.npz
```

Repeat with dimensions 28, 29 or 30 and the corresponding pinned input.
The rebuild verifies the public file's SHA-256, removed-point membership,
the resulting count, distinctness, and the canonical coordinate hash.
`data/provenance.json` records these hashes and the source URLs. No solver,
network access or private search state is needed after obtaining the input.
