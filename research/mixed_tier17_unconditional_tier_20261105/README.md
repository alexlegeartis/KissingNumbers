# Complete full-orbit obstruction for the 17D mixed-tier family

This proves **optimality only within the specified construction family**. It
is not an upper bound for the 17D kissing number and not a new lower bound.
The family is the one defined in `research/mixed_tier17_20261104/README.md`:
the 4,320-point 16D equator and two poles are fixed, the two hemispheres use
the specified flat and short/dual heights and candidates, and all chosen
hemisphere pieces are complete 32-point RM(1,4) sign orbits. Within this
family, each hemisphere contains at most 704 such points, and the known
construction attains `4320 + 2*704 + 2 = 5730` points.

## The new global tier inequality

Let `m` be the number of chosen short tier orbits and `U` the union of their
six-element supports. Let `d` be the number of chosen dual tier orbits.
The exact verifier checks the following relations for every candidate:

- Two short orbits whose supports intersect in at least three coordinates
  conflict, including the two orbits with the same support.
- All dual orbits with the same distinguished coordinate conflict.
- A short orbit conflicts with every dual orbit whose distinguished
  coordinate belongs to its support.

Thus the chosen short supports have pairwise intersections at most two, and
`d <= 16 - |U|`. For each chosen support `S`, form its centered incidence
vector `1_S - (6/16) 1` in R^16. It has squared norm 15/4, and two such
vectors have inner product at most -1/4. Nonnegativity of the squared norm
of their sum gives `m*(16-m)/4 >= 0`, hence `m <= 16`.

If `m >= 3`, the union of any three supports has size at least 12, so
`m + d <= 16 + 4 = 20`. If `m = 2`, the union has size at least 10 and
`m + d <= 8`. If `m = 1`, `m + d <= 11`; if `m = 0`, `d <= 16`.
Therefore **the full short-plus-dual tier has at most 20 orbits**, before
any flat-coset constraints are imposed.

## Every flat-coset count

Write `k` for the number of flat orbits. The exact flat conflict graph is
35-regular on 64 vertices and has least eigenvalue at least -5; the verifier
checks its integer annihilating polynomial. Hoffman's bound gives `k <= 8`.
For the remaining cases:

- `k = 0,1,2`: the new global tier bound gives `k + tier <= 22`.
- `k = 3`: the independently checked weighted clique cover in
  `research/mixed_tier17_three_flat_20261105` gives `tier <= 19`.
- `k = 4`: the independently checked clique covers in
  `research/mixed_tier17_four_flat_20261105` give `tier <= 18`.
- `k = 5`: the existing five-flat certificate gives `tier <= 16`; `k = 6`
  inherits that bound from any five-flat subset.
- `k = 7`: the existing seven-flat check gives `tier <= 1`; `k = 8` inherits
  that bound from any seven-flat subset.

The verifier reconstructs the complete candidate family from the existing
`octads18.json`, checks the three structural conflict claims, checks the
flat-graph polynomial and five-/seven-flat certificates, and invokes the
independent three-/four-flat verifiers. No floating-point solver status is
used as proof. Run from the repository root:

`OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_unconditional_tier_20261105/verify.py`

The conclusion does not cover partial sign orbits, changes to the equator,
other heights, extra candidate families or arbitrary 17D configurations.
Those remain legitimate directions for a new SOTA.
