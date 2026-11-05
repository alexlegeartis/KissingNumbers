# Three-flat obstruction in the 17D mixed-tier family

This is a restricted-family result, not a global upper bound on the 17D
kissing number or a new lower bound. It concerns the same fixed equator and
poles, complete 32-point RM(1,4) orbits, and allowed flat and short/dual tiers
as `research/mixed_tier17_20261104/README.md`.

For exactly three compatible flat cosets in one hemisphere, any compatible
selection of the remaining tier orbits has size at most 19. Consequently
these layers contain at most `(3 + 19) * 32 = 704` points, equal to the
known cap. This remains true without imposing any constraints from the
opposite hemisphere.

The 168 normalized three-flat choices form one verified symmetry class.
For a representative, `fractional_cover.json` gives exact nonnegative integer
weights on cliques of mutually incompatible tier orbits. Every eligible orbit
has total covering weight at least 1,000,000; the sum of all weights is
19,555,639, strictly less than 20,000,000. Thus the weighted clique-cover
inequality gives an independence number strictly below 20, hence at most 19.
The weights and cliques are checked directly with integer arithmetic; the
linear program used to find them is not part of the proof.

The verifier reconstructs the 1,920 tier orbits and the flat and incidence
graphs from the existing `octads18.json`, checks every clique and coverage
weight, and verifies the signed-coordinate symmetry actions and all 168 cases.
Run from the repository root with NumPy:

`OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_three_flat_20261105/verify.py`

This says nothing about zero, one or two flat cosets, partial orbits, changing
the equator, other heights, or constructions outside this family.
