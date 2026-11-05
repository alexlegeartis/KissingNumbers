# Two-flat obstruction in the 17D mixed-tier family

This is a restricted-family result, not a global upper bound on the 17D
kissing number or a new lower bound. It concerns the fixed equator and poles,
complete 32-point RM(1,4) orbits, and flat and short/dual tier candidates
specified in `research/mixed_tier17_20261104/README.md`.

For exactly two compatible flat cosets in one hemisphere, any compatible
selection of the remaining tier orbits has size at most 20. Thus these layers
contain at most `(2 + 20) * 32 = 704` points, equal to the known cap, even
before imposing any opposite-hemisphere constraints.

The 28 normalized two-flat choices form one verified symmetry class. For one
representative, `fractional_cover.json` gives integer weights on cliques of
mutually incompatible tier orbits. Every eligible orbit receives weight at
least 1,000,000, while the total is 20,233,717, strictly less than 21,000,000.
The weighted clique-cover inequality therefore gives an independence number
strictly less than 21, hence at most 20. Every clique and integer coverage
inequality is checked directly; the linear program used to find the witness
is not part of the proof.

The verifier reconstructs the 1,920 tier orbits and flat and incidence graphs
from the existing `octads18.json`, checks the witness, verifies signed-coordinate
symmetries, and covers all 28 normalized cases. Run from the repository root:

`OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_two_flat_20261105/verify.py`

This does not resolve zero or one flat coset, partial orbits, changing the
equator, other heights, or constructions outside this family.
