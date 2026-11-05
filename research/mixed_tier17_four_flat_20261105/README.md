# Four-flat obstruction in the 17D mixed-tier family

This extends the repository's earlier five-flat and seven-flat checks to the
case of exactly four complete flat cosets in one hemisphere. It is a
**restricted-family obstruction**, not an upper bound on the 17D kissing
number and not a new lower bound.

The family is the one described in `research/mixed_tier17_20261104/README.md`:
the 16-dimensional equator and poles are fixed; flat and short/dual tier
points are selected as complete 32-point orbits under RM(1,4). For every
compatible choice of exactly four flat cosets, the remaining eligible short
and dual tier orbits have independence number at most 18. Consequently that
hemisphere has at most `(4 + 18) * 32 = 704` points in these layers, equal
to the known baseline cap. The proof does not require or assume the opposite
hemisphere's compatibility; adding that constraint can only reduce choices.

After translating a flat coset to zero, there are exactly 336 four-coset
choices. Verified signed coordinate permutations partition them into two
symmetry classes of sizes 280 and 56. One representative's eligible tier
orbits admit a 16-clique cover; the other admits an 18-clique cover. Each
clique consists of mutually incompatible 32-point tier orbits, so any
compatible selection contains at most one per clique. The JSON files hold
these two small witnesses. The verifier reconstructs RM(1,4), RM(2,4), every
candidate tier orbit, the flat graph and flat–tier conflicts from the existing
`octads18.json`, checks the cliques by exact integer dot products, validates
all symmetry actions, and checks coverage of all 336 cases. No solver status
or floating-point numerical bound is part of the proof.

Run from the repository root with NumPy installed:

`OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_four_flat_20261105/verify.py`

The obstruction says nothing about fewer than four flat cosets, partial orbits,
changing the equator, new heights or other construction families.
