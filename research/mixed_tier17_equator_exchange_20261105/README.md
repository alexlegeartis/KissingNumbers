# Coupled equator and cap exchanges in 17D

This is a discovery search, not a new lower bound. It uses the known 5,730-point
seed and the order-4 sign subgroup from the related order-4 search. It applies
24 seeded coordinate permutations that preserve that subgroup to the equator,
then decomposes the resulting new points into complete subgroup orbits. The
initial pool contains 13,568 new four-point equator orbits. The entire original
configuration, including both caps, decomposes into 1,614 subgroup orbits of
various sizes.

The search jointly selects new equator orbits and removes their blockers from
any original layer. A shared old orbit is charged only once even if many new
orbits need it removed. For bounded neighborhoods, a MILP maximizes added
points minus removed points and enforces exact candidate–candidate conflicts.
The pool is finite, the neighborhoods are sampled, and solver failures are not
proofs of optimality. Adding new cap orbits is not yet included. Any proposal
must pass the independent full-coordinate verifier and a current literature
comparison before it is called a record.

From the repository root, first prepare the seed as described in the related
`mixed_tier17_order4_20261105` README. Then:

    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/build_pool.py
    OPENBLAS_NUM_THREADS=1 EQUATOR_TRIALS=10000 python research/mixed_tier17_equator_exchange_20261105/search.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/verify.py research/mixed_tier17_equator_exchange_20261105/cache/POTENTIAL.npz

Pool generation independently replays the original seed and exactly checks
that each candidate orbit's point blockers are unions of source subgroup
orbits. The resulting cached pool is under one megabyte; the candidate search
uses SciPy MILP and checks a proposed positive exchange separately.
