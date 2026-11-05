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

A second search builds the exact finite conflict graph containing every source
orbit plus the 2,560 new equator orbits with at most 56 source-point blockers.
Its 4,174 weighted vertices and 112,000 edges support a stochastic weighted
independent-set search over removals and additions in all original layers:

    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/build_full_graph.py
    g++ -O3 -std=c++17 research/mixed_tier17_equator_exchange_20261105/weighted_anneal.cpp -o /tmp/equator_anneal
    /tmp/equator_anneal 170501
    python research/mixed_tier17_equator_exchange_20261105/decode_weighted.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/verify.py research/mixed_tier17_equator_exchange_20261105/cache/weighted_POTENTIAL.npz

This remains a finite candidate pool and stochastic search; it can also remove
old cap orbits but does not add new cap orbits. The full-coordinate verifier,
not the graph or annealer alone, is required to establish a kissing code.
`audit_full_graph.py` compares 2,000 sampled edges and nonedges of this graph
with direct integer Gram checks on the expanded point orbits.
