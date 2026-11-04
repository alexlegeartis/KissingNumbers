# Coupled 17D orbit search

This is a restricted research result, not a new kissing-number lower bound.
It studies one hemisphere of the odd Barnes–Wall construction, with the
equator and poles fixed. The two norm-six body types are signed six-support
vectors and the half-integer dual type. Both are decomposed into complete
32-point sign orbits under RM(1,4). The flat tier is restricted to complete
RM(1,4) cosets in RM(2,4), also of size 32.

The generator obtains 896 short and 1,024 dual candidate tier orbits, plus
64 flat orbits. An exact coupling table decides compatibility against all
32 sign-group elements. For every compatible set of five flat cosets, the
certificate provides sixteen six-cliques covering all 96 eligible tier
orbits. Consequently the mixed tier has at most sixteen orbits, or 512
points. By translating the flat cosets, it suffices to check the 280 sets
containing zero. Every short and dual vector is included under the specified
orbit restriction. This does not bound arbitrary 17D configurations, partial
sign orbits, changes to the equator, or other norm and height choices.

The separate seven-flat check examines all 56 compatible seven-coset sets
containing zero. For each, it checks that the 16 eligible tier orbits form a
clique, so at most one may be selected. Together with the five-flat result,
this gives a precise obstruction for five or more full flat cosets, but does
not resolve the cases of fewer flat cosets or broken sign symmetry.

From the repository root, run with Python 3, NumPy, and a C++17 compiler:

    python research/mixed_tier17_20261104/prepare.py
    c++ -O3 -std=c++17 research/mixed_tier17_20261104/coupling.cpp -o /tmp/mixed17-coupling
    /tmp/mixed17-coupling
    python research/mixed_tier17_20261104/verify_five_flat.py
    python research/mixed_tier17_20261104/verify_seven_flat.py

The scripts write regenerable binary search inputs in this directory. The
proof certificate is `five_flat_clique_certificate.json`. Both verifiers use
integer comparisons, reconstruct flat compatibility independently, and check
each reported clique directly against the 32 sign operations. No numerical
MILP status is used as a proof.

## Partial-orbit discovery experiment

`rank4_prepare.py` and `rank4_coupling.cpp` split the sign group to an
index-two subgroup. This gives 3,840 tier and 128 flat orbits of 16 points.
`rank4_joint_search.py` proposes partially selected flat cosets and solves
the corresponding one-hemisphere tier refill problem. It is a search, not
an exact upper-bound certificate. Even an improved one-hemisphere output
must be combined with a compatible opposite hemisphere and independently
verified before it yields a 17D lower bound.

    python research/mixed_tier17_20261104/rank4_prepare.py
    c++ -O3 -std=c++17 research/mixed_tier17_20261104/rank4_coupling.cpp -o /tmp/mixed17-rank4
    /tmp/mixed17-rank4
    MIXED17_TRIALS=1500 MIXED17_MILP_SECONDS=3 python research/mixed_tier17_20261104/rank4_joint_search.py

The random seed is fixed. The search samples flat sets near one initial
configuration; it is not exhaustive and its current runs have not found
a new record.
