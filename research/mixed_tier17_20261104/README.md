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

## Coupled equator-switch experiment

The `equator_*` programs jointly consider removing groups of odd-Barnes–Wall
equator points and adding six-support tier orbits at both heights. Unlike the
fixed-equator pool, the source scans every six-coordinate support, including
supports outside RM(2,4). The subgroup-orbit seed reconstructs the known
5,730 points, and `verify_seed.py` independently checks all pair products in
integer arithmetic. Each shared equator deletion is charged only once.

The shipped frontier contains 1,376 mirrored candidate orbits with a
single-orbit trade deficit of at most 256 points, after testing individual
orbit compatibility. This is a deliberately restricted finite pool, not all
possible changes to the equator, not a proof that the full search is closed,
and not a new bound. The MILP uses floating-point optimization; any candidate
must be reconstructed and independently checked exactly before it is claimed.
Its objective charges deletions of seed equator, tier and flat orbits and
allows compatible new orbits to share those deletions.

From the repository root, first run `prepare.py` and `coupling.cpp` as in
the first section, then `rank4_prepare.py` as in the partial-orbit section.
The rank-four search and coupling are not prerequisites. Then run:

    python research/mixed_tier17_20261104/equator_seed.py
    python research/mixed_tier17_20261104/equator_prepare.py
    python research/mixed_tier17_20261104/verify_seed.py
    c++ -O3 -std=c++17 research/mixed_tier17_20261104/equator_scan.cpp -o /tmp/mixed17-eq-scan
    /tmp/mixed17-eq-scan
    c++ -O3 -std=c++17 research/mixed_tier17_20261104/equator_graph.cpp -o /tmp/mixed17-eq-graph
    /tmp/mixed17-eq-graph
    MIXED17_MILP_SECONDS=900 python research/mixed_tier17_20261104/equator_milp.py

`equator_anneal.cpp` is an optional heuristic on the generated graph. Neither
a search failure nor a numerical optimization limit proves an upper bound.

`asymmetric_graph.py` expands the same frontier to independent positive- and
negative-height choices, while keeping equator deletion costs shared. This
contains the paired model as a special case. After building the paired graph,
run `asymmetric_graph.py` and `asymmetric_milp.py`; if a positive candidate
is saved, run `verify_asymmetric.py` on it. `verify_switch.py` similarly
checks a paired candidate. The verifiers reconstruct every point and check
all distinct pairs with integer arithmetic. The empty exchange reproduces
the 5,730-point seed in both independent verifier paths.
