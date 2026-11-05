# Joint 17D equator and cap replacement

This experiment allows new orbit pieces in both the equator and the caps,
while permitting the removal of any old subgroup orbit in any layer. It uses
the common order-4 sign subgroup. It is a finite-pool heuristic search, not a
new lower bound or a proof of optimality.

The current pool comprises 1,614 source orbits, 2,560 transformed equator
orbits with at most 56 source-point blockers, and 4,032 new cap orbits with
at most eight old-cap-orbit blockers. All new orbits have four points; source
orbit weights are their exact sizes. The graph has 8,206 vertices and 880,768
exact conflict edges. Its stochastic weighted independent-set search can
simultaneously add equator and cap pieces, remove shared blockers once, and
temporarily lose points. A graph proposal is not a bound: `decode.py` expands
it and the independent full-coordinate verifier checks every pair.

First prepare the seed and order-4 graphs as described in the neighboring
`mixed_tier17_order4_20261105` README, and build the equator pool using
`mixed_tier17_equator_exchange_20261105/build_pool.py`. Then, from the repo root:

    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_joint_layers_20261105/prepare_cap_pool.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_joint_layers_20261105/build_graph.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_joint_layers_20261105/audit_graph.py
    g++ -O3 -std=c++17 research/mixed_tier17_joint_layers_20261105/search.cpp -o /tmp/joint17
    /tmp/joint17 170601
    python research/mixed_tier17_joint_layers_20261105/decode.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/verify.py research/mixed_tier17_joint_layers_20261105/cache/POTENTIAL.npz

The seed replay has exactly 5,730 points. The graph audit samples 2,000 edges
and nonedges and compares them with direct subgroup-orbit Gram calculations;
it is a regression check rather than a full proof of graph correctness. Any
proposed count requires full exact replay and a current literature comparison.
