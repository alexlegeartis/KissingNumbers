# 120-degree equator sums and joint 17D exchanges

This construction search adds an algebraic/geometric candidate family to the
joint equator-and-cap orbit optimizer. For equator vectors u and v of squared
norm 288 with u·v = -144, their sum has squared norm 288. The pool generator
forms such sums between a transformed-equator candidate and a seed equator
point, then closes each sum under the common order-4 sign subgroup. These
integer vectors are proposals; validity relative to the remaining layers must
be established separately.

The initial seeded pool contains 30,792 new subgroup orbits of size two or
four. Their source-blocker counts are computed exactly across the whole
5,730-point seed. The joint finite graph retains the 2,144 sum orbits with at
most 36 source-point blockers, alongside 2,560 transformed-equator and 4,032
cap candidates. Together with 1,614 old subgroup orbits, this gives 10,350
weighted vertices. Both stochastic annealing and bounded shared-blocker MILP
searches may add equator and cap orbits jointly and remove old orbits from any
layer. The graph uses exact integer orbit conflict checks, but its candidate
pool and both optimizers are restricted; neither failure nor a solver status
is a global result.

Prepare the prerequisites in `mixed_tier17_joint_layers_20261105/README.md`.
From the repository root:

    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_root_sums_20261105/build_root_pool.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_root_sums_20261105/build_graph.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_root_sums_20261105/audit_graph.py
    g++ -O3 -std=c++17 research/mixed_tier17_root_sums_20261105/anneal.cpp -o /tmp/root17
    /tmp/root17 170701
    OPENBLAS_NUM_THREADS=1 ROOT_TRIALS=30000 python research/mixed_tier17_root_sums_20261105/lns.py
    python research/mixed_tier17_root_sums_20261105/decode.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/verify.py research/mixed_tier17_root_sums_20261105/cache/POTENTIAL.npz

The audit compares 2,000 graph decisions with direct expanded-orbit Gram
checks. A solver output is only a proposal until the independent full-coordinate
verifier succeeds and the count is compared against current literature. No
new SOTA is asserted here.
