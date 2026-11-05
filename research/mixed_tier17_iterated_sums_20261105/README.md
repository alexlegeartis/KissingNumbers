# Iterated 120-degree sums in the 17D joint search

This experiment repeats the norm-preserving construction u+v for equator
vectors with u·v = -144, now using first-generation sum orbits as one input.
The resulting 20,642-orbit second-generation pool includes 10,149 orbits not
in either the transformed-equator or first-generation pools. The graph adds up
to 2,000 new second-generation orbits with at most 36 source-point blockers to
the previous joint search, yielding 12,350 weighted source/equator/cap/sum
choices. Exact source and mutual conflict checks form the graph; stochastic
and bounded-MILP searches propose exchanges, not automatically valid records.

Prepare the neighboring `mixed_tier17_root_sums_20261105` pool and its
prerequisites first. Then from the repository root:

    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_iterated_sums_20261105/build_pool.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_iterated_sums_20261105/build_graph.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_iterated_sums_20261105/audit_graph.py
    g++ -O3 -std=c++17 research/mixed_tier17_iterated_sums_20261105/anneal.cpp -o /tmp/sums2
    /tmp/sums2 170801
    OPENBLAS_NUM_THREADS=1 ROOT2_TRIALS=30000 python research/mixed_tier17_iterated_sums_20261105/lns.py
    python research/mixed_tier17_iterated_sums_20261105/decode.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_equator_exchange_20261105/verify.py research/mixed_tier17_iterated_sums_20261105/cache/POTENTIAL.npz

The graph audit checks 2,000 sampled decisions against expanded integer Gram
calculations. All proposals require full-coordinate exact verification and
current literature comparison; failure of a bounded search is not a no-go
theorem. No new bound is asserted.
