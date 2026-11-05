# Joint order-eight orbit search in 17D

This experiment breaks each order-16 sign orbit from the existing 17D
construction into two order-8 orbits. It searches changes to both hemispheres
jointly, keeping the 4,320-point equator and the two poles fixed. It is a
search, not a new lower bound or an optimality claim.

There are 256 flat and 7,680 short/dual orbit vertices in each cap. The
same-cap graph has 8,366,496 conflict edges; the two-cap graph also has 8,448
cross-cap edges. A conflict is computed by exact integer sign-orbit
comparisons, including different same- and opposite-height thresholds. The
baseline uses 176 order-8 cap orbits, for 4,322 + 176×8 = 5,730 points.

The local search chooses a group of source orbits to remove, then uses a MILP
to select mutually compatible new orbit pieces while accounting for shared
blockers. It can change both caps together. It explores only sampled
neighborhoods; time limits and failure do not prove impossibility. Any
`POTENTIAL.npz` must be passed to `verify.py`, which independently expands all
selected orbits, equator and poles and checks every integer Gram inequality.
A verified count must also be compared to current literature before claiming a
record. The verifier's baseline replay is not a rediscovery.

From the repository root:

    python research/mixed_tier17_eval_subgroup_20261105/prepare.py
    python research/mixed_tier17_order8_20261105/prepare_negative.py
    g++ -O3 -std=c++17 research/mixed_tier17_order8_20261105/build_graph.cpp -o /tmp/eval3_graph
    g++ -O3 -std=c++17 research/mixed_tier17_order8_20261105/cross_graph.cpp -o /tmp/eval3_cross
    /tmp/eval3_graph && /tmp/eval3_cross
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order8_20261105/verify.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order8_20261105/audit_graph.py
    OPENBLAS_NUM_THREADS=1 EVAL3_TRIALS=30000 python research/mixed_tier17_order8_20261105/search.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order8_20261105/verify.py research/mixed_tier17_order8_20261105/cache/POTENTIAL.npz

Generated graph and proposal files stay in `cache/` and need substantial memory.
The known seed was exactly replayed, and 2,000 random same/opposite-cap orbit
pairs were compared against direct 8×8 integer Gram checks. Neither test proves
that the search will find an improvement.
