# Four-point orbit search in 17D

This experiment further splits each order-16 sign orbit into four pieces,
each an orbit of a common order-4 subgroup. It jointly searches the two caps
while retaining the 4,320-point equator and two poles. This is a restricted
search, not a new kissing bound or a no-go result.

Each cap has 15,872 orbit candidates, with 19,489,536 same-cap conflict
edges and 16,896 opposite-cap edges. The 5,730-point seed has 352 cap orbits.
`verify.py` independently expands every selected orbit and checks all norms,
distinctness and pairwise integer Gram inequalities. It exactly replays the
seed; the seed replay is not a rediscovery. `audit_graph.py` compares 2,000
sampled graph decisions with direct 4×4 coordinate checks.

Run from the repository root after preparing the existing order-16 data:

    python research/mixed_tier17_eval_subgroup_20261105/prepare.py
    python research/mixed_tier17_order8_20261105/prepare_negative.py
    mkdir -p research/mixed_tier17_order4_20261105/cache
    g++ -O3 -std=c++17 research/mixed_tier17_order4_20261105/build_graph.cpp -o /tmp/eval2_graph
    g++ -O3 -std=c++17 research/mixed_tier17_order4_20261105/cross_graph.cpp -o /tmp/eval2_cross
    /tmp/eval2_graph && /tmp/eval2_cross
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order4_20261105/verify.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order4_20261105/audit_graph.py
    python research/mixed_tier17_order4_20261105/prepare_anneal.py
    g++ -O3 -std=c++17 research/mixed_tier17_order4_20261105/anneal.cpp -o /tmp/eval2_anneal
    /tmp/eval2_anneal 170429
    python research/mixed_tier17_order4_20261105/decode_anneal.py
    OPENBLAS_NUM_THREADS=1 python research/mixed_tier17_order4_20261105/verify.py research/mixed_tier17_order4_20261105/cache/anneal_POTENTIAL.npz

The annealer can replace pieces in both caps simultaneously and temporarily
lose orbits. An output is a proposal until independently replayed and compared
with the latest literature. A finite search failure does not prove optimality.
