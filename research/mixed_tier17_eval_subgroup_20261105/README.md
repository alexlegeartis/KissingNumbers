# 17D search with an evaluation-index subgroup

This experiment breaks full RM(1,4) sign orbits into orbits of its index-two subgroup H = {h : h(0)=0}. This subgroup is inequivalent under affine coordinate permutations to the derivative-index subgroup used in earlier experiments. Its order-16 orbit family has 128 flat and 3,840 short/dual tier vertices. The previously proved full-group cap of 22 orbits does not apply here. Same-cap conflicts form an exact integer graph. This is a search for a new construction, not a proof that one exists.

From the repository root, with Python, NumPy and a C++17 compiler:

    python research/mixed_tier17_eval_subgroup_20261105/prepare.py
    g++ -O3 -std=c++17 research/mixed_tier17_eval_subgroup_20261105/coupling.cpp -o /tmp/eval17_coupling
    /tmp/eval17_coupling
    g++ -O3 -std=c++17 research/mixed_tier17_eval_subgroup_20261105/build_graph.cpp -o /tmp/eval17_graph
    /tmp/eval17_graph
    python research/mixed_tier17_eval_subgroup_20261105/verify_candidate.py
    g++ -O3 -std=c++17 research/mixed_tier17_eval_subgroup_20261105/tabu.cpp -o /tmp/eval17_tabu
    /tmp/eval17_tabu

The first verification replays the 5,730-point supplied seed, rather than discovering a new bound. The tabu program attempts to find 45 compatible positive-cap orbits versus 44 in the seed. Its file `cache/potential.txt`, if produced, is only a same-cap candidate. Run `verify_candidate.py --selection research/mixed_tier17_eval_subgroup_20261105/cache/potential.txt` to expand all orbits and independently check every pair against the fixed negative cap, equator and poles. A failure is not a kissing configuration. All computations use exact integers; the verifier uses norm squared 288 and maximum inner product 144.

As of 2026-11-05, this experiment has not produced a verified improvement. Generated files under `cache/` are not needed in version control. The small files in `seed/` contain the supplied equator and opposite-cap orbit representatives; the verifier expands them and checks the complete configuration independently.
