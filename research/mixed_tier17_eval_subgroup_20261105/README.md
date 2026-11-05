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

## Simultaneous two-cap walk

After preparing the graph, `two_cap_walk.py` searches both hemispheres at once. It starts from the known 88 order-16 cap orbits, applies exact conflict constraints, and permits a temporary loss of one orbit while forcing a chosen subgroup orbit without its full-group partner. This helps avoid a neutral walk that simply reassembles complete 32-point orbits. It uses SciPy's MILP solver, so the search itself is heuristic and a solver failure is not an impossibility proof.

    EVAL17_SEED=171448 EVAL17_EPOCHS=1000 EVAL17_TRIALS=6000 python research/mixed_tier17_eval_subgroup_20261105/two_cap_walk.py
    python research/mixed_tier17_eval_subgroup_20261105/verify_two_cap.py research/mixed_tier17_eval_subgroup_20261105/cache/two_cap_POTENTIAL.npy

The two-cap graph has two copies of the same-cap conflict graph, plus 4,224 inter-cap edges. Tier-to-tier pairs are automatically compatible across opposite heights. Opposite flat-to-tier conflicts use the exact threshold 16, and opposite flat-to-flat conflicts occur only for matching subgroup cosets. A proposed improvement must still pass the independent full-coordinate verification, which also checks the fixed equator and poles. A sampled broken-orbit state with 87 cap orbits expanded to 5,714 points and passed that exact check; it is below the seed and is not a new lower bound.
