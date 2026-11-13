# 31D annealed joint design: 251,420

This complete coordinate certificate improves the retained [251,396 joint design](../dim31-joint-design/) by 24 points and the pinned PackingStar 251,254 source by 166. It is a verified lower bound, not a claim of optimality or an exhaustive priority search.

## Idea

Jointly choose the 70 transversal blocks of the Leech–Sloane module, then vary its orientation by exact E8 automorphisms. Price **all exterior deletions**, rather than optimizing only the internal module. An exact bipartite matching selects a compatible union of the source and the proposed module. The frozen result inserts 2,218 points and removes 2,052 from the PackingStar source. Annealing is only a discovery heuristic; the proof does not require its replay or an optimality claim.

## Independent proof

```sh
OPENBLAS_NUM_THREADS=1 python verify.py
```

Requires Python and NumPy. The verifier reconstructs the exchange from the supplied source, checks every added-point interaction using exact int32 arithmetic, and independently verifies the entire final configuration: 251,420 distinct points, squared norm 128 and distinct inner products at most 64. The stored coordinates have denominator 2. Binary64 matrix multiplication is exact here because every integer partial sum is below its exact representability threshold; no numerical tolerance is used.

The full geometry check covers all **31,605,882,490** unordered pairs. Source validity is not assumed for final validity. `data/fullverify251420.json` records the frozen coordinate hash; `final.py` reads the actual artifact cardinality for the repository audit. `data/module251420.npz` preserves the discovery module and is not needed to trust the final proof.

`check_module_recipe.py` independently checks the frozen 2,240-point ingredient, all its integer inner products, and its exact rational orthogonal orientation. It is also invoked by `verify.py`.

`verify_geometry.py` implements the complete bounded-integer pair verification used by `verify.py`.
