# Native partial-design modules in 29D and 30D

13 November 2026. Complete explicit certificates give:

| Dimension | Previous repository | This construction | Gain | PackingStar November 11 | Gain over that source |
|---|---:|---:|---:|---:|---:|
| 29 | 231,776 | **231,814** | +38 | 231,774 | +40 |
| 30 | 248,084 | **248,180** | +96 | 247,992 | +188 |

PackingStar is pinned to [commit 2f3453b](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111). The previous repository certificates and their provenance remain in [November sextets](../dim28-29-november-sextets/) and [30D Leech–Sloane section](../dim30-leech-sloane-section/). These comparisons do not establish exhaustive worldwide priority or global optimality.

## Idea: optimize after taking the section

A coordinate section of a complete fifteen-dimensional design inherits unnecessary extension requirements. We discard those requirements and optimize a native partial quadruple design in the smaller dimension.

Coordinates are grouped into pairs. A transversal quadruple chooses one coordinate from each of four different pairs. Each chosen quadruple contributes all sixteen signed points. It is enough that no triple occur in two chosen quadruples; a complete fifteen-dimensional Steiner system is not required. We keep the other points of the parent section fixed and price conflicts with every exterior source point.

* **30D:** the old fourteen-dimensional local section had 876 fixed points and 30 transversal quadruples, for 1,356 points. A six-pair partial design admits 36 compatible quadruples in the selected orientation: **876 + 36×16 = 1,452** points. Its optimal union with the old certificate adds 384 and removes 288, gaining 96.
* **29D:** the thirteen-dimensional section has 568 fixed points. A five-pair partial design uses 16 quadruples, yielding **568 + 16×16 = 824** points. Whole replacement alone is not competitive. An additional exact E8/tail orientation followed by a bipartite independent-union exchange adds 492 and removes 454, gaining 38. Thus both the native design and retaining compatible pieces of the old configuration matter.

The local eight-dimensional coordinates are embedded on the Golay octad `{0,2,4,5,6,10,11,23}`; the remaining local coordinates occupy the ambient tail. Every local and ambient point has squared norm32 before the final portable scaling.

## Reproduce

Python and NumPy suffice. No network, optimizer, or assumption that the source configuration is valid is needed for verification:

```bash
python final.py
OPENBLAS_NUM_THREADS=2 taskset -c 0-126 python verify.py
# Or select one dimension:
OPENBLAS_NUM_THREADS=2 taskset -c 0-126 python verify.py --dimension 29
```

Adjust CPU affinity for the host and reserve one logical CPU.

`verify_recipe29.py` and `verify_recipe30.py` independently rebuild the fixed parent sections, verify the selected partial designs and rational orientations, and check every local pair using int64. `verify.py` checks hashes, exact source/delta reconstruction, module membership of all additions, and every changed pair by direct int64 coordinate sums. It then verifies the complete final geometry using exact bounded integer products in binary64, without tolerance:

* 29D: **26,868,749,391** unordered point pairs;
* 30D: **30,796,532,110** unordered point pairs.

Portable coordinates have squared norm128 and distinct-pair products at most64. Actual antipodal pairs and unpaired points are counted rather than assumed.

## Provenance and limitations

`data/` contains full final coordinates, source coordinates, norm32 exchanges, the original parent recipe, selected blocks and exact orientations. `discovery/` retains development scripts and logs; those scripts document the original research workspace and are not the portable verifier interface.

The fixed-orientation complete-support MILPs reported matching numerical objective and bound. This is solver evidence about those finite models, not a rigorous global optimality certificate. The resulting configurations are proved independently of the solver by the explicit geometric verifiers. Neither the kissing numbers nor these larger construction families are claimed optimal.

`negative_controls.py` rejects four malformed witnesses: wrong norm, duplicate, illegal angle and wrong normalization metadata. `parent_recipe.py` supplies the shared exact reconstruction of the inherited fixed parent, and `verify_geometry.py` implements complete final geometry. `final.py` reads counts from the actual coordinate arrays. Run verification without Python's `-O` flag; entry points explicitly reject optimized assertion removal.

Discovery scripts retained for provenance are `discovery/native29.py`, `discovery/native30.py`, `discovery/replay_native30.py`, and `discovery/union_native29.py`. They are not required for verification and retain research-workspace imports; use the portable proof commands above instead.

The independently written `independent29/audit.py` and `independent30/audit.py` check all final pairs directly without antipodal grouping, and independently reconstruct the source exchanges. Their receipts and detailed assumptions are documented in the corresponding subdirectory READMEs. Both are included in the repository SOTA runner alongside the main verifier and malformed-witness controls.
