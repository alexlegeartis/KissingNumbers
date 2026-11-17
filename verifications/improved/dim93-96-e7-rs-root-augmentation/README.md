# Three-alphabet constructions with E7 heads and Reed–Solomon/root augmentation

Approved exact lower bounds: **42,951,379,608 in93D**, **64,426,347,212 in94D**, **107,376,020,548 in95D**, and **206,160,530,240 in96D**. This is a construction theorem, not a global optimality or exhaustive worldwide-priority claim.

| Dimension | Attributable external-source baseline | New bound | Source |
|---:|---:|---:|---|
| 93 | 17,179,869,184 | **42,951,379,608** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 94 | 34,359,738,368 | **64,426,347,212** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 95 | 68,719,476,736 | **107,376,020,548** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 96 | 137,438,953,472 | **206,160,530,240** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) |

These comparisons deliberately use external pre-existing work, not our previous constructions. The shortened baselines are derived from the published binary parent, not represented as independently published numerical records in those exact dimensions.

## Dense layer and95D improvement

The [preceding three-alphabet package](../dim93-96-three-crosspolytope-alphabets/) supplies the unchanged93D,94D and96D dense layers and the complete inherited Z8 parent minimum replay. In95D, replace its84 root heads by100 distinct E7 heads, assembled from three isometric38-head sections of the two-row prefix configuration. Every head has squared norm8 and satisfies x0+x4=0, so the head space has dimension7. Within each color its full38-by38 Gram matrix is exactly the source prefix Gram. Across colors, distinct E7 roots have inner product at most4, while the22-row body contributes at most44. Counts are100·8^10=107,374,182,400 dense points.

An orthonormal7D realization replaces(x0,x4) by sqrt2*x0 and retains the other six coordinates. Its exceptional coordinate can reach±2sqrt2. It is explicitly excluded from the RS supports rather than incorrectly treated as bounded by2.

## Compatible sparse layers

Use a length6, dimension4 Reed–Solomon evaluation code overGF16, with five full16-coordinate groups and one restricted group. Each support receives the32 even-parity sign words, with amplitude4. Distinct supports intersect in at most3; each point has squared norm96. All used dense coordinates have magnitude at most2, so dense/RS inner products are at most48. The independent checker enumerates every nonzero cubic and checks the field polynomial and root counts.

Add scaled roots sqrt48*(±e_i±e_j) with both coordinates in one RS group. In95D also include pairs from the excluded coordinate to every other coordinate. Each root support meets an RS support at most once. All root/root, dense/root and sparse/root inequalities are checked exactly, including the exceptional95D coordinate. The dense,6-support and2-support layers are disjoint by their nonzero-coordinate counts. No earlier sparse22 construction or other unverified layer is unioned.

## Complete evidence

Run `python verify.py` (Python3.10+, NumPy). It replays the preceding complete parent/three-alphabet proof, the original E7 head certificate, original RS/root augmentation verifier and a separately authored independent check. Frozen producer and independent packages are preserved unchanged. Their research-stage status remarks are historical. `claims.json` records only the four approved figures; newer candidate constructions are not included.

### Script inventory

- `final.py`
- `three_basis_final93_96_20261117/base/verify.py`
- `three_basis_final93_96_20261117/e7/search.py`
- `three_basis_final93_96_20261117/e7/verify.py`
- `three_basis_final93_96_20261117/verify_augmentation.py`
- `three_basis_final_independent_20261117/check.py`
- `three_basis_final_independent_20261117/input/base/verify.py`
- `three_basis_final_independent_20261117/input/e7/search.py`
- `three_basis_final_independent_20261117/input/e7/verify.py`
- `three_basis_final_independent_20261117/input/verify_augmentation.py`
- `verify.py`

The 94D and 95D rows are a retained construction, superseded by the [E6 and mixed E7 head exchanges](../dim94-95-e6-e7-head-exchange/). Frozen producer and independent proof files below retain their original historical values unchanged.
