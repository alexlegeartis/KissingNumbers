# Three crosspolytope alphabets with exact root-head compression

This package proves **42,949,672,960 in93D**, **64,424,509,440 in94D**, **90,194,313,216 in95D**, and **206,158,430,208 in96D**. It uses a published Z8 parent code in a new explicitly verified three-alphabet union. No claim of a newly discovered parent code, unrestricted optimum or exhaustive worldwide priority is made.

| Dimension | Attributable external-source baseline | New bound | Source |
|---:|---:|---:|---|
| 93 | 17,179,869,184 | **42,949,672,960** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 94 | 34,359,738,368 | **64,424,509,440** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 95 | 68,719,476,736 | **90,194,313,216** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) (shortened, derived) |
| 96 | 137,438,953,472 | **206,158,430,208** | [Duursma–Greferath–Litsyn–Schmidt (2001)](https://doi.org/10.1109/18.923742) |

Comparisons are to attributable pre-existing external work, not our own previous constructions. The shortening baselines are explicitly marked as derived; the original paper supplies the binary(96,2^37,24) code, not each lower-dimensional numerical claim.

## Mechanism

The eight even Gray-sign rows, their odd-parity counterparts, and the eight coordinate vectors ±2e_i are three isometric four-dimensional crosspolytopes. Each row has squared norm4, and every cross-alphabet inner product is±2. Apply the same free rank12 Z8 code of length24 to each alphabet. Within an alphabet, Gram preservation retains minimum squared distance96. Between alphabets, inner product is at most24·2=48. The three copies are disjoint by row parity or zero-coordinate pattern. They therefore give3·8^12=206,158,430,208 points of squared norm96.

For93–95D retain22 rows in88 dimensions and replace the two discarded rows by heads2(±e_i±e_j) in dimensions5,6,7. Three binary cuts partition every edge of the complete head-axis graph. In each color, the row-axis lists are disjoint and have size at most4, so every within-color prefix Gram is preserved exactly. Different colors use different crosspolytope alphabets: their body inner product is at most44, and their distinct root heads add at most4. Every fixed two-symbol prefix has exactly8^10 source words, verified by the explicit invertible ring-coordinate minor. There are40,60,84 heads respectively.

No previous sparse layer is automatically unioned with these constructions. In particular, the third alphabet has entries0 and±2, so compatibility arguments assuming all dense coordinates are±1 must not be reused.

## Source and full evidence

The parent is attributed to [Duursma, Greferath, Litsyn and Schmidt, IEEE Transactions on Information Theory47 (2001),1596–1598](https://doi.org/10.1109/18.923742). Its exact generator bytes and original independent minimum-proof implementation are supplied. The current wrapper replays that inherited parent proof, including its finite Z4 enumeration, and compares the complete result with the pinned source receipt. This is reuse of the established independent source proof, not a new independently written parent solver.

Run `python verify.py` (Python3.10+, NumPy). The complete original producer and independent audit packages are preserved unchanged. Research-stage phrases such as “not published” or “audit requested” describe their historical state. `claims.json` lists the approved current claims, and `final.py` supplies them to the repository registry.

### Script inventory

- `final.py`
- `three_basis_independent_20261117/check.py`
- `three_basis_independent_20261117/input/parent_independent.py`
- `three_basis_independent_20261117/input/verify.py`
- `three_basis_independent_20261117/replay_parent.py`
- `three_basis_z8_93_96_20261117/verify.py`
- `verify.py`
