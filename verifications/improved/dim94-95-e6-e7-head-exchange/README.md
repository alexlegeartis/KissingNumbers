# E6 and mixed E7 head exchanges in dimensions 94 and 95

Exact construction lower bounds: **73,016,281,804 in 94D** and **109,523,504,196 in 95D**. No global optimality or exhaustive worldwide-priority claim is made.

| Dimension | Attributable external-source baseline | New bound |
|---:|---:|---:|
| 94 | 34,359,738,368 | **73,016,281,804** |
| 95 | 68,719,476,736 | **109,523,504,196** |

Both comparison baselines are derived by shortening the published nonlinear binary parent of [Duursma, Greferath, Litsyn and Schmidt (2001)](https://doi.org/10.1109/18.923742). They are external-source constructions, not our earlier results or independently tabulated exact-dimensional records.

## Construction

The unchanged length-24 rank-12 Z8 parent and three crosspolytope alphabets are replayed through the [preceding complete package](../dim93-96-e7-rs-root-augmentation/). Fixing two ring symbols gives a fiber of 8^10 words and an 88-dimensional body. Exact head Gram domination retains the parent distance within each color; across colors, the body inner product is at most 44 and distinct head inner products are at most 4. Every full point has squared norm 96.

In 94D, a Gram-dominated assignment in the E6 hyperplane supplies **68 distinct heads**, with color counts 24, 22 and 22. Within each color the new head Gram is no greater than the original prefix Gram, which suffices to preserve every distance. The frozen 62-head isometric predecessor remains below. The explicit orthonormal six-dimensional realization has all coordinates bounded in magnitude by 2. In 95D, mixed E7 sections supply **102 distinct heads**, with color counts 38, 32 and 32. The seven-dimensional realization has one exceptional coordinate, excluded from the sparse supports exactly as in the preceding package.

The unchanged, independently checked Reed–Solomon and root layers add 1,837,772 points in 94D and 1,838,148 in 95D. Thus the current totals are 68·8^10 + 1,837,772 and 102·8^10 + 1,838,148. The 94D gain over the preceding 62-head certificate is 6,442,450,944 points (9.67715%). All selected head pairs, coordinate isometries, source hashes and augmentation arithmetic are independently checked. Discovery searches and numerical optimization are not proof premises.

## Replay and evidence

Run `python verify.py` (Python 3.10+, NumPy). This replays the complete inherited parent and augmentation verification before both original new-head verifiers and the independent scalar checker. The complete independent archive and its original producer packages are preserved unchanged. Research-stage status remarks in those frozen files are historical. Only the two approved figures are registered here.

### Script inventory

- `e6e7_9495_independent_20261117/check.py`
- `e6e7_9495_independent_20261117/e7_mixed_sections95_20261117/search.py`
- `e6e7_9495_independent_20261117/e7_mixed_sections95_20261117/verify.py`
- `e6e7_9495_independent_20261117/e7_to_e6_94_20261117/search.py`
- `e6e7_9495_independent_20261117/e7_to_e6_94_20261117/verify.py`
- `gram_dominated68/check.py`
- `final.py`
- `verify.py`

## Gram-dominated 94D continuation

`gram_dominated68/candidate.json` contains the exact 68-head assignment. `gram_dominated68/check.py` independently reconstructs alphabets and checks all 2,278 head pairs, coordinates, finite-field supports and augmentation counts. `gram_dominated68/AUDIT.md` records the full construction and limits; `gram_dominated68/inherited_replay.log` is a fresh full parent replay. `verify.py` repeats that inherited replay, the retained original checks and the new check. These data supersede only the 94D count, not the frozen historical input packages. The external-source baseline remains attributed shortening, not an exhaustive worldwide-priority conclusion.
