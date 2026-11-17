# E6 and mixed E7 head exchanges in dimensions 94 and 95

Exact construction lower bounds: **66,573,830,860 in 94D** and **109,523,504,196 in 95D**. No global optimality or exhaustive worldwide-priority claim is made.

| Dimension | Attributable external-source baseline | New bound |
|---:|---:|---:|
| 94 | 34,359,738,368 | **66,573,830,860** |
| 95 | 68,719,476,736 | **109,523,504,196** |

Both comparison baselines are derived by shortening the published nonlinear binary parent of [Duursma, Greferath, Litsyn and Schmidt (2001)](https://doi.org/10.1109/18.923742). They are external-source constructions, not our earlier results or independently tabulated exact-dimensional records.

## Construction

The unchanged length-24 rank-12 Z8 parent and three crosspolytope alphabets are replayed through the [preceding complete package](../dim93-96-e7-rs-root-augmentation/). Fixing two ring symbols gives a fiber of 8^10 words and an 88-dimensional body. Exact head Gram preservation retains the parent distance within each color; across colors, the body inner product is at most 44 and distinct head inner products are at most 4. Every full point has squared norm 96.

In 94D, three E7 sections restricted to an E6 hyperplane supply **62 distinct heads**, with color counts 24, 20 and 18. The explicit orthonormal six-dimensional realization has all coordinates bounded in magnitude by 2. In 95D, mixed E7 sections supply **102 distinct heads**, with color counts 38, 32 and 32. The seven-dimensional realization has one exceptional coordinate, excluded from the sparse supports exactly as in the preceding package.

The unchanged, independently checked Reed–Solomon and root layers add 1,837,772 points in 94D and 1,838,148 in 95D. Thus the totals are 62·8^10 + 1,837,772 and 102·8^10 + 1,838,148. All selected head pairs, coordinate isometries, source hashes and augmentation arithmetic are independently checked. Discovery searches and numerical optimization are not proof premises.

## Replay and evidence

Run `python verify.py` (Python 3.10+, NumPy). This replays the complete inherited parent and augmentation verification before both original new-head verifiers and the independent scalar checker. The complete independent archive and its original producer packages are preserved unchanged. Research-stage status remarks in those frozen files are historical. Only the two approved figures are registered here.

### Script inventory

- `e6e7_9495_independent_20261117/check.py`
- `e6e7_9495_independent_20261117/e7_mixed_sections95_20261117/search.py`
- `e6e7_9495_independent_20261117/e7_mixed_sections95_20261117/verify.py`
- `e6e7_9495_independent_20261117/e7_to_e6_94_20261117/search.py`
- `e6e7_9495_independent_20261117/e7_to_e6_94_20261117/verify.py`
- `final.py`
- `verify.py`
