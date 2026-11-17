# Two-alphabet fibers with compatible sparse layers in 92 dimensions

**Exact lower bound: 26,038,868,280.** This is a direct 92-dimensional construction, not a higher-dimensional embedding with extra poles. No global optimality or exhaustive worldwide-priority claim is made.

| Dimension | Attributable external-source baseline | New bound |
|---:|---:|---:|
| 92 | 8,589,934,592 | **26,038,868,280** |

The comparison is obtained by shortening the published nonlinear binary parent of [Duursma, Greferath, Litsyn and Schmidt (2001)](https://doi.org/10.1109/18.923742). It is a derived external-source construction, not our previous result or an independently tabulated 92D record.

## Four compatible layers

| Layer | Points | Method |
|---|---:|---|
| Dense | 25,769,803,776 | Two Gray crosspolytope alphabets on 22 body rows, with 24 distinct signed D4 heads and 8^10 words per fixed two-symbol fiber |
| Sparse22 | 268,435,456 | Independently enumerated shortened quaternary [22,8,11] support code and binary [22,12,6] sign code |
| Reed–Solomon | 628,864 | Cubic evaluation over F17, six-coordinate supports, 32 even-parity signs, one head coordinate and five body coordinates |
| Coordinate poles | 184 | Both signs on all 92 coordinate axes |
| **Total** | **26,038,868,280** | Exact disjoint union |

Every point has squared norm 96 and every pair has inner product at most 48. The complete independent checker reconstructs the head system, enumerates original code spans before shortening, compares the resulting finite sets, checks all evaluation minors and support counts, and verifies exact cross-layer inequalities. Sparse and dense layers are distinguished by support patterns; no legacy layer is silently included.

## Complete verification and source scope

Run `python verify.py` (Python 3.10+, NumPy). The wrapper first replays the [complete parent verification](../dim93-96-three-crosspolytope-alphabets/), including the inherited Z8 minimum calculation. It then runs every original layer verifier and the separately authored 92D audit. The original producer archive, checksums, independent audit and parent inputs are included unchanged under `research/`, preserving their original relative paths and historical remarks. Only the approved 92D result is registered; the original package's weaker 93–95D variants are retained evidence, not promoted claims.

### Script inventory

- `final.py`
- `research/ring2_92_independent_20261117/audit.py`
- `research/ring2_final92_95_20261117/verify_dense.py`
- `research/ring2_final92_95_20261117/verify_rs17.py`
- `research/ring2_final92_95_20261117/verify_sparse.py`
- `research/three_basis_independent_20261117/check.py`
- `research/three_basis_independent_20261117/input/parent_independent.py`
- `research/three_basis_independent_20261117/input/verify.py`
- `research/three_basis_independent_20261117/replay_parent.py`
- `verify.py`
