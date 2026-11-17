# Four-row partition heads in 92 dimensions

**Exact lower bound: 29,797,911,144.** This direct construction replaces the previous 92D dense and sparse layers. No global optimality or exhaustive worldwide-priority claim is made.

| Dimension | Attributable external-source baseline | New bound |
|---:|---:|---:|
| 92 | 8,589,934,592 | **29,797,911,144** |

The comparison is obtained by shortening the published nonlinear binary parent of [Duursma, Greferath, Litsyn and Schmidt (2001)](https://doi.org/10.1109/18.923742). It is a derived external-source construction, not our previous result or an independently tabulated 92D record.

## Dense construction: 29,796,335,616 points

Retain 20 four-dimensional blocks from the length-24 Z8 code. Replace the first four blocks by signed weight-four heads on 12 axes. Three explicit partitions of these axes into four groups assign each head axis to one deleted ring symbol. The selected transversal-support counts are 81, 24 and 6. Every support has all 16 signs, producing **1,776 heads**, each of squared norm 16.

Within each color, the head Gram matrix exactly equals that of the four removed blocks. Between colors, supports intersect in at most two axes, so head inner products are at most 8. The retained 80-dimensional bodies use three isometric crosspolytope alphabets, with cross-body inner products at most 20·2 = 40. Thus every full point has squared norm 96 and every pair has inner product at most 48.

The first-four-symbol projection is surjective, with exactly 8^8 codewords per prefix. Distinct heads separate the fibers and deleting a fixed prefix is injective. Therefore the dense count is 1,776·8^8 = **29,796,335,616**.

## Compatible augmentation

Add **1,572,864** six-support Reed–Solomon points and **2,664** within-group roots. The RS layer uses cubic evaluations over GF16, one 12-coordinate group, five 16-coordinate groups, and all 32 even-parity signs with amplitude 4. Distinct supports intersect in at most three coordinates. All dense coordinates have magnitude at most 2, so the dense/RS bound is exactly 6·4·2 = 48.

The scaled roots use pairs within one coordinate group, so each root meets an RS support in at most one coordinate. Exact root, RS and dense cross bounds are independently checked. Dense points have at least 24 nonzero coordinates, whereas the new layers have six and two, proving disjointness. The old sparse22 and other legacy layers are **not** included.

## Complete verification

Run `python verify.py` (Python 3.10+, NumPy). It replays the complete inherited Z8 parent verification, then the original producer check and separately authored independent audit. The latter checks all frozen hashes, all 4,096 first-four inputs, every head/Gram block, finite-field arithmetic, all RS supports, sign distances and cross-layer geometry. Original producer, independent and parent packages are preserved unchanged under `research/`, with their archive hashes and original relative paths.

### Script inventory

- `final.py`
- `research/three_basis_independent_20261117/check.py`
- `research/three_basis_independent_20261117/input/parent_independent.py`
- `research/three_basis_independent_20261117/input/verify.py`
- `research/three_basis_independent_20261117/replay_parent.py`
- `research/z8_fourrow92_independent_20261117/audit.py`
- `research/z8_fourrow92_partitions_20261117/search_stronger.py`
- `research/z8_fourrow92_partitions_20261117/verify.py`
- `verify.py`
