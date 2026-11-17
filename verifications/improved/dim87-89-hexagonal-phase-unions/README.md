# Dense ternary fibers, nonlinear hexagonal phases and pure-head roots

Exact lower bounds: **7,266,268,830 in 87D**, **9,590,791,776 in 88D**, and **12,690,155,704 in 89D**. These are construction theorems, not global optimality or exhaustive worldwide-priority claims.

| Dimension | Attributable external-source baseline | New bound |
|---:|---:|---:|
| 87 | 6,218,178,164 | **7,266,268,830** |
| 88 | 6,218,179,920 | **9,590,791,776** |
| 89 | 6,218,181,330 | **12,690,155,704** |

These baselines retain the live table's external-source comparison: [Nebe's Gamma72](https://arxiv.org/abs/1008.2862v3) with the established lower-dimensional orthogonal-sum tails in [common/published.py](../../../common/published.py). They are derived lower bounds, not our earlier constructions or claims of separately tabulated records in these exact dimensions.

## Compatible union

The existing independently audited dense constructions have 42 unit-triangle body blocks, 9, 12 or 16 prefix fibers, and a head of squared norm 3 in dimensions 3, 4 or 5. Each fiber contains 3^18 words; both dense signs give 2·M·3^18 points of squared norm 45.

Replace, rather than retain alongside, all earlier sparse layers with the 37,468 explicitly checked weight-15 Pless words on the selected 42 coordinates and the exact 7,812-word length-15 ternary phase code. This contributes **292,700,016** points. The original author phase file and complete pair-distance replay are in the [phase package](../dim90-92-hexagonal-phase-layers/). A root/triangle product has absolute value at most 3/2, so 15 occupied sparse blocks give dense/sparse inner product at most 45/2. Sparse heads are zero; zero body blocks distinguish sparse points from dense points.

Add the scaled D3, D4 or D5 root system solely in the head space, contributing 12, 24 or 40 points of squared norm 45. Pure-head pairs have inner product at most 45/2. Dense/pure-head products are at most sqrt(135), strictly below 45/2, and sparse/pure-head products vanish. These points are distinct because every dense or sparse point has a nonzero body. No extra antipodal doubling is applied.

## Verification

Run `python verify.py` (Python 3.10+, NumPy). It replays all relevant original dense-head checks and complete nonlinear-phase checks, then the producer consequence and separately written integer/Fraction union audit. Original packages and source premises remain explicit: the published ternary parent minimum 15 is inherited, not newly proved by this union. Frozen producer and independent artifacts are preserved unchanged below.

### Script inventory

- `final.py`
- `hexagon_87_89_consequence_20261117/verify.py`
- `hexagon_87_89_independent_20261117/audit.py`
- `verify.py`
