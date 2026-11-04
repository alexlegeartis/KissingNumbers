# Classified free-zone equators and coupled 13D exchanges

**Current result: no new record.** Eight explicitly generated nonclassical signed equators give exactly verified 1154-point 13D seeds. Joint equator/cap replacement searches have not produced 1155. `manifest.json` separates known-result controls, numerical outcomes, execution status and pinned source/code hashes. Larger finite MILPs remain unresolved at their time limits; none of their solver bounds is a mathematical kissing-number bound.

## Sources and classification correction

Fetched on 4 November 2026, with byte counts and SHA-256 hashes in `sources.json`:

- Takhanov–Yun, [arXiv:2606.03299v1](https://arxiv.org/html/2606.03299v1), full primary HTML.
- The actual [ZE99 data file](https://github.com/tretoef-estrella/kissing-number-13-zinoviev-ericson-reproduction/blob/main/Ze99%201154%20data.md), carrying a transcription attributed to Cohn's coordinate archive. This mirror's attribution was not independently checked against the full original archive here. Its 1154 coordinate rows are independently verified mathematically, including all 48 irrational points; its saved verification claims are not trusted.

The classification must be read with its equivalence relation. [Section V, Theorem 4](https://arxiv.org/html/2606.03299v1#Thmtheorem4) partitions the 1579 signed-permutation classes into 1157 + 278 + 140 from three constructions, plus four classical lifts from egotistic support codes. Those four are **not all classical classes**. [Table I](https://arxiv.org/html/2606.03299v1#S2.T1) gives 17 optimal unsigned support-code classes for n=12. Classical full-sign lifting preserves and reflects support permutation equivalence: signed permutations merely permute supports, and every support permutation lifts. Thus the corresponding count is 17 classical and 1562 nonclassical **under signed permutations**. Construction 1 with both zones classical is an immediate counterexample to calling all first 1575 nonclassical. No classification of the resulting full 13D arrangements up to arbitrary orthogonal isometry is claimed here.

No complete 1579-entry catalog was downloaded. Actual seeds below come from the primary article's explicit [Section V Construction 1](https://arxiv.org/html/2606.03299v1#S5.p3) and [Section VI-A switch formulas H2–H5](https://arxiv.org/html/2606.03299v1#S6.SS1), not guessed representatives or a claimed reconstruction of the entire catalog.

## Actual construction pipeline

A round-robin 1-factorization of K6 supplies five colour classes, each containing three edges. Pair equally coloured edges in two disjoint six-coordinate blocks: 45 supports, with all 16 signs, give 720 points. Add a 48-point signed weight-four code in each six-coordinate free zone. The result has 816 signed weight-four points in 12 coordinates.

The six-zone codes are generated from the classical three-support code using the article's sign-dependent swaps. This wave uses zone indices 0 (classical), 2 (H2), 3 (H3), 4 (H4) and 5 (H5). It does not repeat the parent's already tested H1-only/random-signed-permutation-union experiment. The pairs are:

```
(0,2) (0,3) (0,4) (0,5) (2,2) (2,3) (3,4) (4,5)
```

Each equator is checked for 816 distinct points of squared norm 4 and maximum distinct product 2. Its support multiplicities include partial sign fibers, so it is nonclassical in the paper's signed-coordinate sense. `seed_profiles.json` records support profiles and antipodal-pair counts. These labelled seeds are not asserted mutually nonisometric.

Scale the equator by two and give it height zero. Retain the fetched sign and axis cap rows at heights ±2 and the two poles. Every seed then has 816 + 336 + 2 = 1154 distinct points of squared norm 16. Exact verification checks all **665281** pairs and proves full rank 13 by a nonzero algebraic minor modulo 13. The eight compressed NPZs are `seed_i_j.npz`; their exact receipts are in `results.json`. The fetched source control is separately recorded in `source_exact.json`.

`program_H2.json` is a real executable construction DAG:

1. Build the classical 48-point signed code in six coordinates.
2. Apply `signed_layers.conditional_swap` to produce H2.
3. Use `signed_layers.free_zone_join` to create the 816-point equator from the two zone codes and coloured edges.
4. Use `signed_layers.coupled_lift` to create the 1154-point quadratic-field configuration.

`program_H2_receipt.json` records independent terminal replay. This is a known-result construction control, not an improvement.

## Bounded searches across equator and caps

All searches use separate variables for individual oriented points. **No antipodality constraint is imposed.** An insertion requires removal of every old point it violates, including both equator and cap points. Candidate/candidate conflicts are also enforced. Objective: inserted points minus removed old points. An output of 1155 or more must pass the exact independent quadratic verifier before being reported as an improvement.

| Experiment | Full candidate grammar / seed | Admitted model or trade scope | Result |
|---|---:|---|---|
| `experiment.py` | 125874 integer directions: all signed supports 1–5 across all 13 coordinates, plus body supports 1–5 with last coordinate ±2 | Eight seeds; 408 new candidates with at most five old blockers per seed | Gain 0; finite numerical MILP optimum, not a rigorous negative bound |
| `expanded.py` | 175442 directions, adding all 13728 mixed 3/2/√3 triples, 15840 radical weight-four cap bodies, and 20000 sampled integer 3/2/1/1/1 directions | Three seeds; 7147 candidates pass an eight-blocker cutoff; cheapest 3200 enter each MILP | Gain 0; all three 20-second solves time out, with loose finite-model bounds |
| `small_trades.py` | Same expanded pool | Three seeds; 8888, 8842 and 8707 distinct proposed removal sets from single blocker sets and sampled pair unions, at most 12 removals | Gain 0; these particular sets unlock at most five candidates. This is not all possible removals |
| `unequal_heights.py` | 555552 fixed-norm16 points in five new unequal-height families | Full blocker distributions computed before sampling 300 cheapest per family for joint MILPs; no eight-blocker admission cutoff | Three 1500-candidate MILPs reach numerical optimum gain 0; no rigorous negative bound |

The expanded mixed radical candidates all had more than eight blockers, so **none entered those MILPs**. Their generation and screening must not be described as a successful search of their complete finite exchange graph. The unequal-height follow-up removes that shortcut and explicitly admits representatives from every family, even at higher blocker cost.

All three unequal-height MILPs reached numerical optimum gain 0 in about five to six seconds each. The minimum old-blocker costs were respectively 20, 23, 52, 37 and 30; full histograms are saved. All 1154 old points appeared as removable variables, including equator, caps and poles. These solver outcomes establish neither an exact finite exclusion nor a global upper bound.

The five unequal-height families are:

- Height ±1, body shape one ±2 and eleven ±1: 98304 points, body squared norm 15.
- Height ±1, body shape ±3, ±2, ±1, ±1: 190080 points, body squared norm 15.
- Height ±√3, body shape ±3, ±2: 1056 points, body squared norm 13.
- Height ±3, body shape ±2, ±1, ±1, ±1: 63360 points, body squared norm 7.
- Height ±3, body shape seven ±1: 202752 points, body squared norm 7.

These differ from the parent's already tested 3520 triple-support height-±2 caps and its 3072 rational-rotation cap pool. The new searches still sample finite model subsets and do not cover continuous deformations, all coordinate alphabets, all layer heights or arbitrary coupled replacements.

## Exact arithmetic and reusable operators

`discovery/engines/signed_layers.py` adds three operators, with no dimension-selection table:

- `conditional_swap`: an explicit signed-coordinate predicate and a metric-preserving coordinate swap. The source and resulting integer configuration are independently verified; the metadata alone does not establish injectivity.
- `free_zone_join`: two signed norm-four codes and paired edge colour classes. It constructs disjoint block codes and all signed bridges, then verifies every output pair.
- `coupled_lift`: an integer equator, supplied cap-body coordinates in Q(√3), integer heights, integer equator scale and pole amplitude. It does not infer compatibility from a layer formula: every output pair is checked.

The quadratic backend represents raw points as `A + B*sqrt(3)` with integer coefficient matrices. It requires positive integer squared norms, which may differ between points. For a positive raw dot `p+q*sqrt(3)`, normalized kissing compatibility is equivalent to

```
4*(p+q*sqrt(3))^2 <= norm_i * norm_j.
```

The sign comparison uses integer squares and sign cases; negative raw dots are automatically compatible. Normalized coordinates can require additional square roots, even though raw coefficients and all checked inequalities remain in Q(√3). No interval or floating-point certificate is claimed.

Modulo 13, √3 maps to 4. A nonzero minor under this ring homomorphism proves the corresponding exact algebraic minor nonzero. Full modular column rank therefore proves full real rank; a deficient modular rank is reported only as a lower bound, not an exact rank. The certificate establishes the explicit output geometry. Switching/lift provenance is separately represented by the construction DAG and source receipts.

Output dimension, point/pair/rank-work counts and coefficient-size bounds are checked before coupled output allocation and multiplication. Verification batches small Gram blocks and never builds a full group table. Data hashes identify inputs and never replace the mathematical replay.

## Prior-work scope and reproduction

`KNOWLEDGE.md` §66 decodes the ZE99 layers and describes antipodal pair swaps; §67 and `verifications/closed/dim09-19-record-maximality/` study record-coordinate maximality. The cited 12D discussion already knows the signed-Johnson classification. Those do not establish that all classified-equator 13D coupled exchanges have been searched. The old `research/dim13/` work referred to in §66 is not present in this checkout, so absence of that exact branch cannot be proved from this checkout. The parent explicitly supplied its fresh H1/random-union and cap-pool negative results; this wave avoids treating those as new work.

```sh
OPENBLAS_NUM_THREADS=1 python research/coupled13/experiment.py
OPENBLAS_NUM_THREADS=1 python research/coupled13/expanded.py
OPENBLAS_NUM_THREADS=1 python research/coupled13/small_trades.py
OPENBLAS_NUM_THREADS=1 python research/coupled13/unequal_heights.py
OPENBLAS_NUM_THREADS=1 python -m unittest discovery.tests.test_signed_layers -v
python -m discovery run research/coupled13/program_H2.json --out /tmp/coupled13-H2-receipt.json
```

The final command creates a new receipt and refuses to overwrite an existing destination; direct programmatic reproduction uses `discovery.runner.run_program`. The ten engine tests include high-precision sign-comparison oracles, rejection of duplicates and malformed coefficients, unequal norms, actual H2/free-zone output, a real 1154-point nonclassical seed and output resource preflight. No software installation, commit or push was performed. The remaining blocker is mathematical/search progress: **no exactly verified 1155-point configuration has been found**.
