# Composite cyclic-factor and extension campaign, 4 November 2026

No new kissing construction or record is claimed. This bounded campaign adds exact negative evidence for competitive cyclic modules beyond the earlier single half-degree-factor search. It writes only this directory and does not change the reviewed discovery engine or archive. It is research evidence awaiting any future audited adaptation into the shared construction grammar, not a registered operator or a dimension-partitioned engine.

The uniform grammar uses every monic divisor of `x^n-1` over F2, its unique coprime lift to Zq, the resulting free cyclic module, a possible extra coordinate, and the homogeneous Gray map to binary signs. Here q is 4, 8, or 16, n is odd with 3 <= n <= 47, and binary ambient length `(n+1)q/2` lies in 25–96. In particular this includes products of distinct irreducible factors, factors `x+1`, composite cyclic lengths, and ranks beyond the old half-degree family. The finite scan covers 32 (n,q) pairs and **832 nonzero modules**.

The comparison screen uses dated repository values from the existing baseline snapshot, benchmark catalog and README. Missing exact dimensions use the nearest smaller dimension's bound by isometric embedding (notably 48 uses 46 and 72 uses 71), making pruning conservative. Source files and hashes are recorded and checked. These are repository comparisons, not a fresh global-priority survey. A union of all repeated local binary masks has at most `2^(q/2)` translates; even this generous forecast prunes 315 modules. This cardinality screen does not prune arbitrary nonconstant coordinate masks, outer codes, geometric extra layers, or other assembly operators.

Of the **517 modules** retained by this forecast, the independently replayed primary classification is:

| Exact evidence | Modules | Scope |
|---|---:|---|
| A word of weight w with `w+q/2 < ceil(N/4)` | 417 | Every value of one extra ring coordinate fails for this full parent, including nonlinear extension rules |
| Explicit words reject every scalar-sum extension | 91 | All q rules `s * sum(coordinates)` fail; other linear functionals remain open |
| Six inconsistent linear-extension constraints | 2 | Every Z4-linear extra-coordinate functional fails |
| Sparse/residue screen remains unresolved | 7 | Four Z4 Golay-family modules, two Z8 Golay-family modules, and one small full Z16 module |

“Unresolved” means this rejection screen did not exclude the parent; it is not a new candidate or a certified construction. The six Golay-family survivors include known controls and were not subjected to a new distance proof here. The Z8 DGLS two-coset count **137438953472** remains below the full **142162362560** repository comparison. No gain to its extra geometric layers was found or claimed.

## Two useful obstructions

The rank-15 length-31 Z4 modules have parent cardinality `4^15 = 1073741824`, exceeding the recorded 64D comparison **335496064** if distance 16 could be achieved. For each of the two reciprocal modules, six explicit base words have Lee weight 14. An extra Z4 symbol must therefore equal 2 on each word to reach 16. Linearity imposes six equations `f(message_i)=2`. The supplied coefficients combine their left sides to zero modulo 4 and their right sides to 2. Thus `0=2 (mod 4)`, excluding all `4^15` linear functionals without enumerating them. This does not exclude nonlinear extension assignments. The rank-16 parents have weight-13 words, excluding every possible value of an added coordinate.

All ten initially surviving length-45 Z4 modules (six rank 18, four rank 17) contain an exact weight-16 word. Their binary length would be 92 and target distance 23, while an additional Z4 symbol contributes at most 2. Thus even arbitrary extra-coordinate values fail. Residue enumeration found these witnesses after the sparse-message scan missed them.

## Retracted unrestricted repeated-mask bound

**The unrestricted 262144-point upper bound and the statement that this branch was closed are retracted. This branch remains unresolved.** The cross-distance dynamic program computes minimum distance including coincident points. A minimum of zero proves overlap; it does not establish an inadmissible positive distance after deduplication. Consequently the supplied incompatible-subspace argument only bounds **pairwise-disjoint** translate unions, not general unions of their point sets.

The reviewer supplied a concrete counterexample to the claimed component incompatibility: for scalar 1 and repeated mask `0x0f0f0f0f`, each parent has 4096 words, their intersection has 2560 words, and their union has 5632 distinct words with minimum distinct distance 8 (15856896 pairs checked). This must be retained as a regression case for any future overlap-aware compatibility operator. See `overlap_counterexample.json`.

The recorded 262144 value is now explicitly scoped to pairwise-disjoint translate unions only. It does not close the unrestricted repeated-mask branch. The primary count of seven unresolved parent-screen controls is unchanged. The earlier 417 arbitrary-coordinate, 91 scalar-sum, and two linear-extension obstructions do not use this flawed overlap argument.

## Reproduce

From repository root, using installed software only:

```sh
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/rings/composite_screen.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/rings/refine_residues.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/rings/linear31_obstruction.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/rings/repeated_mask_dp.py
python research/discovery_campaign_20261104/rings/replay.py
python research/discovery_campaign_20261104/rings/replay_repeated_masks.py
python research/discovery_campaign_20261104/rings/overlap_regression.py
```

The four searches took approximately 2, 15, 4, and 1.4 seconds. The two independent standard-library replays took approximately 0.2 and 2.1 seconds. No installation, push, publication, giant point allocation, or external messaging was performed. The first three receipts supply generators, witnesses and bounded search provenance; the dynamic-program receipt supplies full cross-distance tables and incompatible subspaces. Each search records a source hash, and refinement pins its input receipt. Hashes detect changes and do not authenticate mathematical truth.

Both replays recompute the mathematical facts and raise explicit exceptions. They do not import NumPy, SymPy, or search modules and remain active under `python -O`. The main replay checks the exact 32-family domain and unique coverage, binary irreducibility, factor products, ring lift identities, module ranks, recorded reference values, every supplied word and homogeneous weight, scalar exclusions, and elementary six-word contradictions. The second reconstructs the entire 16-by-256 distance table and verifies the subspace coset bounds for pairwise-disjoint translate unions only. The main PASS receipt is `independent_replay.json`; `repeated_mask_independent_replay.json` explicitly reports `PASS_FOR_DISJOINT_TRANSLATES_ONLY`. The overlapping-union regression is independently replayed in `overlap_counterexample.json`.

Remaining research work is to leave the excluded full-parent/extension families: new residue/torsion modules, selective subcodes, nonlinear switching, nonconstant masks, outer-code assembly, or extra geometric layers. No task permission or environment blocker remains.
