# Bounded follow-up and certification review, 4 November 2026

**Outcome: no new kissing-number record or new binary-code claim.** The new experiment excludes a precise 36D construction family, with compact exact witnesses and an independent standard-library replay. It does not establish a kissing-number upper bound. The earlier numerical 28D output remains numerical; this review does not promote it to a certified lower bound.

## What was tried

The existing campaign searched Fourier/quadratic-phase representations, reconstructed the 10D and 14D controls, and synthesized the published DGLS code using factorization, Hensel lifting, a Gray map and compatible translates. Its tested ring lengths were 7 and 23 and its extension coordinate was always minus the coordinate sum. Merely increasing the Fourier sampling budget is not justified by its noncompetitive higher-dimensional outputs.

The follow-up instead tested the two degree-eight binary factors of `1+x+...+x^16`, their unique coprime monic lifts to Z4, and a free rank-nine cyclic module of length 17. Generalized Gray embedding after one additional Z4 coordinate has binary length 36 and cardinality `4^9 = 262144`. To yield a kissing sign layer, its distance must be at least `ceil(36/4) = 9`.

1. Enumerate every message for each factor and every scalar-sum extension `s * sum(coordinates)`, `s = 0,1,2,3`: eight construction programs, each with 262144 words. The minimum distance is 7 at scalar zero and 8 for every nonzero scalar. None qualifies as a kissing sign layer in R^36. Do not count these 262144-word binary codes as valid 36D kissing constructions.
2. Generalize the extra coordinate to **every Z4-linear functional on the nine message coefficients**. This is the full set of 262144 linear one-coordinate extensions for each fixed module. A base word of homogeneous weight 7 needs extension symbol 2; a base word of weight 8 needs a nonzero extension symbol. Only 11 such word constraints for factor 0, and 10 for factor 1, eliminate every extension.
3. Replay the finite obstruction independently using only Python integers, direct dot products, the homogeneous-weight tuple `(0,1,2,1)`, and `itertools.product`. Recompute every witness word's weight and examine every extension. Both modules have zero survivors. The replay imports neither NumPy nor the construction/search module.

The scalar search also screened cyclic length 11, but it has no individual irreducible half-degree factor over F2 in this grammar. Consequently there was no length-11 module run, and no 48D output. An empty grammar branch is not mathematical nonexistence evidence.

The exclusions apply to the **entire free module** with one **linear** extension coordinate. They do not exclude subcodes, nonlinear extension rules, puncturing, switching, other base modules, longer extensions, or other 36D configurations. The registered repository comparison in 36D is 501687, which already exceeds the entire 262144-word parent cardinality. Thus even a successful parent extension would be a construction-control result rather than an improvement; a competitive direction also needs assembly or a larger parent.

## Deliverables and reproduction

- `discovery/benchmarks/catalog.json`: versioned controls and dated comparison snapshots, including 26D, 28D, 32D, 36D and 96D. Each comparison gives its source, evidence and scope; it is not a global-priority oracle.
- `discovery/benchmarks/followup_extensions_20261104.json`: all eight generated matrices, lift histories, complete homogeneous-weight histograms and measured distances; the search source is embedded for audit.
- `discovery/benchmarks/all_linear_extensions_20261104.json`: both fixed generators, exact obstruction words, elimination counts, search/replay sources and parent-module hash.
- `discovery/benchmarks/replay_extensions.py`: portable independent finite-family replay.

From the repository root:

```sh
python discovery/benchmarks/replay_extensions.py
```

Expected: PASS, 262144 extensions checked per factor, 11 and 10 constraints, zero survivors. This replay proves the exclusion directly from the supplied module matrices and witness words; trusting the saved `status` is unnecessary. It does not need to re-factor the polynomial to prove this narrower claim.

The original experiments ran without installing software, enabling Actions or pushing changes. Their arithmetic is exact integer arithmetic. NumPy batches were small enough that products cannot overflow int64: the largest module multiplication sums nine products of integers in 0..3. The independent replay uses unbounded Python integers.

## Comparison and source review

Ganzhinov's primary [Highly symmetric lines](https://arxiv.org/abs/2207.08266) was re-read on 2026-11-04 and confirms the known 510-point 10D and 1932-point 14D controls. The campaign's exact coordinate receipts check distinctness, norms, every pair and exact rank. The 10D certificate has 18 weighted coordinate columns but rank 10; ambient coordinate count must not replace its dimension proof.

The DGLS authors' [university publication record](https://experts.illinois.edu/en/publications/a-zsub8sub-linear-lift-of-the-binary-golay-code-and-a-nonlinear-b/) was rechecked on 2026-11-04. It confirms the published `(96,2^37,24)` code and coset union, [DOI 10.1109/18.923742](https://doi.org/10.1109/18.923742). Reconstructing its ring lift from binary factors is a useful recovery control; it does not make the code a new discovery. The repository's two independent exact parent replays provide the distance evidence independently of abstract metadata. The stronger repository 96D comparison is **142162362560**, including its additional layers; comparing only with the old linear `[96,33,24]` code would be misleading.

Fresh retrieval of the [Cohn table](https://cohn.mit.edu/kissing-numbers/) failed in this follow-up: the web tool returned HTTP 502 and standard-library HTTP returned 503. This is source availability, not an approval denial. The catalog retains the already captured 2026-11-04 table and its recorded source hash in `research/group_search_20261104/baselines.json`. It records **26D 201136, 28D 221988, 32D 352800** as dated comparisons, without claiming a new complete literature survey. The repository source review `common/source_review_20261029_27_31.json` additionally supplies pinned primary-release evidence for the 27–31D comparisons.

The **36D 501687** value is supported by the current `completion/verification.json` and `completion/coordinate_verification.json` in the Coxeter–Todd package. The older `phase/verification.json` says 501633 and must not be substituted for the current completion. The catalog distinguishes repository certificates, published controls and captured external tables.

## Review of the unified architecture

The useful interface is a typed construction recipe with replayable evidence, not a collection of saved success flags. Core serialization and verification dispatch already establish that boundary. Certification should always recompute facts from construction data and compare claimed cardinality, dimension and distance. A stored `verification: exact` field is a receipt, not authority.

The following distinctions are necessary when extending the engines:

- **Metric versus geometry:** a ring module needs injectivity/cardinality and distance, followed by an all-symbol-pair Gray isometry check. A binary sign embedding then transfers `distance >= N/4` to the spherical inequality. Parent distance and union distance differ; each translate union needs cross-distance and disjointness evidence.
- **Exact dimension versus sufficient ambient dimension:** integer coordinate certificates can compute exact rank; implicit sign codes generally certify a construction in the stated ambient dimension without proving full span. Positive diagonal coordinate weights preserve rank and must enter every norm and inner product.
- **Numerical versus exact group claims:** rounded orbit closure, eigenspaces, phases, unitary tests and MILP compatibility graphs remain numerical. They do not prove the group finite or certify a spherical code. Independent exact coordinate reconstruction or verified interval/algebraic arithmetic is needed before accepting a lower bound.
- **Group compression:** representative-row cross-angle checks require genuinely transitive orbits under the same action and preservation of actual vector phases. Projectively canonicalized lines cannot be assigned arbitrary independent phases. A finite graph optimizer's bound applies only to that graph and tolerance, not to kissing numbers or all seed choices.
- **Scalability:** budgets should guard allocations before materialization. Large codewords and orbit unions should retain implicit algebraic descriptions; a `2^37` cardinality is not an instruction to create that many point rows. Proof dispatch must verify every compressed child, not accept untrusted child facts.
- **Provenance:** preserve source recipe, generated choices, seed/RNG, bounded grammar, rejected branches, numerical tolerances, search cutoff and independent checker identity. Known-result reconstruction, restricted rediscovery, numerical candidate and newly certified improvement should be separate labels.

The two implemented branches share artifact/proof dispatch, but this alone does not implement arbitrary group/code composition. In particular, double cosets, subgroup-class enumeration, higher-multiplicity character spaces, continuous relative phases, finite/Galois-ring generalization, partial signed weight-enumerator recurrences, nonlinear switching, general concatenation and geometric sections are future operators. Supporting one DGLS reconstruction does not establish coverage of every construction of similar sophistication. The finite extension obstruction supplied here is research evidence, not yet a registered construction operator.

## Next bounded experiments

1. Add a linear-functional extension operator with explicit low-weight-word obstruction certificates. The completed search makes a useful rejection control; a subcode/expurgation operator could evade it but should predict a competitive final cardinality before spending a larger budget.
2. In the code branch, test nonconstant translate masks and structured outer codes against exact local distance matrices, then attempt a distance certificate stronger than the conservative rowwise bound. A failure of the bound need not reject a valid actual union. Compare the final geometric objective, not just isolated parent distance.
3. For genuinely competitive high-dimensional candidates, prioritize new residue/torsion modules, concatenation or switching with compressed certificates and a cardinality forecast exceeding the strongest recorded comparison. The unique lift of an unchanged coprime factor pair is not an optimization parameter.
4. In the representation branch, add a small supplied catalogue of low-degree exceptional/sporadic representations and verified subgroup-character projectors, rather than re-running the same sampled Fourier family. Preserve numerical status until a separate exact/interval checker succeeds.

No pending authorization or environment blocker remains. New SOTA remains an open research objective; fresh comprehensive source comparison and an exact competitive candidate are both absent.
