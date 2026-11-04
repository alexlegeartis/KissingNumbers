# Exact construction coverage, dimensions 1–24

The executable audit has exactly24 rows. Current exact target recipes cover
**1–10,13,14,16–18,22–24**. The remaining dimensions **11,12,15,19–21** have
explicit gaps. These are existing lower bounds, with no new record claim.
The9D route executes runtime support search; its separate held-out benchmark
is documented in `discovery/docs/CONSTRUCTION_SEARCH.md`. The12D search reaches
historical840, below current841. Large17D/18D and Leech controls use explicitly
declared resource profiles; their self-contained receipts are compressed.

The source is [Cohn's table](https://cohn.mit.edu/kissing-numbers/), successfully retrieved on **2026-11-04T12:36:52.845533+00:00**. The preserved HTML is `research/discovery_campaign_20261104/baselines/cohn.html`, SHA-256 `94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e`. Earlier browser/network failures in the incoming audit are historical observations; this later successful response supplies every target. Within this range, only dimensions 1, 2, 3, 4, 8 and 24 have proved exact kissing numbers; the other targets are lower bounds.

The [campaign matrix](../../research/discovery_campaign_20261104/coverage_1_24/README.md) gives all 24 targets, per-row source citations, recipes, proof receipts, provenance, resource limits and gaps. `discovery/benchmarks/coverage_1_24_20261104.json` is an executable-checked export of that matrix. `root_controls_20261104.json` contains the eight static root fixtures; `coverage_controls_20261104.json` contains the other replayed routes. The default CLI comparison catalogue is deliberately unchanged and remains partial, so successful exact verification need not imply that default `compare` has a baseline for that dimension.

The native `lattice.root_system` operator constructs A, D and E root systems at runtime. Recipes `roots01.json` through `roots08.json` reach the sourced targets 2, 6, 12, 24, 40, 72, 126 and 240. The incoming `root_controls.py` benchmark generator supplies a second route: it constructs deterministic formula coordinates and writes ordinary `geometry.integer_code` recipes (`roots_a1`, `roots_a2`, `roots_d3`–`roots_d5`, `roots_e6`–`roots_e8`). The coverage checker requires exact equality between each static fixture and regenerated recipe. Executor ancestry is still recorded as imported coordinate data, with verified formula generation recorded separately. This prevents a supplied certificate from acquiring an unsupported discovery claim.

The root generator uses A1 antipodes, A2 differences of three coordinate vectors and D roots with two nonzero signed entries. Doubled E8 roots comprise 112 signed vectors with two entries ±2 and 128 vectors with eight entries ±1 and even minus parity. The E7 and E6 controls select the sections x6+x7=0 and additionally x5+x7=0 (zero-based indices). The independent verifier checks distinctness, common norm, every pair inequality and exact rank. Ambient coordinates are reported separately: the static A2 fixture uses three coordinates and E6/E7 use eight. Native weighted root recipes can use minimal rank coordinates.

An explicit research-recipe allow-list also includes `research/coupled13/program_H2.json`, which executes registered switching, free-zone joining and coupled lifting to reproduce known13D1154. It retains imported ancestry from supplied classical equator and cap coefficients; independent exact replay checks all665281 pairs in Q(sqrt(3)) and full rank13. This is a known construction, not a new record. Planner configurations lacking runner nodes and the UNVERIFIED96D coset forecast receive no exact coverage.

The 5D Gram control is a supplied exact 40-point Gram certificate. The E8 section chain generates 126 points of rank 7 and 72 of rank 6 while retaining eight coordinates. The 10D and 14D recipes reconstruct known exact 510- and 1932-point constructions. Separate numerical 10D subgroup-character search is not promoted to exact autonomous synthesis. The 16D `groups.real_stabilizer_sector` recipe generates an implicit 4320-point affine-flat/quadratic-sign family and replays a symbolic count and separation proof; it materializes zero spherical points. The 256-point octacode recipe remains a partial 16D control.

Existing repository prototypes receive no credit until an actual registered
recipe and exact proof run exists. Historical19D11948 and21D29768 are below the
captured targets12460 and31560. The18D8358 package now has a native weighted
quadratic adapter. The known17D5730 construction is independently reconstructed
from Boolean quadratic/RM formulas. A24D exact orbit recipe generates196560
points from four algebraic actions; closure, transitivity and196559 seed
comparisons certify all19317818520 unordered pairs. Exact source subsets then
produce23D93150 and22D49896. The same signed-pair-normal search gives only
21D27720 and20D17400, below31560 and19448. These shortfalls remain explicit.

The external author's original 14-file patch is preserved byte-for-byte at `research/discovery_campaign_20261104/coverage_1_24/integration/coverage-benchmarks-20261104.patch`, SHA-256 `7ac4796cc14c476424e3a09888e499c077844828c9921b7508c8e81e384bd741`. Its source commit `0a6152a9ea6caf8b7865e2f41c2097a838566e2a`, author Alex Kravatsky, date and reviewed base are recorded in `integration/provenance.json`. That patch preserves original construction-attribution notes and external measurements. Active files replace its stale source availability, native-operator and 16D coverage statements with current independently replayed evidence.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --write
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/reconcile_views.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --check
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/reconcile_views.py --check
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/test_matrix.py
OPENBLAS_NUM_THREADS=1 python -m unittest discovery.tests.test_root_controls
```

The matrix checker rejects changed source bytes, wrong dimensions, omitted eligible outputs, modified recipes/receipts, false import/generation labels and fixtures that do not regenerate. The view checker requires exact equality with the checked canonical matrix and an unchanged comparison catalogue. No network requests, installation, record writes, commit or push are part of these commands.
