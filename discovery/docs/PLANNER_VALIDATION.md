# Planner validation and discovery deficits

This is the imported validation report from the planner patch's earlier base.
Its 81-test count and source fingerprints describe that workspace, not the
current integrated engine. Patch provenance and current type-integration checks
are recorded separately in
[planner_integration_20261104.json](../benchmarks/planner_integration_20261104.json).
The measurements below remain historical controls and failed quality gates;
integration tests do not turn them into competitive search successes.

Measured 2026-11-04 in the existing Python 3.12 environment, with
`OPENBLAS_NUM_THREADS=1`. No dependencies were installed. The machine-readable
[measurement report](../benchmarks/planner_synthesis_20261104.json) preserves
source fingerprints, catalogue digests, invocation limits, execution budgets,
ledger summaries and the two group-search output summaries.

The complete discovery suite passed **81 tests in 38.619 seconds**, including
16 planner tests. The repository audit (`python audit.py --repository-only`)
also passed after refreshing the line-ending file census for the eight new
files. Planner tests cover independent tiny
DAG enumeration, parameter combinatorics, two interruption points, exact resume
ordering, concurrent-writer rejection, changed-source/budget rejection,
resource-stop continuation, retry history, adaptive quality ordering, real
workers, brute-force binary-code checks and ordinary receipt verification.

## Generated coset subsets, then spherical embedding

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery plan discovery/recipes/planner_cosets8.json \
  --state-dir /tmp/planner-cosets8-complete --max-steps 25000 \
  --max-evaluations 200 --time-limit-seconds 120
```

The planner generated `cosets.build -> cosets.union -> geometry.sign_embed`
programs from types and available operations, without supplied graph slots or
input edges. It generated each zero-containing subset of the eight even
syndromes of the extended binary Hamming code. This uses the actual subset
domain and quotient-distance operators; it is not a list of prebuilt spherical
configurations or a replay of the E8 construction.

| Measurement | Result |
|---|---:|
| Generated stored prefixes/programs | 258 |
| Intermediate, non-target prefixes | 130 |
| Evaluated spherical candidates | 128 |
| Completed independent exact checks | 128 |
| Non-target leaves pruned at the node bound | 16,769 |
| Summed worker time, including interpreter startup | 43.731 s |
| Best exact cardinality in ambient dimension 8 | 128 |
| Known E8 target | 240 |
| Exact deficit | 112 |
| Fraction of target | 53.333% |
| `target_quality_passed` | false |

The retained [best receipt](../benchmarks/receipts/planner_cosets8_best.json)
re-verifies with `python -m discovery verify`. Its binary code has length 8,
128 words and minimum distance 2. The spherical object stays implicit. As with
the existing sign embedding, the dimension statement is ambient, without an
exact-span claim. The smaller test family also independently enumerates all
256 binary words and compares each candidate's cardinality and minimum distance.

This family cannot meet 240 by changing the coset subset. Any length-8 sign
configuration with the required Hamming distance at least 2 contains at most
128 words: partition the binary cube into 128 pairs that differ only in their
first bit, and at most one word from each pair is allowed. The best family
member reaches this sign-only ceiling. More search in this grammar is therefore
insufficient to close the target deficit.

There is a concrete missing construction step: assemble these 128 sign points
with the 112 support-two points in the same eight-dimensional space. At norm
8, use all vectors with two coordinates equal to ±2 and the others zero.
There are `4 * C(8,2) = 112` of them. Their mutual inner products are at most 4,
and their inner products with sign points are at most 4. The components have
different supports, so are disjoint. Their union is the known 240-point E8
configuration. The generic `geometry.orthogonal_union` places components
in separate spaces and cannot express this same-space assembly, especially
with an implicit sign component. A general compatibility/union backend and its
certificates would make this construction step searchable. The imported planner
patch did not implement that backend or claim discovery of E8. The current
integrated engine separately has `gluing.union` for explicit integer/Gram sources
and supplied cross inner products. It still does not accept the implicit sign
component used by this measured grammar, so this historical deficit persists.

## 28D orbit search: reproduced failure to reach target quality

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery plan discovery/recipes/planner_groups28.json \
  --state-dir /tmp/planner-groups28-complete --max-steps 100 \
  --max-evaluations 10 --time-limit-seconds 110 --per-program-seconds 50 \
  --budget max_pairs=30000000 --budget max_materialized_points=100000
```

The planner generated the representation-to-search edge and two search seeds
from an integer domain. The representation was fixed to `p=3, rank=3, parity=1`;
its choice was supplied, not discovered. Each search used 12 trials, a 1,600-line
cap and the existing finite phase-grid/operator defaults.

| Seed | Candidate points | Status | Closure failures | Materialization limit |
|---|---:|---|---:|---|
| 20261104 | 6,720 | numerical | 61 | reached |
| 20261105 | 6,720 | numerical | 65 | reached |

Both generated candidates completed the independent numerical check. Summed
worker time was 19.669 seconds. Neither is an exact-certified lower bound.
The best exact output table is empty, and the target-quality gate is false.
Against the **dated repository target of 221,988**, the observed deficit is
**215,268 points**, and the candidate reaches **3.027%** of the target. The
comparison is the supplied repository snapshot, not a fresh literature survey.

Both MILPs optimized the same 30-block restricted graph, with 387 conflicts and
finite numerical graph bound 6,720. An optimum on that graph cannot explain away
the missing 215,268 points: much of the construction space was never represented.
The worker diagnostics explicitly retain the materialization flags and closure
failures, so finished outer enumeration is not confused with complete internal
search.

These results expose a noncompetitive search family. They do not establish that
another seed, a larger budget, all representations, or every group construction
fails. Useful next instruments include richer representation/character spaces,
higher multiplicity parameters, orbit compression that avoids the observed
materialization losses, and general nonorthogonal mixed-component assembly.
The planner exposes these as missing/omitted capabilities and keeps their
unsearched scope explicit. It cannot synthesize mathematical operations absent
from its catalogue, and it is not yet a universal record-finding engine.

## Handoff limits

This change delivers actual typed DAG synthesis, parameter candidate generation,
quality-guided ordering, resumable state and competitive-quality diagnostics.
Its measured searches remain controls and failures against the quality targets;
**no new SOTA was found**. A successful record experiment still needs a richer
construction grammar, competitive search evidence, an independent exact
certificate and a fresh comparison review.
