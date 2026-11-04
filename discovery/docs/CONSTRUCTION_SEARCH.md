# Target-independent construction search

The three operators in `engines/construction_search.py` generate and assemble
candidate configurations. They use the existing DAG, budgets, typed planner and
independent exact coordinate verifier. All 32 registered operators have planner
input/output signatures. No dimension-to-record dispatch or target-count table
is used by the new construction operators.

`search.support_packing` searches four-coordinate supports such that no triple
is repeated. An optional permutation, described by cycle lengths, compresses
supports into cyclic orbits. It discards internally conflicting orbits, solves
the weighted packing problem and checks the proposed discrete incidence before
expanding all 16 signs and adding the coordinate axes. Search inputs contain
dimension, cycles and a time budget. They contain no support design, record
coordinates or expected output count. A deterministic feasible packing is kept
if the numerical optimizer fails to supply a better incumbent.

`selection.integer_pool` accepts any finite common-norm integer point pool in
a positive integral diagonal metric. It deduplicates points, constructs every
conflict with exact integer arithmetic, and selects points individually. There
may be any number of mutually conflicting cohorts. Optional initial indices
supply a feasible incumbent; all those points can be removed. The numerical
MILP is a change of variables between insertions and deletions. Its rounded
incumbent is checked against the complete exact conflict graph and then every
output norm, pair, distinctness and exact rank is independently recomputed.
Solver bounds are diagnostic and never certify optimality or impossibility.

`search.contact_exchange` generates `x+y` whenever two source points satisfy
`<x,y> = -norm/2`. The new vector has the same norm. It deduplicates these
proposals and jointly selects from old and new points with the pool optimizer.
This is an algebraic shell-completion mechanism; it does not guarantee a valid
union. A candidate limit is explicit and incomplete generation is reported.

Pool shape, dimension, pair count, temporary Gram entries and integer sizes are
bounded before corresponding allocations. Support enumeration and output have
separate preflights. Large integer pool products fall back to Python integers.
Numerical optimization can affect cardinality and runtime across installations;
the exact mathematical validity of a returned coordinate artifact is replayed
without trusting the optimization receipt.

## Held-out measurement

The benchmark executes 200 programs across dimensions 1 through 24. Its search
command has no comparison-source argument and loads no incumbent coordinates.
A separate evaluation command reads the captured table and independently
replays every saved result. Formula controls and actual support search have
separate result fields. This finite catalogue is intentionally a measured
baseline, not a claim of universal competitive coverage.

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery.benchmarks.blind_rediscovery search \
  --seconds .4 --out /tmp/new-blind-search
OPENBLAS_NUM_THREADS=1 python -m discovery.benchmarks.blind_rediscovery evaluate \
  /tmp/new-blind-search \
  --source research/discovery_campaign_20261104/baselines/cohn.html \
  --metadata research/discovery_campaign_20261104/baselines/retrieval.json
```

The saved run matches targets in 1–9 and 16. Actual support search matches
4,7,8,9; 9D306 is a held-out rediscovery without record coordinates. Other
matches include classical root or stabilizer formula controls. Search also
reaches historical 12D840, below current 841. The 24D result is only 7040 versus
196560, a substantial failure of this small catalogue. Existing exact target
replay coverage is broader: 1–10,13,14,16. Replaying a known target and finding a
target through held-out search are different measurements.

## Growth diagnostic and competitive trials

The evaluator parses all 48 neighboring lower and upper bounds, treating a
`colspan=2` exact-value cell as both bounds. For each n it computes exact
rational lower and upper growth factors and their quotient:

`score(n) = (U[n]/U[n-1]) / (L[n]/L[n-1])`.

A value above one means the upper/lower gap expands in that step. A centered
two-step quotient is also saved. This is a search-order heuristic, not evidence
that a larger construction exists. The largest snapshot scores are 31D and
37D: 1.421043 and 1.410351. The source is the SHA-verified Cohn capture retrieved
2026-11-04 at 12:36 UTC. Fresh browser/executor requests in this continuation
failed with HTTP 502/503; no fresh worldwide-priority verification is claimed.

The continuation generated and searched 13D contact pools, ran a registered
14D contact-exchange DAG, and generated actual 31D core–tail reflections and
auxiliary contact sums. There is no improved certified lower-bound witness.
The [campaign report](../../research/sol_discovery_20261104/README.md) records
finite-pool certificates, unresolved numerical searches and reproduction.

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest \
  discovery.tests.test_construction_search discovery.tests.test_blind_benchmark
OPENBLAS_NUM_THREADS=1 python -O -m unittest \
  discovery.tests.test_construction_search discovery.tests.test_blind_benchmark
```

These tests cover held-out support construction, current versus historical
targets, complete small-pool comparison to brute force, genuine incumbent
replacement, overlap removal, three-cohort assembly, weighted geometry,
arbitrary precision, shape/budget rejection and exact-table parsing.
