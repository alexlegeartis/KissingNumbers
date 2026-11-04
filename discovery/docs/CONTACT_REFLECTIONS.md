# Verified contact-reflection search

`search.contact_reflections` reuses
[`research/contact_reflection_20261104.py`](../../research/contact_reflection_20261104.py).
The prototype's default standalone search is unchanged. Its `return_pool=True`
interface exposes the same exact generator and blocker lists to the engine.

For equal-norm vectors with squared norm N, a contact pair with product N/2
gives x-y of squared norm N. A pair with product -N/2 gives x+y. The generator
also considers both negatives, removes source overlaps and duplicate candidates,
and optionally filters by an explicit blocker cutoff. This expands the existing
negative-contact-sum operator, especially for asymmetric sources.

The engine adapter accepts an explicit integer spherical code with a positive
integral diagonal metric, encoded by `weights`. It independently replays source
geometry and positivity before invoking the prototype. The standalone prototype
also supports symmetric integer metrics under its documented positive-definiteness
and source-validity assumptions; this adapter does not silently transfer those
assumptions to an exact engine result.

Every source point can be removed. Shared incumbent deletions are charged once.
When a checked matching assigns every generated candidate a distinct mandatory
source blocker, no subset can gain within that generated pool. This presolver
does not use candidate conflicts or numerical dual bounds. Otherwise the engine
constructs the complete exact candidate-conflict graph and calls its existing
shared-blocker selector. Finally it checks every output norm, point distinction,
pair inequality and exact rank with the independent coordinate verifier.

Parameters are `seeds`, `pair_limit`, `candidate_limit`, `blocker_limit` and
`seconds`. Defaults use at most the first 200 source indices, 100,000 directed
seed/source pairs, and 256 generated candidates. A supplied seed list must have
distinct valid indices. A cutoff, incomplete seed list or generation cap defines
a restricted search. No timeout, numerical result or matching in that restricted
pool establishes global optimality. The receipt records limits, truncation,
source and prototype hashes, selected points, shared deletions and any matching.
Generation provenance survives independent terminal verification in a DAG.

Dimension, materialized-point count, generation/screening/pair work, seed count
and solver time are bounded by the existing execution context. The prototype
uses proven safe int64 arithmetic or Python integers. The adapter always returns
a completely reverified ordinary coordinate artifact and never declares novelty.

The strict regression uses norm 6 and the identity metric:

```
(-1,-2,1), (-1,1,-2), (1,1,2),
(-1,2,1), (2,-1,-1), (-2,-1,-1).
```

The old sum-only operator produces no new candidate. Reflection search removes
the third point and inserts `(2,1,1)`, `(1,-1,2)` and `(1,2,-1)`, reaching eight
points. Every norm and pair is checked in integer arithmetic. This is a useful
discovery regression, not a record claim. Further tests cover a weighted hexagon,
arbitrary precision, exact blocker matching, invalid cached receipts, resource
preflights, truncation, empty seeds and retained DAG provenance.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  python -m unittest discovery.tests.test_contact_reflections
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  python -O -m unittest discovery.tests.test_contact_reflections
```
