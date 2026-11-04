# Typed construction-program planner

The planner synthesizes construction programs from an operator catalogue. A
caller supplies finite parameter domains and a node bound; the planner chooses
which operators to apply, their parameters, and which earlier results supply
each input. This differs from `search`, whose caller supplies graph slots and
edge alternatives, and `run`, whose caller supplies a complete DAG.

The goal is competitive new constructions. A correctly enumerated language or a
valid configuration below the target is not discovery success. The ledger has
separate `enumeration_complete`, `target_quality_passed`, candidate-quality and
exact-quality diagnostics. Neither completion nor a target match proves novelty.

## Run, resume and inspect

No installation is needed. Use the existing Python environment and dependencies.
The following search generates all 128 zero-containing subsets of eight even
syndromes of the extended Hamming code, and synthesizes their spherical images:

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery plan discovery/recipes/planner_cosets8.json \
  --state-dir /tmp/coset-planner --max-steps 25000 --max-evaluations 200 \
  --time-limit-seconds 120
python -m discovery plan-status /tmp/coset-planner
```

`--state-dir` must be new unless `--resume` is supplied. To try smaller batches:

```sh
python -m discovery plan discovery/recipes/planner_cosets8.json \
  --state-dir /tmp/coset-batches --max-evaluations 8
python -m discovery plan discovery/recipes/planner_cosets8.json \
  --state-dir /tmp/coset-batches --resume --max-evaluations 32
```

The ledger's `best_exact_by_dimension` names a `program_id`. Export it for the
ordinary independent verifier:

```sh
python -m discovery plan-export /tmp/coset-planner --program-id ID --out /tmp/candidate.json
python -m discovery verify /tmp/candidate.json
```

The existing comparison command still requires a dated catalogue entry; the
planner does not add entries or register records. All writes are local. There
are no network, publication, Actions, installation or push operations.

## Catalogue and finite parameter domains

```json
{
  "schema_version": 1,
  "max_nodes": 3,
  "strategy": "adaptive",
  "output_kinds": ["implicit_spherical_code"],
  "catalogue": [
    {"op": "cosets.build", "params": {"parity_check": [[1, 1, 1, 1]]}},
    {"op": "cosets.union", "domains": {"syndrome_labels": {
      "subsets": {"of": [1], "include": [0]}
    }}},
    {"op": "geometry.sign_embed"}
  ],
  "quality_targets": [{"dimension": 4, "cardinality": 24, "label": "known control"}],
  "unavailable": [{
    "name": "geometry.general_compatible_union", "family": "mixed",
    "reason": "Implicit sign components cannot yet feed the explicit coordinate/Gram gluing backend"
  }]
}
```

`params` fixes literal values; `domains` varies disjoint parameter names. Every
axis has exactly one of these forms:

- `{"values": [v1, v2]}`: literal alternatives, retaining order and duplicates.
  Array-valued parameters are individual elements of the outer list.
- `{"integer_range": {"start": 1, "stop": 10, "step": 2}}`: Python range
  semantics, with exclusive `stop` and optional nonzero `step` (default 1).
  Integer ranges are decoded directly without allocating all their elements.
- `{"subsets": {"of": [2,4,6], "include": [0], "min_size": 0, "max_size": 3}}`:
  generated subsets of `of` with fixed `include` elements. Size limits count
  chosen `of` elements; `include` does not count towards them. Elements must be
  distinct and the two lists disjoint. Sizes ascend by default;
  `largest_first: true` reverses size order. Within a size, combinations follow
  the order of `of`. Combinations are unranked without materializing the family.

Parameter names are sorted, with the final axis varying fastest. Entry order
and input-node tuple order are stable. Mathematical preconditions and unknown
operator parameters are checked at execution, not guessed by the planner.

Compilation precedes the execution timer and has explicit hard ceilings:
256 catalogue entries, 256 parameter axes across the catalogue, 4096 literal
alternatives per axis, 512 total `of`/`include` items per subset axis, and
1,000,000 aggregate subset-work units. A subset axis costs
`(len(of) + len(include) + 1) * (max_size - min_size + 1)` units; the work check
occurs before enumerating sizes or calling `math.comb`. Range integers and
variant counts have a 4096-bit ceiling, including each entry's Cartesian
product. Ranges such as `range(10**50)` remain lazy. Exceeding a compiler ceiling
is an explicit invalid-request error before checkpoint creation, not an
exclusion of the underlying mathematical family.

`operator_catalogue()` exposes all registered operations and their input/output
kind contracts, without a dimension selector. Existing contracts are in
`planner.py`; new backends can expose `OUTPUT_KINDS` beside `INPUT_KINDS`.
Each input position accepts either one kind or a finite list of alternatives,
matching the runner's contract. A producer whose possible outputs intersect that
list is bound once, even if several kinds match. This makes `gluing.union`
reachable from both `spherical_code` and `gram_spherical_code` producers, including
mixed pairs. Conditional producer outputs remain possibilities until execution.

`planner_supported` and `unsupported_contract_reasons` expose missing or invalid
input/output contracts in both `operator_catalogue()` and the `operators` CLI.
Selecting an operator with an unsupported contract fails explicitly. Unimplemented
capabilities are diagnostic declarations, never executable placeholders.
All registered but unselected operations appear as `registered_but_omitted`.

The current type system uses artifact kinds, not theorem-level refinements.
For example `groups.search` can return either `spherical_code` or `search_result`.
The planner keeps potentially compatible edges for either return; the runner
checks the actual kinds at runtime. A kind match also does not guarantee the
stronger mathematical conditions required by an operation: an orthogonal union
currently requires explicit integer coordinates, for example. Such failures
remain unresolved evaluations, not mathematical exclusions.

The current fallback contracts include roots, exact Gram codes, sections,
coordinate/Gram gluing, implicit stabilizer sectors, both selection operations,
and all three signed-layer operations. Independently verified coordinate,
implicit, Gram and quadratic spherical outputs can enter the exact ranking.
Quadratic output support also applies to the finite scheduler and `compare`.
This registration admits operators to the search language; it does not discover
their parameter values or certify every statically typed composition. In
particular, the explicit gluing backend does not accept implicit sign artifacts
and requires a supplied cross Gram matrix.

## Search and feedback

Starting from an empty graph, the planner appends every admitted operation with
all finite parameter variants and all positionally compatible tuples of earlier
nodes. Repeated inputs are allowed. Every nonempty prefix whose last kind can
match `output_kinds` is evaluated; other prefixes can still feed later
operations. At the node bound, non-target leaves are counted and pruned because
they cannot reach a requested output within the declared scope.

The finite language contains topologically labelled DAGs up to `max_nodes`.
Different node orders, duplicated nodes, unused earlier nodes and equivalent
constructions remain distinct. This deliberately avoids unjustified equivalence
or mathematical pruning, but can be expensive. No failure prunes descendants.

`breadth_first` is deterministic by generated program id. `adaptive` first
evaluates pending requested outputs and gives a newly generated prefix an
initial expansion opportunity. This initial descent has length at most
`max_nodes`. Suspended prefixes then receive alternating turns: oldest last
productive expansion first, or optimistic kind distance followed by observed
quality and stable id. A productive expansion stores one child and advances
the parent's persistent scheduling turn, yielding to other prefixes. The root
retains all constructor alternatives in its cursor. Quality feedback uses a
candidate's cardinality divided by a declared dimension-specific quality target;
without targets it uses `log2(cardinality + 1)`. Counts without a matching target
receive zero feedback when targets are present. Scores propagate to children.
Numerical and unverified counts can guide ordering, but only independently
exact-verified spherical outputs enter `best_exact_by_dimension` or pass an exact
quality target. Ordering can delay branches under resource bounds; nothing is
silently removed on the basis of dimension or score. This is an initial search
heuristic, not a learned construction policy or a claim of search efficiency.

This replaces the imported policy that generated all root alternatives first.
In the regression with ten representation constructor alternatives and a
two-node representation-to-search path, `max_programs=5` now reaches candidate
evaluations while preserving the remaining constructor cursor. Resuming with a
larger cap reaches all ten paths in the same order as an uninterrupted run.
The corresponding depth-three regression also reaches candidates under the
same five-program cap, but its first candidates can have three nodes: fresh
descent is not a shortest-program guarantee.
The queue does not promise a useful result before every cutoff. Non-target
leaf variants still consume steps without advancing the productive-expansion
clock, failed operators can consume their time caps, and very small limits can
stop before any candidate. Those remaining limits stay explicit in the ledger.

Dimensions are ordinary operator parameters or derived facts. The caller's
`RunContext` budgets still apply, including the explicit default dimension cap
of 256. A cap/failure is recorded and is never proof of impossibility. Code,
ring, group and geometry paths use the same synthesis mechanism. Their ability
to connect depends on available type-compatible operations; declaring a mixed
capability does not implement it.

## Checkpoints and interruption semantics

The single-writer `planner.sqlite3` database stores specification/source/contract
identity, all generated prefixes, each prefix's expansion cursor, per-entry
binding counts, candidate receipts, every evaluation attempt and invocation
limits. SQLite transactions with full synchronization commit a new graph and
the parent's cursor advancement and scheduling turn together. Checkpoint schema
version 2 stores each prefix's expansion turn and the shared expansion clock;
resume preserves adaptive interleaving as well as graph order. A nonblocking file lock rejects a
second writer. The OS releases the lock if the process exits.

An evaluation's `running` state commits before its worker starts. Its result
and current program state commit together. If interrupted before that commit,
resume marks the previous attempt `interrupted`, requeues the program, and
repeats its computation. A completed computation whose transaction did not
commit may therefore be repeated. There are no external publication side
effects to deduplicate. No committed frontier entry is lost or skipped.

Resume requires identical catalogue, node bound, targets, strategy, source
fingerprints, contracts, operator budgets and work directory. Use a new state
directory for a changed search scope or implementation. Invocation limits can
change: step/evaluation counts, stored-program count, wall times, program bytes
and receipt bytes. `--resume --retry-unresolved` explicitly requeues previously
failed, timed-out or resource-limited candidates while retaining their attempt
history. Ordinary resume does not silently retry completed failed attempts.

The source pin covers the discovery Python modules through the shared runner's
fingerprint manifest. The planner checks source hashes and contracts at loop
boundaries and after each worker. A completed worker receipt must also match
the saved source manifest. A detected change discards that candidate result,
keeps the attempt unresolved, and stops with
`source_or_contract_changed_during_invocation`. This detects ordinary edits;
it is not a hermetic execution environment, complete environment snapshot or
signed provenance. Checkpoints are local experiment data; mathematical
proof still comes from the exported receipt's independent verifier.

## Coverage and resource ledger

The ledger distinguishes:

- Generated, completed, failed, timed-out, resource-limited, pending and running
  programs, including unresolved reasons and interrupted attempt history.
- Typed extension opportunities and type-incompatible raw bindings for each
  catalogue entry at expanded prefixes; prefixes without compatible inputs;
  non-target leaves pruned at the declared node bound.
- Pending evaluations, unexpanded prefixes, counted immediate extensions, and
  prefixes whose extension counts are not computed. Descendant counts are
  reported as unknown until expanded, not invented as zero.
- Omitted registered operators, explicitly declared missing capabilities,
  parameters outside the finite domains and programs beyond the node bound.
- Operator-reported sampling/materialization limits, solver diagnostics,
  output restrictions and proof scopes. A finished numerical group-search
  worker does not establish exhaustive coverage of its mathematical family.
- Target cardinalities, observed candidate counts/statuses, exact counts,
  fractions and deficits. A completed search can still fail every quality test.

`enumeration_complete` means the declared finite frontier has been exhausted
and every target program attempted. `status: complete_within_declared_scope`
additionally requires that those evaluations returned normally. It says nothing
about exactness, competitive quality, internal numerical-search exhaustiveness,
novelty or unrestricted optimality. `target_quality_passed` is separate and is
null if no quality targets were declared.

CLI exit codes are 0 for completed enumeration/evaluations, 2 for incomplete
search and 1 for invalid input or storage/output errors. A quality deficit is
visible in the ledger even with exit 0; use `target_quality_passed` as the
separate quality gate in experiments.

`--max-steps` bounds frontier/evaluation actions, including non-target work.
`--max-evaluations` bounds worker launches per invocation. `--max-programs`
bounds persisted prefixes. Per-worker timeout is capped by remaining total
execution time; a disposable interpreter is killed at timeout. `max_program_bytes`
and `max_receipt_bytes` bound serialized program and retained receipt sizes.
Existing operator budgets bound supported allocations and proof work. SQLite,
serialization and process startup add overhead outside the execution deadline.
There is no hard OS RSS or filesystem quota: these limits are not promises
against an operating-system OOM or disk-full condition. Large catalogues can
produce many redundant graphs; checkpointing does not cure combinatorial growth.

## API and tests

```python
from pathlib import Path
from discovery.core import RunContext
from discovery.planner import plan_search, read_ledger, export_receipt

ledger = plan_search(spec, RunContext(Path('.')), '/tmp/new-state',
                     max_steps=1000, max_evaluations=20)
ledger = plan_search(spec, RunContext(Path('.')), '/tmp/new-state',
                     resume=True, max_evaluations=100)
```

`compile_catalogue`, `extension_counts` and `extension_at` expose the pure
synthesis functions. Tests compare a complete tiny language with an independent
Cartesian graph enumerator, subset domains with `itertools.combinations`,
uninterrupted runs with resumed runs, real worker limits with unresolved ledger
states, and generated coset codes with exhaustive binary-word enumeration.
They interrupt both graph generation and an evaluation after its durable start,
check writer exclusion and source/budget pins, verify adaptive quality feedback,
and exercise export through the ordinary verifier.

The realistic family measurements and structural deficits are recorded in
[PLANNER_VALIDATION.md](PLANNER_VALIDATION.md). They are construction-search
controls and diagnostics, not new record claims.
