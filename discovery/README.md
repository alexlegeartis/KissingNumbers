# Mathematical discovery workbench

An experimental, composable foundation for searching and checking kissing
configurations. It is separate from the repository's registered record pipeline.
**The work so far reproduces known constructions and records bounded negative
searches. It has not found a new SOTA.**

The Sol continuation adds three target-independent operators: support-orbit
packing, arbitrary integer point-pool selection across any number of cohorts,
and negative-contact-sum generation with joint selection. A subsequent extension
adds exterior-square representations, exact rational orthogonal orbits and
weighted quadratic-coordinate input. Joint character-coset and sign-orbit
assembly now replace sparse bases and dense caps together. There are now 37
registered operators. The separate 1–24 held-out benchmark matches ten targets
with formula controls and four through actual support search, including 9D306.
The 13D/14D/31D record attempts have no improved certified witness. See
[construction search](docs/CONSTRUCTION_SEARCH.md) and the
[continuation evidence](../research/sol_discovery_20261104/README.md).
The [large known controls](../research/orbit_controls_20261104/README.md)
now include17D5730,18D8358,22D49896,23D93150 and24D196560 through the same
registered engine, with explicit resource profiles and independent proof replay.
The [coupled-family campaign](../research/character_families_20261104/README.md)
adds a rational19D12460 control with a strict auxiliary gap. Its17D/19D/21D
construction searches have no certified improvement. The
[character-family operators](docs/CHARACTER_FAMILIES.md) retain whole-family
selection and numerical shape deformation as distinct proposal mechanisms.

Run from the repository root with the existing Python environment; no installation
is required. Dependencies are NumPy, SciPy and SymPy. The code makes no network
calls, installs nothing, and never pushes, publishes, or edits the record registry.

```sh
python -m discovery operators
python -m discovery run discovery/recipes/ganzhinov10.json --out /tmp/known10.json
python -m discovery verify /tmp/known10.json
python -m discovery compare /tmp/known10.json
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s discovery/tests -t . -v
```

`run` refuses to overwrite an existing receipt. Pick a new output filename for a
repeat. `verify` exits 0 only for an independently exact-verified result, 2 for a
numerical/unverified result, and 1 for invalid input. `compare` re-verifies first
and compares only to a dated baseline; exceeding it is not a worldwide priority
claim. Library versions and numerical tie-breaking can affect sampled searches.

## Construction programs

A version-1 JSON program is a directed acyclic graph. Each node has an `id`, an
allow-listed `op`, positional input node ids and operator-specific `params`.
Forward references are allowed; missing inputs, cycles, type mismatches, unknown
operators and unknown parameters are rejected. Recipes cannot import or execute
arbitrary code.

```json
{
  "schema_version": 1,
  "nodes": [
    {"id": "factor", "op": "rings.binary_factor", "params": {"n": 7, "index": 0}},
    {"id": "lift", "op": "rings.hensel", "inputs": ["factor"], "params": {"modulus": 4}},
    {"id": "ring", "op": "rings.extend", "inputs": ["lift"]},
    {"id": "binary", "op": "rings.gray", "inputs": ["ring"]},
    {"id": "sphere", "op": "geometry.sign_embed", "inputs": ["binary"]}
  ],
  "outputs": ["sphere"]
}
```

This reconstructs the known 256-word, length-16, distance-6 code and keeps its
spherical sign image implicit. It is a control, not a record.

An `Artifact` has a kind, a data recipe, computed facts, and a verification
receipt. `RunContext` carries resource limits. `OPERATORS` and `INPUT_KINDS` are
the module extension points. A new backend also needs an independent verifier
registered in `certification.py`, proof-scope documentation and corruption tests.
No unimplemented instrument is registered as a placeholder operator.

## Implemented instruments

| Module | Working capabilities | Important restrictions |
|---|---|---|
| `engines/groups.py` | Supplied unitary generators and Fourier/parity representations, sampled subgroup-character eigenlines, actual vector orbits, relative phases, compressed cross-orbit screening, finite block MILP; exact Ganzhinov controls | Numerical discovery; finite phase grid; sampled subgroups; orbit and output caps; no automatic certification of approximate output |
| `engines/rings.py` | Binary factorization, coprime Hensel lifts, cyclic free modules, extensions, Gray isometry, exhaustive small-code distance, two independent DGLS proofs, compatible repeated-mask unions | Half-degree irreducible factors; moduli 4/8/16; small enumeration caps; restricted greedy translate selection |
| `engines/cosets.py` | Full-rank binary parity checks, exact syndrome leader distances, separate nonzero kernel-distance DP, specified coset unions | At most 16 syndrome bits and bounded transitions; not a generic high-distance 96D code solver |
| `engines/layers.py` | Native compressed support/sign-layer verification; known 352800-point 32D control; parametrized supported family | Dense binary span must fit enumeration budget; sparse support and exceptional-line pair budgets; no general layer search optimizer |
| `engines/geometry.py` | Exact integer and rational Gram input, sign embeddings, orthogonal unions and exact central hyperplane sections | Generic implicit mixed-layer gluing remains incomplete; bounded heterogeneous gluing is being independently reviewed |
| `certification.py` | Exact weighted integer Gram checks, distinctness, norm and rank; independent backend proof replay | Explicit coordinates need a bounded all-pairs check; large objects require a supported compressed proof |
| `runner.py`, `__main__.py` | Typed DAG runner, CLI, per-node provenance and timing, hashed receipts, dated comparison | No distributed scheduler, persistent search cache or autonomous theorem invention |

The `groups.ganzhinov14` control replays the E7 frame construction. It does not
claim autonomous discovery of the U(3,3) subgroup characters. The 10D numerical
search does recover the required line systems from generators without record
coordinates. Its separate exact control reconstructs the algebraic result.

## Proof obligations shared by all branches

- Establish distinctness and cardinality. For implicit families, prove injectivity
  or count disjoint cosets; a sum of component sizes alone is insufficient.
- Establish norms and a sufficient ambient dimension. Integer certificates also
  compute exact span rank. Sign images explicitly report ambient dimension;
  they do not claim a computed full rank.
- Establish every within-component and cross-component separation condition.
  Group representative checks use a common unitary action. Code unions use
  within-code distance and cross-coset bounds. Orthogonal sums use zero cross dots.
- Preserve exactness through the embedding. Binary length N and distance d give
  norm N and maximum inner product at most N−2d, requiring 4d≥N for kissing.
- Keep provenance distinct from mathematical validity and novelty.

Saved `verification.status` values are not proof. `verify` recomputes the relevant
checks, rejects false geometric facts, and never promotes numerical coordinates
to exact. Integer accumulation uses a proven overflow bound or Python integers.
The legacy DGLS routines use assertions, so their adapter refuses optimized
Python (`-O`) rather than silently disabling those checks.

Receipt hashes detect accidental modification of the program, artifacts and
events. They are not signatures or authenticated provenance. Source hashes record
the implementation used; mathematical verification runs with the current code.
A valid coordinate certificate alone does not establish the historical process
by which its coordinates were discovered.

## Dimension-independent interfaces and practical scaling

Dimensions are parameters and certificate facts, not a hard-coded 5–16 search
range. The default dimension budget is 256, which accommodates constructions
through 96D. This does **not** mean every implemented algorithm is competitive or
tractable in every dimension.

The group branch compresses a cross-orbit check from |A||B| pairs to one seed
against B for each constituent transitive orbit. It still enumerates the orbit;
there is no general double-coset backend yet. A final numerical coordinate output
is checked in bounded blocks.

The syndrome branch enumerates 2^r quotient states instead of 2^(n−r) kernel
words. It must separately certify the shortest nonempty zero-syndrome support:
the empty word's distance zero cannot certify the kernel's minimum distance.
Its 16-check cap cannot represent every useful high-distance code. Large
high-distance examples need stronger algebraic proof backends such as DGLS.

The DGLS branch keeps **137438953472** length-96 binary words implicit. It verifies
the synthesized ring generator using two independent native proofs, each checking
2^24 lower-ring words and 759 small coset exclusions. No 2^37-word array or giant
Gram matrix is allocated. This validates compressed composition, not a new
96D construction. The full existing 96D layered record remains larger.

The `layers.known32` control independently verifies the existing 352800-point 32D
construction using 131072 dense codewords and 1826076 compressed support/sign
checks. It allocates no spherical point array. `layers.signed_supports` verifies
the same explicitly defined layer family with dimension parameters through 96;
this is a restricted construction grammar, not a generic mixed-layer optimizer.

Run the larger group trial with explicitly sufficient budgets:

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery run discovery/recipes/search28.json \
  --budget max_materialized_points=100000 --budget max_pairs=30000000 \
  --out /tmp/search28.json
```

Hard algorithm ceilings and user budgets are both enforced. Budget exhaustion is
recorded or rejected; it is never reported as a mathematical impossibility.
Large unverified ring recipes remain structural candidates. High-degree group
representations can have enormous orbits even when their dimension is modest.

## Benchmarks and continued search

Recipes include `ganzhinov10.json`, `ganzhinov14.json`, `discover10.json`,
`search26.json`, `search28.json`, `octacode.json`, `dgls96.json`, `ring_sign32.json`,
`syndrome_cosets.json`, `known32.json`, and `orthogonal_union.json`. Known-control cardinalities
must not be confused with competitive new records: for example the ring32 sign
layer is only 8192 points, below the recorded 352800 in that dimension.

The dated [catalog](benchmarks/catalog.json) separates benchmark expectations from
comparison records. [Research evidence](docs/RESEARCH.md) records an additional
exact negative search: all single-coordinate linear extensions of two fixed Z4
modules are excluded from meeting the 36D kissing threshold by small witness sets.
Replay it with `python discovery/benchmarks/replay_extensions.py`.

The earlier exploratory files and execution evidence remain in
`research/group_search_20261104/`. Production operators do not import those
prototype files. The engine is additive; `audit.py`, `run_all.py`, and the
registered construction packages retain their existing behavior.

## Remaining research work

Priorities include richer group catalogues and subgroup projectors; higher
multiplicity parameter spaces; continuous phase solving; double-coset proof
compression; generalized concatenation, switching and gluing; character-weighted
partial weight-enumerator recurrences; and generic section/mixed-layer search.
The DGLS recurrence is not an unsigned disjoint-support count. These are missing
instruments, not executable stubs or features already delivered.

New records remain the objective. Neither a successful known-result benchmark,
a finite MILP optimum, nor a failed bounded search establishes a new SOTA or an
unrestricted optimum. A candidate must pass exact verification and a fresh
comparison review before it enters the record pipeline.

## Measured end-to-end controls

`benchmarks/run_controls.py` executes the 10D, 14D, 352800-point 32D and implicit
DGLS 96D recipes through the shared runner and re-verifies terminal results. It
saves separate receipts plus timings and cumulative process peak memory:

```sh
OPENBLAS_NUM_THREADS=1 python -m discovery.benchmarks.run_controls --outdir /tmp/discovery-controls
```

Choose a fresh directory for every run. `benchmarks/replay_extensions.py` is the
separate pure-integer replay of the bounded 36D negative search.

Measured test, benchmark and search outcomes are in [docs/VALIDATION.md](docs/VALIDATION.md).

An independent pre-handoff review and its fixes are documented in [docs/FINAL_REVIEW.md](docs/FINAL_REVIEW.md).

## Subsequent additive increments

The reviewed foundation and finite scheduler were pushed separately as
`e25c2924f377a1061b6c60429c950d34df7515e9`. Newer working-tree increments are
documented in [root systems](docs/ROOT_SYSTEMS.md),
[Gram certificates and sections](docs/GRAM_SECTIONS.md), and
[object/contract schemas](docs/OBJECT_CONTRACTS.md). Their exact coverage is
reported by `research/discovery_campaign_20261104/coverage_1_24/build_matrix.py`;
a working recipe is not automatically a blind discovery or a new record.

The full existing96D layered construction, not just its DGLS ingredient, has a
specialized dual-checker benchmark:

```sh
python -m discovery.benchmarks.run_full96 --outdir /tmp/new-full96-replay
```

It replays both existing compressed proofs and requires agreement on all five
layers and142162362560 points. This is a known-construction benchmark bridge,
not a generic layer-synthesis operator.

## Typed construction planner

`python -m discovery plan` synthesizes bounded construction DAGs from a finite
operator/parameter catalogue, choosing operators and their input edges. It
checkpoints the frontier in SQLite for interruption recovery. `plan-status`
inspects the ledger and `plan-export` exports a normal receipt for independent
verification. See [typed planner](docs/TYPED_PLANNER.md) for the interface,
contracts, resource limits and exact scope.

Current contracts include coordinate, implicit, Gram and quadratic spherical
outputs and the runner's alternative input kinds, including mixed explicit
gluing. Missing contracts are visible in `operators` and rejected explicitly.
Quality targets affect ordering and reporting without filtering dimensions.
The imported search measurements remain failed quality gates: 8D 128 versus
240 and numerical 28D 6720 versus 221988. Enumeration and passing tests do not
establish competitive discovery. See [historical validation](docs/PLANNER_VALIDATION.md)
and the separate [integration checks](benchmarks/planner_integration_20261104.json).
