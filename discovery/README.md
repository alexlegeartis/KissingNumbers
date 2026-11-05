# Mathematical discovery workbench

An experimental, composable foundation for searching and checking kissing
configurations. It is separate from the repository's registered record pipeline.
**The generic engine now generates a certified 43D improvement from 2,919,782
to 2,919,785 points.** This seeded body/head search independently reproduces the
separate 43D task's discovery. Its full replays retain the documented P48n and
complete-family population inputs; there is no held-out or worldwide-priority
claim.

The [algebraic continuation](../research/algebraic_switches_20261104/README.md)
finds **209,825 points in 27D** (+18) and **251,196 in 31D** (+8). Complete
explicit geometry and independent direct-integer changed-pair proofs pass.
A reusable [shared-blocker selector](docs/BLOCKER_EXCHANGE.md) assembles partial
orbits and charges shared incumbent deletions once.

The Sol continuation adds three target-independent operators: support-orbit
packing, arbitrary integer point-pool selection across any number of cohorts,
and negative-contact-sum generation with joint selection. A subsequent extension
adds exterior-square representations, exact rational orthogonal orbits and
weighted quadratic-coordinate input. Joint character-coset and sign-orbit
assembly now replace sparse bases and dense caps together. A finite body/head
operator finds shared-deletion exchanges. A [partial sign-orbit operator](docs/PARTIAL_ORBITS.md)
uses exact group kernels while selecting individual points and changing all
incumbent layers jointly. The [contact-reflection adapter](docs/CONTACT_REFLECTIONS.md)
reuses the shared prototype, checks the source independently, and verifies every
output pair. There are now 41
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
The [high-dimensional campaign](../research/highdim_characters_20261104/README.md)
records the 43D improvement, negative 33D–35D/39D/40D experiments and independent
full verification. See [body/head exchanges](docs/BODY_HEAD_EXCHANGE.md).

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
| `engines/cosets.py` | Full-rank binary parity checks, exact syndrome leader distances, separate nonzero kernel-distance DP, specified coset unions | At most 19 syndrome bits and bounded transitions; not a generic high-distance 96D code solver |
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

### Binary translates of nonlinear block images

`search.binary_translates` proposes binary shifts of a ring code's block
image. It works in the binary image itself: it does not assume that a Gray
map turns ring addition into XOR. A mixed-integer model chooses one permitted
binary block per coordinate. An adversarial oracle seeks two ring messages
whose images produce a forbidden cross-layer distance. Every discovered pair
adds an exact, globally valid constraint on subsequent proposals. A local
oracle is incomplete: surviving proposals are **unverified**, not distance
certificates. For sufficiently small ring codes the oracle can instead
exhaustively check all message pairs; even then, internal distance, cardinality,
and the spherical embedding remain separate obligations.

`recipes/dgls96_binary_translate_probe.json` applies this search to the Z8
Golay generator in the Duursma–Greferath–Litsyn–Schmidt construction. Its
output shifts `u` have even parity in every four-bit block. XORing `u` with a
fixed one-bit-per-block vector gives a prospective third layer: its distance to
the unshifted even-block layer is automatically at least 24, but its distance
to the existing odd-block layer still requires a complete proof. The local
oracle only finds counterexamples; it does not supply that proof. No record
follows from running the recipe. The search receipt independently replays its
negative witnesses and any explicitly exhaustive cross checks, while remaining
unverified as a spherical construction.

Candidate and reference ring offsets may be specified separately. This lets the
search split a ring code into affine pieces before applying binary shifts,
without incorrectly treating ring offsets as XOR shifts. For example,
`recipes/dgls96_index_two_probe.json` considers one 2^35-word subcode against
both affine halves of the original 2^36-word code. If a valid extra shifted
half-layer were ever independently proved compatible with the two existing
full layers, the five half-layers would contain 171,798,691,840 words. The
recipe searches for such a shift; it does not establish that one exists.

`recipes/dgls96_mixed_parity_probe.json` removes the even-block restriction
and tests arbitrary binary shifts of the index-two subcode against all four
existing affine half-layers. This is a wider search, not an assertion that an
additional layer exists.

`search.block_isometries` goes beyond binary shifts. It jointly selects
user-specified local maps, checking that each map preserves every block-level
Hamming distance before it is admitted. Coordinates can be assigned to groups
that share one choice. The search proposes all group choices simultaneously,
then uses explicit bad word pairs to add exact constraints. The DGLS recipes
`dgls96_block_isometry_probe.json` and `dgls96_subcode_isometry_probe.json`
use the 384 isometries of the eight even-parity four-bit words. The first
applies them to the full 2^36-word image; the second works on an index-two
subcode. Neither recipe proves inter-layer compatibility or a new bound.

`search.signed_fiber_exchanges` constructs and optimizes pointwise signed
antipodal-line exchanges around a compressed signed-support construction.
Selected full fibers can independently choose even or odd sign parity, and
individual old and new lines can be removed or added together. The returned
construction is replayed by the exact compressed verifier; its search
provenance reports the numerical solver's status so that a time limit is not
mistaken for an exhaustion proof. Candidate pools and flexible fibers bound
the scope of a search; an unchanged result is not a global impossibility.

The signed-fiber search uses proved support-level cuts as well as individual
conflict constraints. A set of signed octads on nine coordinates has at most
64 antipodal lines, and any octads sharing a fixed six-set also have at most
64. These cuts reduce the solver's search space without ruling out changes to
other supports or layers. `dim33_parity_exchange_probe.json` is a reproducible
bounded search on the 33D source, not a global optimization certificate.
A finite candidate pool can also be excluded from positive gain when its
supports match injectively to distinct conflicting incumbent fibers; this
structural result is recorded separately from a numerical solver status.
`min_new_lines` and `return_neutral` permit exact same-cardinality outputs
for subsequent searches. Such outputs are not automatically claimed to be
nonisometric or even distinct without a separate comparison.

`search.ring_orbit_assemblies` searches a different kind of algebraic
exchange. For a free Z8 code C and a diagonal 1/5 scaling plus shift, both
codes decompose into orbits of the common subgroup 2C. The operator seeks
explicit cross-orbit conflict witnesses, expands them by exact common
symmetries, and optimizes the resulting bipartite graph. Any matching of
verified conflict edges gives a rigorous upper bound for the specified
whole-orbit pool. Conversely, its independent set in the **sampled** graph is
only a proposal: missing conflicts, parent-code distance and the final sphere
embedding must still be checked. `dgls96_ring_orbit_replay.json` seeds a
specific search with 64 witnesses and replays the 4,096-edge matching; it is
a restriction for that transformation, not a global no-go theorem.

`search.group_orbit_seeds` jointly optimizes seeds for several orbits
of an exactly closed finite rational orthogonal group. It uses all relative
group products rather than assuming different orbits are independent. Its
numerical output is a proposal, and its orbit cardinality is only a
forecast. An optional common stabilizer subgroup restricts all seeds to its
exact rational fixed space; unintended additional stabilizers can shrink orbits
and are conservatively treated as conflicts by the numerical objective. `groups.certify_multi_orbit` reconstructs the group, deduplicates
each orbit and checks every relative product using exact integer arithmetic,
including unequal seed norms. `group14_search_strict504.json` demonstrates a
504-point exact orbit below the known 14D benchmark; the 13D and 14D multi-
orbit probes are exploratory. The seed-search operator does not yet support a different stabilizer per orbit.

`search.coset_group_orbits` constructs another rational orthogonal
representation directly from a finite group's action on cosets of a supplied
subgroup. For a square subgroup index it removes the trivial direction using
an exact rational change of basis, and can take several copies plus trivial
coordinates. It then jointly searches orbit seeds in the resulting
representation. `group16_coset_rep_orbit.json` independently certifies a
504-point 16D orbit, well below the 4,320 benchmark; the nine-orbit recipe is
a search probe, not a certified bound. The operator also supports a supplied
supplemental representation after checking all source-group relations, so
inequivalent representations can be combined. The
15D example combines an 8D coset representation with a 7D representation and
certifies an orbit of 504 points, below the 2,564 benchmark. Nonsquare-index
standard representations are not yet synthesized by this operator.

`search.mixed_stabilizer_orbits` searches several orbits of the same group
while allowing each seed its own stabilizer subgroup. It restricts each seed
to an exactly computed rational fixed space, but its numerical output is only
a proposal. `group14_mixed_stabilizer_benchmark.json` supplies 42 known seed
orbits of sizes 56, 28 and 14 and exactly replays the known 1,932-point 14D
configuration. Because the seeds are supplied, this is a regression benchmark,
not an independent rediscovery. `group14_mixed_stabilizer_extension_probe.json`
adds a generic orbit and jointly searches for an extension; it is not a proof
of a new bound. Different group orbits can still have more symmetry than their
chosen stabilizers, and the promotion step checks actual distinct counts.

`search.partial_group_assemblies` permits each orbit to use a different
rational group, including the identity group for independently movable points.
It can therefore break a selected orbit while leaving other parts symmetric.
An optional positive integer diagonal metric supports weighted-coordinate
constructions; every group must preserve that metric, and promotion checks
all materialized point pairs exactly, including unequal norms. Because this
operator builds a dense numerical Gram matrix, its Gram-entry budget scales
quadratically with the proposed point count.

The `group14_partial_orbit_benchmark.json` and
`group15_weighted_partial_benchmark.json` recipes replay supplied seeds for the
known 1,932- and 2,564-point configurations, respectively. They are calibration
checks, not discoveries. The companion `*_partial_extension_probe.json`
recipes add a singleton while breaking blocker orbits; their output is
numerical and must not be treated as a valid bound without exact promotion.
Use `--budget max_gram_entries=4000000` for the 14D recipes; the 15D recipes
need `--budget max_gram_entries=7000000 --budget max_nodes=150`.

`search.quadratic_group_assemblies` extends partially symmetric search to
coordinates in Q(sqrt(3)). It can start from algebraic seeds, numerically move
them, and produce provisional algebraic seeds for exact promotion. The
independent verifier compares all inner products using integer sign tests in
the quadratic field, including unequal algebraic norms; larger numerators
use a budgeted Python-integer fallback to avoid overflow. The supplied
`group13_quadratic_benchmark.json` exactly replays the known 1,154-point
configuration from supplied coordinates; the related extension probe is a
numerical search, not a new record. Both require `--budget max_nodes=1200`.

`group16_partial_benchmark.json` replays a supplied 4,320-point 16D code;
`group16_partial_extension_probe.json` splits selected blocker orbits and
searches for one additional point. Run these with `--budget max_nodes=250
--budget max_gram_entries=20000000`. The latter is an unverified numerical
search, and the former is not a rediscovery.

Exact rational diagonal generators may be written compactly as
`{"diagonal": [...], "denominator": 1}` instead of a dense matrix. The same
orthogonality and metric-preservation checks apply.

`search.cavity_repacking` jointly moves and replaces a selected neighborhood
while holding the remaining points fixed, without requiring any symmetry. It
supports weighted integer source coordinates, arbitrary numerical replacement
seeds, and exact promotion with unequal norms. The 17D cavity probe starts
from a supplied 5,730-point source, removes 308 points, and searches for 309
replacements. Its output is a numerical proposal, not a 5,731-point bound; run
with `--budget max_nodes=350 --budget max_pairs=17000000`. Exact promotion also
rechecks every fixed-fixed pair, so the source's validity is not assumed.

The joint optimizer evaluates the symmetric moving-point Gram matrix in place
and avoids explicit triangular-index and second-weight matrices. It still has
quadratic time and memory cost in the number of simultaneously moving points;
set resource budgets explicitly for large whole-configuration searches.

`search.adaptive_cavity` takes a proposed point, identifies its current
blockers, jointly optimizes their replacement and the new point, and can
release additional near-contact source points in subsequent rounds. The
17D adaptive probe demonstrates this automated neighborhood expansion. Its
results are numerical proposals, and it neither certifies fixed-fixed pairs
nor registers a new bound without exact promotion.

For unequal integer norms, the cavity verifier now uses integer square-root
lower bounds and checks only ambiguous pairs with arbitrary-precision integer
arithmetic. This avoids a quadratic table over all distinct norms while
retaining exact separation checks; the full set of pairs is still covered.

Finite-field discovery: `finite_field.support_pool` enumerates coefficient grids for permutation polynomials over small binary fields, images affine subspaces and filters the resulting support masks. `finite_field.octad_exchange` then searches individual antipodal signed lines in both sign parities, with individual source-line removals. Both return proposals, not automatic record claims. A verified 32D seed and a runnable, independently checked workflow are in `research/gf32_polynomial_octads_20261105/`.
The field-support generator accepts multiple explicit coordinate injections, so GF(2^m) images can occupy noncontiguous coordinates in an ambient space such as 96D. This is an embedding mechanism, not a claim of complete coverage or of an improved bound in that dimension.

Representation discovery: `search.spectral_graph` constructs exact rational or real-quadratic spectral projectors from the eigenspaces of a supplied symmetric graph, searches unions of sectors and optionally antipodal points, and independently certifies any resulting Gram code. It currently supports integer spectra or spectra contained in one real quadratic field, and is limited by the exact matrix and solver budgets; it is not a general character-table or arbitrary algebraic-number representation engine.
Polynomial maps need not be globally bijective: with `allow_nonpermutations`, the field generator retains an image only when the restriction to the particular source flat is injective. This enlarges the candidate grammar without treating a collapsed image as an octad.

Nonlinear trace-code discovery: `galois.kerdock` synthesizes Z4 generators from Teichmuller trace over GR(4,3) or GR(4,5), and exactly certifies their Gray-image distances by Lee-weight enumeration. `galois.evaluation_search` searches order-two or arbitrary non-Teichmuller evaluation points, `galois.find_translates` constructs exactly verified binary translate unions, and `galois.translate_search` bridges to the general translate proposal engine. See `research/galois_kerdock_20261105/` for limitations and reproducible controls. This is not a reproduction of the 96-bit DGLS parent.
For nonlinear binary parents too wide for exhaustive ambient-cube search, `galois.random_translates` samples shifts but verifies their cross-distances against the complete binary XOR difference set. Its output is an implicit, exactly replayable spherical construction; it does not assume that Gray images are binary-linear.

Pointwise projective exchange: `search.projective_exchange` builds an integer-exact conflict graph for supplied integer vectors with different norms under a positive diagonal metric. It numerically searches for simultaneous source deletions and candidate additions; `geometry.certify_projective_exchange` independently replays the proposed result exactly. It is a finite-pool search, and time limits or failures are not nonexistence proofs. See `research/e7_projective13_20261105/` for a reproducible 13D experiment.
Quadratic pointwise exchanges: `search.quadratic_exchange` extends the finite-pool search to integer-coefficient vectors in Q(sqrt(3)) with positive rational squared norms, including unequal norms. Its independently checked output is a `quadratic_spherical_code`; see `research/quadratic_exchange13_20261105/`. Solver infeasibility and time-limit messages are not exact certificates.
Both `search.quadratic_exchange` and `search.projective_exchange` accept `decompose: true` for a one-point target: it separates the candidate/source conflict graph into connected components and seeks a positive exchange in each. This exact decomposition changes search efficiency, not proof requirements; numerical solver conclusions remain provisional.

Whole-family construction search: `search.support_recombination` applies supplied coordinate permutations to an existing constant-weight support code, then recombines complete families with exact five-subset packing constraints and retained protected supports. `codes.certify_support_recombination` certifies the resulting abstract support code; a separate check must still establish all spherical layers and the current record. See `research/family_recombination_20261105/` for 32–35D experiments.
`search.support_coupled_recombination` additionally lets a support-family search remove individual accessory supports at a specified cost while optimizing a net score. Its certificate proves support overlap constraints only; the caller must independently establish that the chosen weights and all spherical layers yield a valid kissing configuration.

`search.binary_subcode_cosets` samples binary subcodes, constructs the complete syndrome-distance graph, and searches for more mutually compatible cosets than the original code uses. `codes.certify_binary_subcode_cosets` promotes a successful proposal through independent exact quotient BFS and nonzero-kernel distance DP. The bounded backend now supports up to 19 syndrome bits; search timeouts or skipped large candidate graphs do not prove a code bound. See `research/dense_subcode_cosets_20261105/` for the 32–35D dense-layer target.
`codes.binary_span_cosets` exactly verifies a union of explicit affine binary subspace cosets by enumerating the underlying span and each cross-coset difference, under explicit work budgets. It can certify small-kernel constructions even when the full quotient is too large to enumerate.
