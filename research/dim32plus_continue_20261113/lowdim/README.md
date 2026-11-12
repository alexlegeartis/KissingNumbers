# Bounded low-dimensional continuation, 13 November 2026

## Verified result

This campaign produced **501,822 points in dimension 36** and **505,638 in
dimension 37**, increasing each preceding repository witness by **nine points**.
The final proof package is
[`dim36-37-orbit-exchange`](../../../verifications/improved/dim36-37-orbit-exchange/).
Its default verifier freshly replays the pinned preceding construction, rebuilds
the entire current base from exact short-shell enumeration, and directly checks
every affected pair. The full replay passed in 48.1 seconds on one assigned CPU.

The compact exchange removes the three phase orbits (23,77,23,u,1), u=0,1,2,
and adds (23,105,23,u,1) and (38,77,12,u,2), again for all three u. It replaces
nine points by eighteen. All 90 existing base-triangle switches and all 636
latitude-body orbits survive. The complete 37D retention check is part of the
proof, rather than an assumed consequence of the 36D count.

The package's exact new checks cover 7,266,672 added/base products, 1,766,106
added/final-equator products, and 34,344 added/latitude-body products. Their
maxima are 6, 6, and 5, against the allowed product six at squared norm twelve.
Five corruption controls pass. Restoring a deleted orbit gives an exact bad
product of ten despite retaining the correct norm and updating the claimed count.

`delta36.json` and `delta36_verification.json` retain the discovery-stage compact
delta and its direct coordinate check. The final package additionally pins all
dependencies and freshly replays the complete unchanged baseline proof.

## Searches and bounded negative outcomes

### Current partially switched Coxeter–Todd base

The admission preparation uses the actual current 90 individual triangle
switches. Each base triangle contains all three root phases, so exact admission
of a candidate phase orbit reduces to pairwise label tables even when only some
large vectors of an owner were switched. The script freshly enumerates the short
shells, constructs these tables, and checks that the incumbent is admitted.

There are 722,728 admitted product cells. Preserving the 636 latitude orbits
leaves **6,444,168** candidate equatorial phase orbits.

| Search | Candidate vertices | Conflict edges | Bounded work | Best result |
|---|---:|---:|---|---:|
| At most six incumbent blockers | 103,903 orbits | 2,026,072 | 100,000,000 moves, at most 600 seconds | 32,706 orbits |
| Shared blocker sets of size at most four | Same six-blocker pool | Recomputed local conflicts | 13,295 candidate regions in one pass | 32,706 orbits |
| At most eight incumbent blockers | 178,548 orbits | 10,963,489 | 78,972,928 moves in 600 seconds | 32,706 orbits |
| Individual points, with order-three symmetry released | 535,644 points | 32,890,467 | 17,832,960 moves in 600 seconds | 98,118 points |

The first improvement appeared after 151 moves in the six-blocker search.
Discarding neutral conflict components reduced the stored exchange from 45
removals/48 additions to the three removals/six additions above. The later
searches did not improve it. These are bounded heuristic outcomes, not upper
bounds on any full construction family. The shared-blocker routine also limits
large local subproblems and search nodes; its negative output is not an exhaustive
independent-set certificate.

### Additional 37D latitude search

The new equator was used to rebuild the lift pool against the actual switched
base. This produced **168,344** admissible three-point body orbits of norm pattern
(4,4,0). No orbit could be added by itself. The pool of candidates with at most
four blockers contains 891 orbits, including the retained 636, and 8,265 conflict
edges. A 120-second search made 13,324,288 moves and retained a best of 636.

A separate 60-second numerical MILP call asked for at least 637 orbits in that
891-orbit pool. HiGHS reported infeasibility. This is recorded as a restricted
numerical search outcome, not a formal mathematical upper-bound certificate.
The 37D result remains **505,638**.

`screen_other37_latitudes.py` establishes exact obstructions for several further
integer-body shapes while retaining the current base and latitude layer:

- All **18,289,152** bodies of pattern (4,6,0), in every block order, conflict
  with an unchanged (6,6,0) base point. Their minimum forced product is seven,
  above the limit six. The 224 norm-six owners in each block pair cover every
  norm-six coset, so this includes each individual body, not only whole cells.
- All **61,236** single-block norm-eight bodies have product eight with a retained
  base point. Every triangle switch preserves its large norm-eight block.
- Single-block norm-four and norm-six bodies at the corresponding larger
  heights conflict with a retained latitude point. The smallest possible maxima
  against the occupied root directions are respectively two and three.

These obstructions do not cover arbitrary real bodies, changed base ownership,
different latitude systems, or unexamined short-shell patterns.

### Dense fourteen-coordinate heads in 38D

The pilot used full fourteen-coordinate sign heads at body norm fraction 1/4
and the two unused orthogonal body lines 4(e22 ± e23), with zero-based indices.
These points satisfy the retained Fano-grid and Hamming cross-layer bounds.
The cap and flat-axis requirements are the obstruction.

All 8,192 antipodal head lines fail the uniform cap/axis screen. Using the actual
Leech cap-class assignments instead, the cheapest head on the two tested body
lines forces respectively **124** or **114** additional cap-point deletions,
before deleting its own cap line. Each head also conflicts with 12 or 14 flat
axes, while one head/body pair contributes four new points. This did not yield
a candidate. Shared deletions among larger jointly modified families were not
exhausted, and no general 38D ceiling is claimed.

No 38D or 39D improvement was obtained. Their current constructions use different
geometries, so the 36D exchange was not transferred by an unsupported inheritance
or monotonicity assumption.

## Independent extended-normal review

The fully independent Gamma72 review is packaged in
[`dim68-69-gamma72-extended-moments`](../../../verifications/improved/dim68-69-gamma72-extended-moments/).
Its `verify_independent.py`, `independent_verification.json`, and
`INDEPENDENT_REVIEW.md` use only package-local data. The exact reviewed bounds
are **471,886,208 in 68D** and **746,325,051 in 69D**.

The reviewer directly generated all 162,969 possible D7 profiles, reconstructed
the 1,681-profile D4 action by explicit generator traversal, and directly summed
all 18 moments and objective weights. All 60 outer orbit coefficients, 84 known
roots, rational dual inequalities, and geometric pair types agree. Both production
`--enumerate` replays passed. Sixteen adversarial controls passed against the
research precursor, and all thirteen controls passed against the final package.

This proof retains the stated published Gamma72 identity, integrality, minimum,
shell-population and spherical 11-design premises. The review does not enumerate
the ambient minimal shell. The additional three normal coordinates strengthen
the population constraints and remain inside the 68D body.

## Independent 67D affine-sign review

The middle-dimensional lane produced a frozen joint extension from 355,090,822
to **355,095,446 points in dimension 67**, adding **4,624** points in 49 proper
affine sign classes. A fully independent reviewer was written here and then
adapted into
[`dim61-67-compatible-sign-families`](../../../verifications/improved/dim61-67-compatible-sign-families/)
as `verify_subcode67_independent.py`. The final mathematical input is the
package's `data/subcode67.json`, with SHA-256
`90a421e678dd2f61f9fd668c374152145a6937dcfa4cbd4f5107825a44b5fe93`.

The independent method filters all 65,536 possible binary words against each
listed parity check. It does not import the discovery code, production verifier,
row-reduction routine, or shared-parity proof. It verifies the full sizes and
minimum distances of all changed classes, reconstructs all 40,428 shortened
parent supports, and directly multiplies all affine sign pairs whenever support
intersection exceeds eight.

The exact review passed:

- All 357 replacement classes retain 2,048 points and minimum sign distance four.
- The 49 new classes have sizes 512 (one), 256 (six), 128 (seven), 64 (nineteen),
  32 (thirteen), and 16 (three). Their minimum distances are four in 45 classes
  and six in four classes.
- All 15,072,183 changed-old support interactions have intersection at most eight.
- All 2,069,956 new-support interactions are covered. The 758 exceptional pairs
  have intersection nine in 238 cases and ten in 520 cases. Direct enumeration
  checks **115,998,720** sign-pair inner products, with maximum **32** at squared
  norm **64**. The other new interactions are certified by support overlap.

The retained dense, four-support, and axis layers admit all new and reordered
classes by their coordinate-amplitude bounds. Distinct supports and distinct
affine words prove population uniqueness. The independent extension audit
inherits the unchanged parent construction; the final package's main verifier
separately replays the baseline and all parent support pairs.

Eight package-local rejection controls include incorrect population and source
pins, a distance-two sign family, a full class reduced to half size, duplicated
support, and phase changes that preserve all individual norms and populations.
The two phase controls each exhibit an exact bad inner product **40 > 32**.
The same geometric collision is rejected under optimized Python. The default
reviewer is read-only; `--report` is required to write a receipt. Its pure
`review67(path)` function is available to the package's main verifier.

## Script roles

| Script | Role |
|---|---|
| `prepare_current36.py` | Fresh exact base-admission tables and complete equator/lift candidate pools; optional new-equator delta. |
| `filter_phase.cpp` | Bounded incumbent-blocker filtering and immediate insertions; adapted from the retained completion search. |
| `refill_phase.cpp` | Deterministic-step phase-orbit exchanges with a wall-time limit; adapted from the retained completion search. |
| `shared_trades.cpp` | Simultaneous exchanges through shared blocker sets of size at most four. |
| `refill_points.cpp` | Individual-point search that releases the common order-three symmetry. |
| `verify_delta36.py` | Discovery-stage direct integer-coordinate proof of every affected equator-exchange pair. |
| `lift_milp.py` | Bounded numerical MILP discovery on the exact finite latitude conflict graph. |
| `screen_other37_latitudes.py` | Exact rejection of the specified additional-height integer-body shapes. |
| `screen_dense14.py` | Exhaustive dense-head screening and actual cap/axis deletion accounting for the two tested 38D body lines. |
| `review_d7.py` | Independent explicit-profile review of the research Gamma72 certificates. |
| `gamma_rejection_controls.py` | Sixteen malformed-input and geometric controls for the research Gamma72 checker. |
| `review_subcode67.py` | Discovery-stage independent exhaustive sign and affected-pair audit for the frozen 67D joint extension; the final package copy adds eight controls and deterministic receipts. |

## Reproduction and retained data

The complete positive proof is reproduced by the final verification package.
For optional search reproduction, compile the C++ sources with a C++17 compiler
and `-O3`, run `prepare_current36.py`, and invoke the filter/search programs from
the generated pool directory. The exact recorded commands and outcomes are
represented by `current36/pool.json` and its `filter*`/`*search` logs.
`--mode lifts --equator-delta delta36.json --output lifts37` prepares the final
latitude continuation. `lifts37/` contains its pool metadata, search logs and
numerical MILP report.

The campaign used CPUs 4–5 with BLAS/OpenMP/MKL thread counts set to one and
reserved CPU 8. Use affinity appropriate to the host and leave a CPU free.
Time-limited trajectories are not proof inputs.

Large `.bin` candidate pools, selected-set exports, compiled executables,
bytecode caches and redundant NumPy diagnostic arrays are ignored. All mathematical
certificates and positive proof receipts are retained; the search pools can be
regenerated from the scripts and pinned repository inputs.
