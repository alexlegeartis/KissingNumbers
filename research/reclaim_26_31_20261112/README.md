# Reclaiming 26–31D: November source review and construction campaign

**12 November 2026. Two small improvements are established: 222,120 in 28D
and 231,776 in 29D. Substantial improvements across 26–31D remain open.**

The comparison uses PackingStar's release labelled November 11, pinned to
commit [`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
The starting KissingNumbers registry is commit
[`f1be2cc3d2ae2e5d7c0ebab447c3533e43e9c18d`](https://github.com/alexlegeartis/KissingNumbers/tree/f1be2cc3d2ae2e5d7c0ebab447c3533e43e9c18d).
These are exact source identities, not an exhaustive worldwide priority survey.

| Dimension | Repository before this campaign | PackingStar release | Registered repository result after campaign | Difference from PackingStar |
|---:|---:|---:|---:|---:|
|26|201,136|201,180|201,136|−44|
|27|209,825|209,862|209,825|−37|
|28|221,988|222,104|**222,120**|**+16**|
|29|231,697|231,774|**231,776**|**+2**|
|30|247,972|247,992|247,972|−20|
|31|251,200|251,254|251,200|−54|

The full coordinates and proofs for the two improvements are in
[`../../verifications/improved/dim28-29-november-sextets/`](../../verifications/improved/dim28-29-november-sextets/).
The 28D exchange removes 40 PackingStar points and inserts 56; the 29D exchange
removes six and inserts eight. The retained configurations are credited to
PackingStar. The discovery uses the existing corrected Golay sextet union
method, with fresh inputs and explicit new exchanges.

Both final configurations were checked completely: squared norm 128,
distinctness, and signed inner product at most 64. The verification covers
all 24,668,536,140 unordered pairs in 28D and 26,859,941,200 in 29D, using
actual antipodal pairs only where both points occur. All unpaired points are
included. Bounded integer products and every accumulation are exactly
representable in binary32; no geometric tolerance is used. A separate direct
int64 path rechecks 12,438,664 and 1,854,200 changed-point comparisons and
reconstructs the pinned source hashes. Four malformed witnesses are rejected.

## What the larger experiments established

The main obstacle is compatibility with the incumbent's nonzero-tail layers.
Compatibility with the Leech core alone substantially overstates the net gain.
Many individually expensive candidates share old blockers, so individual
deletion prices also understate the potential of a coordinated replacement.
Both effects must be modeled together.

| Experiment | Outcome | Scope |
|---|---|---|
| Fresh sextet unions | +16 in 28D; +2 in 29D | Exact complete final witnesses; finite discovery frames and candidate cutoff |
| Literal union of our previous result and PackingStar | No improvement in all six dimensions | Exact injective blocker witnesses for every repository-exclusive point |
| Classical scaled-Leech free completion | No new point in all six dimensions | Complete fixed-source family with existing D_k tails; deletions not included |
| Amplify the 56 new 28D points onto all target D_k tails | No improvement in all six dimensions | All 7,104 candidates have distinct verified incumbent blockers; no blocker cutoff |
| Explicit three-core/five-tail layer | 2,688 mutually compatible points | Best sampled orientation deletes 9,649 old points; whole-fiber subsets in that orientation have exact maximum gain zero |
| New norm-20/norm-12 radius | 501,760 candidates; minimum 16 blockers | Selected 6,000-point MILP timed out without an improving incumbent; unresolved |
| Quarter-angle midpoint tetrahedra | 143,360 four-point bundles; minimum 65 deletions | Sampled bodies, all coordinate tetrahedra; no joint selection claim |
| Shrunk midpoint triangles | Minimum 334 deletions per bundle in the 31D pilot | Fixed planes and phases; exact 280-head bound per phase |
| Norm-8/norm-24 latitude | Exact no-gain certificate in 30D | Entire 15,360-point finite family; a distinct 31D family has capacity at most 192 |
| Whole signed-octad redesign | 928 compatible supports, 218,948 points in 31D | Far below the incumbent; conditional targets are not achieved results |
| Barnes–Wall module replacement | 251,244-point alternative 31D seed | Net −10 versus PackingStar; changed-pair proof inherits unchanged source validity |
| Refill the Barnes–Wall seed | No free Leech/classical points; no positive contact selection at cutoffs two or four | Finite tested closure; no maximality assertion |
| New large-tail coherent fibers | 102,400 fibers priced; selected 4,000-fiber model unresolved | Whole fibers on one normal octad; no optimum or gain claimed |
| Individual signs on mixed core supports | 726,528 points priced; no positive exchange in a 19,110-point adaptive pool | Three bounded heuristic schedules; unresolved, not an exclusion certificate |

The 251,244-point Barnes–Wall seed is useful changed geometry, but it is not a
new lower-bound record and is not substituted into the registry. Forcing a
2,340-point module instead of the old 2,282-point module offers a nominal gain
of 58, but the best tested orientation requires 68 additional external
deletions. All newly affected products pass exact checks. Its later contact
closure generated 241,358 distinct candidates from 1,281,954 incidences; the
428 candidates at blocker cutoff four gave numerical maximum gain zero.

The later mixed-support individual-point pilot reached minimum deletion costs
of four for the three-core/five-tail family and two for the two-core/six-tail
family, substantially below the whole-fiber costs. However, the selected
19,110-point pool has 4,105,565 mutual conflict edges. Neither ordinary nor
hotter annealing nor a continuation of the deletion penalty found a positive
exchange. The best saved nonempty exchange had gain −1; the empty exchange
preserves the source. All saved witnesses were checked with direct integer
products. This pilot makes the broader search reproducible without implying
that a thousand-point gain is near.

## Package map

| Directory | Contents |
|---|---|
| [`source_and_layers/`](source_and_layers/) | Current-source schemas and hashes, exact six-source union obstructions, complete classical residual scans, successful sextet discovery reports, large-tail pricing and optimization pilots |
| [`modules/`](modules/) | New tail alphabets, exact 30D incidence certificate, 31D cube capacity argument, Barnes–Wall module search, forced-seed deltas and contact closure |
| [`midpoint_radii/`](midpoint_radii/) | Five-plus-three support pilot, quarter-angle tetrahedra, exact radical 26D sextet search and Leech/classical refills |
| [`theory/`](theory/) | Explicit 2,688-point layer, signed-support redesign, exact triangle bound and complete exchange-amplification matching certificates |
| [`proof_review/`](proof_review/) | Independent source/delta and synthetic brute-force coverage review of the new 28D/29D proof package |
| [`STRATEGY.md`](STRATEGY.md) | Concrete larger-gain formulations, candidate families, acceptance gates and limitations |

Search programs and reports are separated from registered mathematical
certificates. Solver timeouts, forced seeds below the incumbent, numerical
optima and valid auxiliary layers are not promoted into lower-bound claims.
Read each subpackage's README for its exact scope and source requirements.

## Reproduction and resources

The complete positive certificate needs only Python and NumPy:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 \
  python verifications/improved/dim28-29-november-sextets/verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0-1 \
  python verifications/improved/dim28-29-november-sextets/negative_controls.py
```

Run these commands from the repository root. Short research certificates have
separate read-only replay commands in their READMEs; they require NumPy and
exact source hashes where applicable. Discovery uses SciPy or a C++ compiler
in some subpackages. Large transient candidate matrices and complete copies
of external coordinates are omitted when they can be regenerated from the
pinned source and the supplied code.

`python research/reclaim_26_31_20261112/verify_manifest.py` checks the complete
campaign file inventory and SHA-256 hashes. It establishes file identity,
not the mathematical conclusions; replay the relevant certificate programs
for those conclusions.

The execution host had nine available logical CPUs and about 9.7 GB memory.
The campaign reserved one CPU and bounded the combined active worker count.
No GPU or separate large-memory research server was available in this session.
This resource statement describes the experiments actually performed; it does
not imply that an unresolved large optimization has been exhausted.

The registered manuscript remains at its documented earlier snapshot. The
current repository registry, source comparison and proof runner include the
two new results. Historical certificates retain their construction counts;
current-registry columns are refreshed explicitly where needed.
