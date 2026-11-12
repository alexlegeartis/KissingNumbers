# Source audit and new-layer experiments, 12 November 2026

The current PackingStar source is pinned to
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6).
Its README labels the release **11 November 2026**. The commit timestamp is
10 November 2026, 19:32 UTC; the release label and commit timestamp are reported separately.
The six coordinate files are in `25D-31D/new_bounds_20261111/`.

All source file SHA256 values are in `reports/source_audit.json`.
The older repository certificates compared here are pinned by file hash and
their original commit in `reports/repository_inputs.json`. Later repository
improvements do not change the meaning of these comparisons.

## Actual current structures

Scale each vector to squared norm 32. The columns below classify the squared
norm of its first 24 coordinates; the remaining coordinates are called tails.

|Dimension|Total|Core norm 32|Core norm 24|Core norm 16|Pure tail|
|---:|---:|---:|---:|---:|---:|
|26|201180|192136|8832 +16 quadratic|192|4|
|27|209862|183538|26024|288|12|
|28|222104|172478|48450|1152|24|
|29|231774|163648|66006|2080|40|
|30|247992|150240|92792|4900|60|
|31|251254|148688|95852|6630|84|

The norm 24-core layers use tail vectors `2(±e_i±e_j)`, the D_k roots with
`k=dimension-24`. The norm 16-core layers use norm 16 tails: coordinate axes
`±4e_i` and signed four-coordinate vectors with entries ±2. In particular,
the latest 31D source uses D7 tails; the older paper's E7 triangle construction
does not describe its current midpoint layer.

The 26D archive contains 201,164 integer points and 16 points over Q(sqrt(3)).
All 3,218,624 radical/rational pairs and 120 radical/radical pairs pass an exact
integer sign comparison. This source inspection checks every norm, but is
**not** a rational/rational all-pairs certificate. The ordinary public
PackingStar verifier uses a floating tolerance and does not directly read
this newer NPZ schema; it is not the proof used here.

## Six literal source unions: exact no-gain certificates

For each dimension, let X be PackingStar's source and Y the older repository
certificate. Every point in Y\X has an explicitly stored, distinct blocker
in X with inner product above the kissing threshold. Hence selecting any m
points from Y\X forces deletion of at least m points of X:

`|Z| <= |X|` for every feasible `Z subset X union Y`.

This is a literal finite-union statement. It does not exclude other images,
new coordinates, larger pools, or a global improvement. It does not depend
on a numerical optimization tolerance. `source_union.py --verify-only`
checks injection, complete coverage, input hashes, norms and every matched
angle directly with integers (and Q(sqrt(3)) signs in 26D).

|Dimension|Repository-exclusive points|Matched distinct blockers|
|---:|---:|---:|
|26|156|156|
|27|37757|37757|
|28|2262|2262|
|29|1872|1872|
|30|5188|5188|
|31|116336|116336|

The small index witnesses and complete reports are in `matching/`; source
coordinate arrays are not duplicated here.

## Complete classical-lift residual scans

The probe enumerates every missing Leech vector L and every existing D_k
tail t, testing `(sqrt(3)/2 L,t)`. A present Leech point L itself blocks its
lift, while all other retained Leech vectors automatically clear it.

Against a rational point `(a,b)/den`, exact compatibility is

`L·a <= floor(2(16*den-t·b)/sqrt(3))`.

The nonnegative right side is computed with integer square roots, without
tolerance. The 26D radical points receive a separate exact check. The source
Leech shell is rebuilt from the repository's Golay generator, and every
actual core point must belong to that shell.

|Dimension|Missing Leech vectors|Tail directions|Candidates tested|New compatible points|
|---:|---:|---:|---:|---:|
|26|4424|4|17696|0|
|27|13022|12|156264|0|
|28|24082|24|577968|0|
|29|32912|40|1316480|0|
|30|46320|60|2779200|0|
|31|47872|84|4021248|0|

In 26D the rational part permits exactly the 16 radical points already in
PackingStar's archive. These negative outcomes concern only free insertions
in this fixed classical family; deletions and changed radii remain allowed
in the broader research problem.

## Positive fresh sextet searches

Fresh searches from the November source gave 222,120 in 28D (+16) and 231,776
in 29D (+2). The 28D final difference is 56 additions and 40 deletions;
29D uses 8 additions and 6 deletions. These are small improvements, not the
requested thousand-point-scale advance.

The complete positive proofs and coordinates are maintained in
[`dim28-29-november-sextets`](../../../verifications/improved/dim28-29-november-sextets/).
They are deliberately not duplicated in this research directory.
The 28D run completed 35 sextet frames. The 29D and 30D runs completed 8 frames
each; 30D had no gain in those 8 frames. All used blocker cutoff 32, so a
negative sweep is restricted to that cutoff. The reports preserve every
completed frame and gain.

## Larger new-tail radii: exact setup, unresolved optimization

The coherent-fiber pilot introduces signed ±2 octads split between core and
tail in the ratios `(3,5)`, `(2,6)` and `(1,7)`. They clear every Leech core
point and every pure D7 root. They also clear every flat `(±1)^24` midpoint
body, whose complete dot product is at most 14 in norm 32 scaling.

At a fixed tail, core vectors of squared norm 4r must have pairwise products
at most `4r-16`. Positivity of the squared sum bounds a fiber by
`floor(4/(4-r))`: respectively 4,2,1 points for r=3,2,1. Tetrahedra, antipodal
pairs and single points attain those capacities individually. In 31D the
sum of the three nominal capacities is

`4*32*C(7,5) +2*64*C(7,6) +128*C(7,7) =3712`.

This capacity calculation does **not** assert simultaneous compatibility
or a net gain after incumbent deletions.

The normal-octad pilot uses core coordinates `{0,2,4,5,6,10,11,23}` and all
tails of these three types. It generates 102,400 coherent fibers containing
353,280 individual candidate points. All 80,210 relevant incumbent points
compress exactly to 4,586 projections on these 8 core coordinates and 7 tails.
Multiplicities are retained as deletion weights: if a candidate blocks one
such profile, every incumbent with that projection must be deleted.

All fibers were priced. At deletion cap 96, 99,147 remain. Minimum fiber costs
are 53 deletions for 4 new 3+5 points, 12 for 2 new 2+6 points, and 21 for 1 new 1+7
point. Thus the nominal capacity is not accessible through cheap individual
fibers. The joint model permits shared deletions and includes all cross-fiber
incompatibilities. A 4,000-fiber budget preserves tail-direction coverage;
the selected model has 2,045 deletion profiles and 348,307 incompatibility edges.

**Optimization remains unresolved.** HiGHS returned no result after roughly
300 seconds despite a 90-second requested solver limit, and the run was
interrupted. No improving incumbent was returned. This is neither an exact
no-gain certificate nor a new lower bound. The packaged reproducer isolates
the optimizer in a subprocess and enforces an external wall deadline.

Two restrictions must be relaxed before rejecting this direction: this pilot
keeps each 4/2/1-point fiber whole, and limits its core supports to one octad.
Mixed core triples such as `{11,23,j}` can have cheaper individual points.

## Individual signs beyond the normal octad

A further bounded pilot removed the whole-fiber restriction. It priced
726,528 individual points: all signed 3+5 points on 42 mixed core triples,
all signed 2+6 points on the 276 core pairs, and all signed 1+7 points on
the 24 singleton supports. The triple pool contains every `{11,23,j}` and
20 additional mixed triples selected using the complete-sign-bundle prices.
Every appropriate signed tail direction is present.

Individual minimum deletion costs were 4 for 3+5, 2 for 2+6, and 21 for 1+7.
Their respective medians were 32, 34, and 77. After preserving alternatives
for every signed tail and repricing by blocker reuse, the adaptive pool had
19,110 points, 5,735 possible incumbent victims, 302,518 old/new conflict
incidences, and 4,105,565 new/new conflict edges. Pricing and graph construction
took about 12 seconds on the development executor.

Three bounded heuristic schedules returned no positive exchange. Ordinary
annealing reached at most 14 selected points with 184 deletions, hotter
annealing 46 with 624, and continuation of the deletion penalty from 0.02
toward 1 reached 262 with 1,667. These maxima are trajectory diagnostics;
the best saved nonempty exchange still had gain -1. Every stored best witness
was replayed directly with integer dot products against the full source.
The result remains **unresolved**, not a no-gain certificate. This short
search did not reach a competitive large replacement.

The unrestricted 24-core-coordinate grammar has 11,381,760 points and
473,381,256 incumbent conflict incidences. The counts follow exact
signed-subset generating functions on the source coordinate shapes. A single
uint32 adjacency direction occupies about 1.89 GB before offsets. Full
pricing is plausible; a monolithic model of all candidate pairs is not.
Further work should grow incidence neighborhoods and candidate columns
jointly, with explicit partial signs, signed-tail capacity constraints, and
exact deletion sharing. The present pool is too narrow to reject that route.

## Reproduction

Use the repository's CPU reservation rule. The example CPU range below is
for the 9-CPU development executor and reserves CPU 8; adapt it on other hosts.
Install NumPy and SciPy, and clone PackingStar adjacent to KissingNumbers
or supply `--packingstar /path/to/PackingStar`. Input hashes are mandatory.
Run from this directory:

```sh
taskset -c 0-7 env OPENBLAS_NUM_THREADS=2 python inspect_sources.py
taskset -c 0-7 env OPENBLAS_NUM_THREADS=2 python source_union.py --output matching --verify-only
taskset -c 0-7 env OPENBLAS_NUM_THREADS=2 python classical_residual_scan.py --dim 31 --output /tmp/reclaim-residual31
taskset -c 0-7 env OPENBLAS_NUM_THREADS=2 python sextet_sweep.py --dim 28 --frames 35 --output /tmp/reclaim-sextet28
taskset -c 0-7 env OPENBLAS_NUM_THREADS=2 python large_tail_pilot.py --core octad --cap 96 --max-groups 4000 --solver-seconds 90 --wall-seconds 120 --output /tmp/reclaim-large-tail
taskset -c 0-1 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 timeout 300s python partial_sign_pilot.py --output /tmp/reclaim-partial-signs --report /tmp/reclaim-partial-signs/report.json --wall-seconds 285
```

To rebuild matching witnesses instead of checking them, omit `--verify-only`
and choose a separate output directory. The matching verifier itself checks
only O(N) matched pairs; rebuilding the full 31D cross-conflict graph is more
expensive. Large candidate pools and downloaded sources are not committed.
The partial-sign pilot also requires a C++17 compiler. Its default heuristic
uses continuation; `--search-mode anneal` selects the cycling-temperature
schedule. Wall-clock limits can change proposal counts between machines.
