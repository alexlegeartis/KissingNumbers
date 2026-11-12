# New tail alphabets and Barnes–Wall module transplants

12 November 2026. This is a bounded research campaign against PackingStar
commit [`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
**This branch of the campaign has not produced a new record.** Its strongest
alternative 31D seed has **251,244 points**, ten below the 251,254-point input.
The explicit replacement, exact changed-pair proof and executable searches are
retained so that subsequent work can start from the changed geometry.

The normal octad is `I = {0,2,4,5,6,10,11,23}`, with zero-based indices.
Unless a file states otherwise, the discussion uses squared norm 32 and maximum
distinct-pair product 16. The complete forced-seed files and deltas use twice
these coordinates, squared norm 128 and product limit 64.

## Results and scope

| Experiment | Explicit scope | Result |
|---|---|---|
| 29D norm-24 core bodies with new norm-8 tails | 1,792 signed six-subsets of I, 256 sampled incumbent bodies, all 160 tails of shape `(±2,±1,±1,±1,±1)` | 327,680 candidates priced; minimum five old blockers. The 4,000-point selected pool remains unresolved. |
| 30D norm-8 E8 bodies with norm-24 tails | 240 E8 roots on I × 64 tails `(±2)^6` | 15,360 candidates. Exact five-blocker certificate and a capacity argument exclude gain from this entire finite family. |
| 31D analogous norm-8/norm-24 latitude | 240 bodies × 448 signed six-coordinate tails | 107,520 candidates, minimum six blockers. The 583 candidates with at most ten blockers gave numerical optimum zero. The full selection problem was not solved. |
| BW15 module orientations | 1,536 sampled exact orientations; 250,000 subsequent annealing steps | Best full-module external deletion cost fell from 114 to 68. No tested two-copy union improved the incumbent. |
| Forced best BW15 replacement | Remove 2,006 actual old points and add 1,996 actual new points | 251,244 points, net −10. All changed inequalities checked exactly. |
| Contact closure of the new BW geometry | All 1,996 newly inserted heads against all 251,244 seed points | 1,281,954 sum incidences and 241,358 distinct new points; cutoff two admits eight and cutoff four admits 428. Both numerical selections return zero gain. |

The source, hashes, selection scope and numerical status are in `reports/`.
No unfinished optimization or sampled family is described as a global ceiling.
The first 64 BW orientation trials were a prefix of the 1,536-trial sweep, not
64 additional distinct trials.

## 1. Why change the tail alphabet?

Earlier midpoint searches paired norm-24 core bodies with norm-8 tails
`(±2,±2)`. When there are at least five tail coordinates, the same norm also
admits `(±2,±1,±1,±1,±1)`. There are 160, 960 and 3,360 such signed tails in
dimensions 29, 30 and 31 respectively. These candidates lie outside the old
two-root-tail grammar.

For a fixed core body, the script groups incumbent points by their exact tail
vector and counts core inner products in an integer histogram. A candidate
tail then obtains its complete deletion cost by looking up the thresholds
`16 - <new_tail, old_tail>`. Small integer products are exactly represented in
binary32; admitted blocker lists are separately reconstructed by direct integer
coordinate sums.

In the 29D pilot, no candidate has at most four blockers. There are 316 at
cutoff eight, 20,910 at cutoff sixteen, and 37,310 at cutoff thirty-two. The
4,000 cheapest candidates have 41,084 blocker incidences, 855 distinct old
blockers and 357,106 mutual conflicts. A maximum matching covers only 841 of
them, so there is substantial deletion sharing. This does **not** imply an
improvement: new/new conflicts must also be respected. The large MILP was
interrupted during presolve and is unresolved.

`tail_shell_pilot.py` also implements the full 6,720-vector E8 norm-24 body
shell, with shapes `(±2)^6`, `(±4,±2,±2)` and
`(±3,±3,±1,±1,±1,±1,±1,±1)` with the appropriate even sign parity.
That expanded option was implemented but was not run in this campaign.

## 2. A complete 30D obstruction for one new latitude

Consider the 15,360 candidates `(h,t)` where h is one of the 240 norm-8 E8
roots on I and `t ∈ {−2,2}^6`, so `||t||² = 24`.
The certificate gives **five distinct incompatible incumbent points for every
candidate**. All these old blockers have tails supported on four coordinates,
with values ±2, and the verifier checks `<h, old_core> ≤ 8` for every h and
every listed old blocker.

For a compatible selected set A, fix a tail t. Two core bodies at that tail
must satisfy `<h,h'> ≤ 16−24 = −8`. Since both have squared norm eight,
there are at most two, and if there are two they are antipodes. An old blocker
can conflict only when its four tail signs agree with t: otherwise its tail
product is at most eight and its total product is at most sixteen. There are
exactly four such t values, from the two remaining free signs. At each t,
at most one selected core body has positive product with that old core,
because two selected bodies there are antipodes.

Consequently one deleted old point can account for at most four selected
candidate/blocker incidences. Every selected candidate has five distinct
listed blockers, so

`5 × added ≤ 4 × deleted`.

No positive exchange can come from this finite family. The statement allows
arbitrary incumbent deletions and individual, non-antipodal choices of new
points. It does not cover new core shapes, different radii, other tail vectors,
or mixtures with other new latitudes.

`verify_latitude30.py` independently regenerates all 15,360 candidates and
checks the 76,800 supplied forbidden products using integers. The minimum
witness product is twenty, above the allowed sixteen. It also checks all
conditions used in the four-incidences argument. The compact certificate is
`data/latitude30_five_blockers.npz`.

## 3. The 31D latitude has a cube incidence bound

The 448 tail vectors with six ±2 entries and one zero are exactly the edge
midpoints of the seven-dimensional cube with vertices `{−2,2}^7`.
At any cube vertex, two distinct incident edge midpoints have product twenty;
two points with the same midpoint have tail product twenty-four.

Thus the core bodies of all selected points incident to one cube vertex,
each of squared norm eight, have mutual products at most −4. If there are m,

`0 ≤ ||sum h_i||² ≤ 8m − 4m(m−1) = 4m(3−m)`,

so m is at most three. Every selected point is incident to two cube vertices.
There are 128 vertices, hence the entire new latitude contains at most
`128 × 3 / 2 = 192` points.

This bound permits arbitrary real norm-8 core bodies; the E8 root catalogue is
not needed. It is a capacity bound for this tail alphabet, **not a bound on
net exchange gain or on the kissing number**. The full 107,520-point E8-body
selection problem remains open in this package.

## 4. BW15 replacement and the 178 external obstructions

The coordinate section of the standard 4,320-point Barnes–Wall configuration
in 16 dimensions has 2,340 points in 15 dimensions. The implementation
regenerates it as 420 roots `(±4,±4)` and 15 affine-hyperplane supports in
`F₂⁴`, each carrying 128 even sign patterns `(±2)^8`. It checks every pair of
the 2,340-point module by integer arithmetic.

Splitting off an eight-coordinate E8 equator gives 240 core points, 84 pure
tail points and 2,016 norm-16/core, norm-16/tail heads. The incumbent module
on I plus the seven tails has 2,282 points: the same 240 and 84, but only
1,958 heads. A complete module swap therefore has nominal gain 58.

Only 178 incumbent points outside this module can obstruct the swap:

- 108 midpoint points whose projection onto I has squared norm sixteen.
- 70 other head points whose projection onto I has squared norm eight.

The scripts check finite integer product bounds covering every other old
point. The actual union conflict graph has 2,136 relevant old vertices:
these 178 exceptions plus the 1,958 old module heads. Each BW orientation
uses exact E8 root reflections, signed permutations, and tetrad isometries;
tail transformations are signed permutations. A maximum bipartite matching
selects an optimal union within that tested pair of configurations.

The 1,536-orientation sweep found no gain. Optimizing the external deletion
count over 50,000 moves reduced it from 114 to 76. A separate 200,000-move
warm restart reduced it to 68. Forcing that full module gives

`251254 − 1958 − 68 + 2016 = 251244`.

Twenty new heads coincide with old heads, so its canonical coordinate-set
delta removes 2,006 and adds 1,996. The preceding 251,236-point alternative
is retained as a distinct geometry. These are research seeds below the
incumbent, not registered improvements.

`force_bw_seed.py` verifies all **506,528,064** source/new products with exact
bounded-integer binary32 arithmetic and all **2,031,120** new/new unordered
pairs with int32 arithmetic. It checks norms, uniqueness and the exact delta.
Unchanged source pairs still require the source's full geometric proof.

## 5. Reproduction

Requirements are Python, NumPy and SciPy. Obtain the pinned PackingStar
coordinate files using the source URLs above; the scripts take explicit
paths and do not download anything. Choose output directories that do not
already exist. Run compute jobs with an affinity that leaves a logical CPU
free, and set BLAS/OpenMP thread counts according to the available machine.

From this directory, with `PACKINGSTAR_SOURCE` pointing to the pinned
`25D-31D/new_bounds_20261111` directory:

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python verify_latitude30.py \
  "$PACKINGSTAR_SOURCE/D30_247992_coordinates.npz"

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python force_bw_seed.py \
  "$PACKINGSTAR_SOURCE/D31_251254_coordinates.npz" \
  data/bw15_module68.npz --out /tmp/bw15-forced-replay

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python contact_refill.py \
  /tmp/bw15-forced-replay/forced31_251244.npz \
  /tmp/bw15-forced-replay/exchange.npz \
  --out /tmp/bw15-contact-replay --cutoff 2 --seconds 60

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python tail_shell_pilot.py \
  "$PACKINGSTAR_SOURCE/D29_231774_coordinates.npz" \
  --out /tmp/new-tail29 --pool-limit 4000

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python bw15_transplant.py \
  "$PACKINGSTAR_SOURCE/D31_251254_coordinates.npz" \
  --out /tmp/bw15-orientations --trials 1536
```

`optimize_bw_orientation.py --help` documents the random seed, temperature,
restart interval and optional initial module. Stored reports record the actual
budgets. `select_price_pool.py` and `solve_tail_pool.py` allow coherent candidate
pool experiments from an existing price table; large MILPs may remain unresolved.

The package keeps small modules, coordinate deltas, the 30D certificate and
compact reports. It omits the large transient priced pools and full forced
configurations, which `force_bw_seed.py` reconstructs from the pinned input.

## Prior-search audit

The previous [competitive 28D campaign](../../discovery_campaign_20261104/competitive28/README.md)
already excludes improvement from its specific 24,178-point multi-sextet pool
by an injective blocker assignment; 24,114 of those points have zero tail.
The prior [whole-layer campaign](../../exceptional_tail_replacement_20261108/README.md)
closes fixed-tail head capacities in 27D and 28D and records the E6/E7 spinor
gate failures. Its final 3,762-point 31D pool has an exact restricted no-gain
certificate. Those statements do not exclude the new radii, tail alphabets
or forced BW geometry explored here.

The historical 6,720-point 28D Weil result was a small restricted construction,
not a miscount of the 221,988-point configuration. Correcting duplicate phase
admission did not make it competitive. This campaign uses the full current
PackingStar coordinates as its starting point.
