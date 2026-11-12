# Additional midpoint radii and exact 26D sextet search

These experiments did **not** produce a new kissing-number lower bound.
They use the pinned November 11 PackingStar coordinates described in
[`../source_and_layers/`](../source_and_layers/), or the explicitly reconstructed
Barnes–Wall research seeds in [`../modules/`](../modules/).
All counts below concern finite, stated candidate families. A timeout is not
a proof of no improvement.

## Norm-20 core, norm-12 tail: signed five-plus-three supports

In squared-norm-32 coordinates, combine a core vector with five entries ±2
and a tail vector with three entries ±2. The 31D pilot used all 1,792 signed
five-subsets of the normal octad and all 280 signed three-coordinate tails:
**501,760 points**. The input contained no existing norm-20 core bodies, so
the requested optional sample from that source layer contributed no extra bodies.

Every old point was eligible for deletion. Exact tail-group histograms priced
every candidate; the minimum was **16 old blockers**. The 6,000 cheapest
candidates had 128,552 old/new conflict incidences, 2,268 distinct old blockers
and 310,096 new/new conflicts. A matching covered only 2,256 candidates,
showing why a per-point deletion count alone does not settle the joint problem.

The shared-deletion MILP returned no improving incumbent within its 150-second
solver limit. Its numerical remaining gain upper bound was 177; that is an
unresolved search statistic, not an exact mathematical bound. The full compact
price distribution and numerical result are in `data/five_support_report.json`
and `data/five_support_optimization.json`.

This radius had already appeared in earlier seed-restricted searches. The
new scope here is the current PackingStar source and a larger shared-deletion
pool; it is not credited as discovery of the five-plus-three alphabet.

## Quarter-angle midpoints with tetrahedral fibers

For Leech vectors u and v with squared norm 32 and inner product eight,
the midpoint h=(u+v)/2 has squared norm twenty. Four tails of squared norm
twelve with mutual product minus four give four valid points above that body.
For every other Leech vector w, both w·u and w·v are at most sixteen, so
w·h is at most sixteen. Only the two parent Leech points can obstruct this
four-for-two replacement within the Leech core.

Compatibility with the incumbent's nonzero-tail layers is the difficult part.
The pilot sampled 2,048 distinct midpoint bodies, half from pairs of missing
Leech points and half from the full shell. All 70 coordinate tetrahedra were
priced for each body: **143,360 four-point bundles**. The cheapest bundle
required **65 incumbent deletions**; the median required 354. The 128 cheapest
bundle prices were independently replayed with direct integer coordinate sums.
No cross-bundle selection or improved final configuration is claimed.

The complete statistical report, sampled parent labels and direct blocker
lists are in `data/quarter_midpoint_report.json`. The large temporary prices
can be regenerated with the fixed seed in `quarter_midpoint_tetra.py`.

## Exact radical sextet continuation in 26D

`sextet_quadratic26.py` represents all points as A+B√3 with integer A and B
at squared norm 128. It transports both parts through the same corrected
sextet isometry, compares the signs of a+b√3 with integer squares, and solves
the admitted bipartite union graph. The **35 completed frames** at blocker
cutoff 32 produced no gain over 201,180. The result excludes only profitable
unions in those admitted pools. It does not exclude larger blocker cutoffs,
other frames, new candidates or joint combinations of several images.

The frame-by-frame transcript is `data/quadratic26_sextet_report.json`.

## Refills after forcing a Barnes–Wall module

The alternative 31D seeds have 251,236 and 251,244 points, below the 251,254
PackingStar source. Their 47,872 missing Leech core points produced **no free
insertions**. For the better seed, 48,943 distinct nonzero-tail core projections
were checked; early rejection evaluated 804,008,448 exact integer products.
The complete classical `(√3/2 L,t)` family with D7 tails also produced zero
new compatible points for both seeds. The saved logs and reports record these
restricted scans, not maximality of either seed. The separate module package
records the subsequent negative-contact closure.

`leech_refill.py` regenerates the complete Leech shell and independently checks
any surviving insertion. `classical_refill.py` uses a supplied regenerated shell,
verifies that the equator is contained in it, and uses exact integer-square-root
thresholds. Their conclusions inherit the stated finite-family geometry;
neither is a full proof of its input seed.

## Reproduction

Use Python, NumPy and SciPy, and reserve at least one logical CPU. The examples
use two worker threads; adapt the affinity to the machine. From this directory,
set `PACKINGSTAR_SOURCE` to the pinned repository's
`25D-31D/new_bounds_20261111` directory. Use fresh output paths.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python export_leech.py /tmp/reclaim-leech.npy

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python five_support_pilot.py \
  "$PACKINGSTAR_SOURCE/D31_251254_coordinates.npz" --out /tmp/reclaim-five-support

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 RESEARCH_THREADS=2 timeout 190s taskset -c 0-1 \
  python solve_five_support.py /tmp/reclaim-five-support \
  --source "$PACKINGSTAR_SOURCE/D31_251254_coordinates.npz" --seconds 150

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python quarter_midpoint_tetra.py \
  --source "$PACKINGSTAR_SOURCE/D31_251254_coordinates.npz" \
  --leech /tmp/reclaim-leech.npy --out /tmp/reclaim-quarter-midpoints

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python sextet_quadratic26.py \
  --source "$PACKINGSTAR_SOURCE/D26_201180_coordinates.npz" \
  --frames 35 --cutoff 32 --out /tmp/reclaim-radical-sextets
```

After reconstructing a forced Barnes–Wall seed using `../modules/force_bw_seed.py`:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python leech_refill.py \
  /tmp/bw15-forced-replay/forced31_251244.npz --out /tmp/reclaim-leech-refill

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 python classical_refill.py \
  --dim 31 --source /tmp/bw15-forced-replay/forced31_251244.npz \
  --leech /tmp/reclaim-leech.npy --output /tmp/reclaim-classical-refill.npz
```

The external `timeout` is deliberate: a solver's internal time limit may not
interrupt every presolve operation. A terminated run has no mathematical
conclusion. The optimizer reproducer imports the existing shared-deletion
model in `research/exceptional_tail_replacement_20261108/optimize_rational.py`;
it checks any proposed positive changed pairs independently before writing
coordinates. The stored numerical pilot result remains the result of the
original executed run, not a promise of identical solver timing on every host.

`data/manifest.json` pins the files in this subdirectory. Recorded absolute
paths inside reports identify the original execution environment; every
executable entry point accepts portable paths.
