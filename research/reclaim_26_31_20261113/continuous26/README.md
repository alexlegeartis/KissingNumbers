# Continuous and module-transplant pilots in 26–27 dimensions

13 November 2026. **The initial 26–27D pilots found no improvement. A later
30D section construction gives 248,084 points, improving the PackingStar
input by 92 and the preceding 248,020-point campaign result by 64.**
Its exact construction and independent changed-pair replay are preserved in
`ls15_section30/`; a separate complete final-configuration replay passed in
the main verification campaign. The frozen inputs are PackingStar commit
`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`, with 201,180 points in 26D and
209,862 in 27D. Their files are expected under the neighboring `PackingStar`
checkout. The other campaign results are outside this directory's scope.

## Positive follow-up in 30D

The algebraic track constructed an exact 2,564-point 15D Leech–Sloane module
with a 240-point E8 equator. Setting its last tail coordinate to zero gives
a 1,356-point 14D section. All seven choices of a deleted tail coordinate
have that count. Lower-dimensional coordinate sections do not exceed the
incumbent normal-sector counts: their maxima are 308, 348, 552 and 728 in
dimensions 26, 27, 28 and 29 respectively.

The successful 30D section has 240 E8 points, 1,056 heads, and 60 pure tails.
An E8 Weyl transformation and signed tail permutation reduce its exterior
deletion cost to 32. It replaces the 932 incumbent heads and keeps the
shared 300 equator/pure-tail points:

```text
247,992 - 932 - 32 + 1,056 = 248,084.
```

The improvement uses partial head fibers. Its 252 distinct norm-16 tail
vectors have occupancies:

| Heads per tail | Number of tails | Total heads |
|---:|---:|---:|
| 2 | 192 | 384 |
| 10 | 48 | 480 |
| 16 | 12 | 192 |
| Total | 252 | 1,056 |

The twelve coordinate-axis tails carry sixteen heads each. Among the fifteen
four-coordinate supports, three have sixteen tails with ten heads each,
and the other twelve have sixteen tails with two heads each. Thus restricting
every selected fiber to its full size would miss this construction.
All final coordinates are integers before normalization; this improvement
does not require irrational inner products.

`ls15_section30/best_module.npz` stores the exact 14-by-14 transformation
numerator `R` and denominator four. Start with the canonical 15D module
from `../algebraic/ls15_exact/module.npz`, retain rows whose coordinate 14
is zero, drop that coordinate, and calculate `section @ R.T / 4`.
Every division is exact, and `R @ R.T = 16 I`. The omitted tail index is
six, with zero-based indexing among the seven tails.

The minimal norm-128 exchange, after cancelling 300 unchanged points,
removes 964 points and adds 1,056. `finalize_section30.py` independently
rebuilds the final point set and verifies all **260,861,568**
new-to-retained-source products and **557,040** new/new pairs by direct
int64 arithmetic. Every squared norm is 128 and every distinct inner
product is at most 64. `ls15_section30/verification.json` records the exact
source and coordinate hashes. The complete final witness is
`ls15_section30/coordinates30_248084.npz`, and the minimal delta is
`ls15_section30/exchange_minimal.npz`.

A final literal union with the earlier 248,020-point BW14 construction
gave no further gain. The 800 BW-exclusive points have distinct incompatible
blockers among the 864 LS-exclusive points. The exact injection is saved in
`ls15_section30/bw_union_matching.npz`; its scope is only that two-source
union. The bounded orientation search completed 553,282 iterations in
120 seconds, reaching its best exterior cost after about fourteen seconds.

## Main completed experiments

| Experiment | Executed scope | Outcome |
|---|---|---|
| Continuous 26D module deformation | 96 searches for 17 replacements of the 16 loose radicals; eight searches freeing the complete 381-point attempted module | No feasible numerical candidate; best complete-module cosine excess was about 0.0981 |
| Fixed-module insertion polytope | 1,024 supporting-LP starts in ten dimensions | Largest sampled squared radius 24; insertion needs at least 32; no global maximum claim |
| Rational tail angles | 32 off-axis Pythagorean directions with denominator 65, paired with norm-16 E8 heads | 640 new exterior-admissible points; numerical MILP returned 380 on the 1,020-point union |
| Rational latitude variants | 15 rational radius choices and 36 rational tail angles; incumbent blocker cutoff eight | Only 352 unscaled new candidates survive; numerical MILP returned 380 on the 732-point union |
| Rational cavity vertices | 2,048 starts after removing the sixteen radicals | Sixteen sampled vertices of norm squared 40; their radial contractions fit only eight points in that finite candidate pool |
| Move the complete 380-point module to other Golay octads | All 759 octads, one deterministic coordinate ordering each | Every integer-capacity two-copy minimum cut has gain zero; pooled 2,004 distinct cheap candidates are retained for possible joint work |
| Exact 438-point E8/D3 module in 27D | 256 E8/tail signed-axis orientations | Every tested two-copy union has gain zero; minimum full-module exterior cost 8,000 |
| Same 438-point module with general rational tail rotations | 1,024 quaternion rotations plus E8 reorientations | Every tested two-copy union has gain zero; minimum full-module exterior cost 8,384 |

The two-copy conclusions concern only the explicit module and source union
for each tested orientation. They do not cover simultaneous mixtures of
several transplanted modules or arbitrary deformation of the retained source.
MILP solver optima are numerical scheduling evidence, not exact upper-bound
certificates. The point clouds saved by the discovery programs are not new
kissing-number certificates.

## Exact 438-point module

`lambda11_transplant.py:module()` constructs the code in eleven coordinates:

* 240 E8 roots of squared norm 32 in the first eight coordinates;
* 192 points with a norm-16 head and a norm-16 tail;
* six pure tail points `±4 sqrt(2) e_i`.

The heads use three cross-polytopes with sixteen points each. The first is
the coordinate cross-polytope `±4 e_i`. The other two use signed rows of
`2 H4` on the two tetrad partitions
`0123|4567` and `0145|2367`. Their distinct-frame products have absolute
value at most eight. Assign one head frame to each of the three tail
coordinate planes, using all four tails `2 sqrt(2)(±e_i ±e_j)` in that
plane. This yields `3*16*4=192` caps.

Same-tail distinct heads have product at most zero. Different tails in
one plane have product at most zero. Different-plane tails have product
at most eight, and distinct-frame heads have product at most eight.
Each cap therefore has product at most sixteen with every other cap.
E8-root/cap and pure-tail/cap products also have bound sixteen. The program
checks all squared norms and all 95,703 unordered pairs using exact signs
in `Q(sqrt(2))`. Thus the standalone 438-point code is valid. This package
does not require or assert an independent identification with a particular
laminated lattice.

`lambda11_438_module.npz` is the exact small module. It does **not** by
itself improve the 27D lower bound. A nominal replacement of the current
422-point normal sector would gain sixteen only if there were no external
conflicts. The complete exterior pricing is the decisive failed condition.

## Fixed-exterior compression

The normal octad is `{0,2,4,5,6,10,11,23}`. Appending the tail coordinates
gives the ten-dimensional 26D normal module and eleven-dimensional 27D module.
`inspect_normal.py` reconstructs the exact exterior projection multiplicities.
For 26D they are:

| Squared profile norm | Distinct profiles | Multiplicity of each |
|---:|---:|---:|
| 0 | 1 | 4,320 |
| 8 | 240 | 512 |
| 16 | 2,300 | 32 |

The module contains 380 points. Sixteen algebraic points have two contacts,
whose projections to the tangent space are opposite. Their tangent contact
span therefore has rank one. The other 364 module points have numerical
tangent rank nine. Low contact rank motivated the experiments; it does not
imply that an additional point fits.

The octad transplant search uses **all** projected incumbent profiles with
their actual multiplicities. Its bipartite network has a unit-capacity edge
for each new point and an old-profile deletion cost equal to that profile's
multiplicity. Candidate/profile incompatibilities use exact
`a+b sqrt(3)>16` signs in 26D, or `a+b sqrt(2)>16` signs for the 27D module.
Integer maximum flow computes each finite union's best possible gain.

`transplant759.npz` contains the 2,004 collected 26D candidates with at most
four incumbent blockers, and complete old-point blocker indices. They are
all rational: 1,444 have core norm 32, 496 core norm 16, and 64 core norm 24.
Their **combined** optimum was not established here. Individual two-copy
certificates do not automatically exclude joint mixtures.

## Arithmetic correction made during discovery

An early radius-pool generator multiplied large scale factors by an `int16`
head array before assigning to an `int64` output. That intermediate overflow
created invalid points. The erroneous pool and dependent outputs were
discarded. The current generator widens operands before multiplication and
checks every exact squared norm **before** exterior screening. The consumer
also checks the supplied norms before using the pool. The completed radius
report above is from the corrected run. The continuous, cavity, off-axis
denominator-65, and octad-transfer experiments did not use the defective
large-scale arithmetic.

## Reproduction

Python, NumPy, and SciPy suffice. Reserve one available CPU and respect the
repository's combined worker budget. Representative commands from the
repository root (adjust the CPU numbers to the host):

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 4-5 \
  python research/reclaim_26_31_20261113/continuous26/inspect_normal.py --dimension 26

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 4-5 \
  python research/reclaim_26_31_20261113/continuous26/octad_transplant.py \
    --trials 759 --collect --cutoff 4 --out /tmp/transplant759.json

OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 4-5 \
  python research/reclaim_26_31_20261113/continuous26/lambda11_transplant.py \
    --trials 1024 --rotate --out /tmp/lambda11_rotated1024.json
```

The complete octad sweep took about eight minutes on two logical CPUs.
The 1,024 rotated 27D module tests took about fifty seconds. These runtimes
describe the completed finite experiments, not exhaustion of the remaining
research problem.
