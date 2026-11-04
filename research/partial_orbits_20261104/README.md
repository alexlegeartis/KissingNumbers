# Coupled partial-orbit and rational-image continuation

This post-merge campaign adds a reusable discovery operation and searches
17D, 27D, 29D, 30D and 31D. **It has no further improved lower bound.** The
certified 27D/31D/43D improvements already merged into main remain unchanged.
No 12D search was run in this continuation.

## Individual phases and changing all layers

The registered [`search.partial_sign_orbits`](../../discovery/docs/PARTIAL_ORBITS.md)
operator selects individual points from diagonal sign-group orbits. It retains
internally conflicting orbits, deduplicates stabilizers, and uses exact signed
relative-phase kernels. Every incumbent point is removable; shared deletions
across height, equator, root and flat layers are charged once. Valid shared-blocker
clique inequalities strengthen the numerical relaxation. A deterministic point
reservoir can explore different partial-orbit frontiers under a fixed pool budget.

The 17D seed is reconstructed from the archived binary code and coset choices,
then independently checked through the engine's full coordinate verifier.
There are 5,730 points, diagonal metric `(1,...,1,2)`, squared norm 288 and
inner-product limit 144. The candidate grammar includes both height signs:

| Body in the first 16 coordinates | Height | Number of points |
| --- | --- | ---: |
| Six entries of magnitude 6 | ±6 | 1,025,024 |
| One entry of magnitude 9 and fifteen of magnitude 3 | ±6 | 2,097,152 |
| Sixteen entries of magnitude 4 | ±4 | 131,072 |
| One entry of magnitude 12 and two of magnitude 6 | ±6 | 26,880 |
| Eight entries of magnitude 6 | 0 | 3,294,720 |

All support and sign choices in these shapes are generated: 6,574,848 points
before removing source overlap. The last two shapes expand the earlier tier-only
catalogue. Source symmetry reduces screening to 211,300 new orbit representatives
and **1,210,749,000** exact representative/source products. The cheapest new
points have eight source blockers; this is a finite-catalogue statement for this
seed. At cutoff eight, 636 orbits yield 20,160 distinct points.

The full eight-blocker pool has **2,979,840** candidate conflicts and 36,088
shared-blocker clique cuts. Its bounded numerical optimization returns no
improving incumbent. The known 5,730-point seed is retained as a valid fallback.
This pool remains unresolved; a time limit is not an exclusion certificate.
A 512-point control also returns 5,730. Three independently sampled 2,048-point
frontiers are recorded separately. Pool limits, random seeds, solver status,
source hashes and exact output facts are in `evidence/`.

Earlier development attempts with numerical presolve ran past their requested
solver limits and were stopped. `evidence/development_attempts.json` preserves
their status. Presolve is now disabled for this operator, and research commands
also use an outer wall-time limit. Completed larger development runs also gave
no improving incumbent; they do not close their catalogues.

## Rational tail and mixed core/tail images

The preceding image search required integral images. `rational_tail_search.py`
instead proposes exact Householder maps with normals of shapes `(2,1,1)`,
`(2,2,1,1)`, `(3,1,1,1)` and six ones where the dimension permits them.
Their denominators are 3, 5, 6 and 3. Coordinates are scaled exactly, never rounded.
Each trial jointly selects points from the source and its isometric image, using
the existing generic two-copy selector and its checked matching/vertex-cover
certificate. There is **no blocker cutoff**.

Four tail trials in each of 27D, 29D, 30D and 31D, plus four mixed core/tail trials
in 27D, give no gain: respectively 209,825; 231,697; 247,972; and 251,196 points.
Mixed trials transform the entire auxiliary subset and allow source deletions
from every layer. Exact finite optimality applies only to each tested two-copy
union. Different maps, several simultaneous images and other alphabets remain open.

`antipodal_search.py` separately tests the full source/antipode union. The new
antipodal candidates number 143, 147, 480 and 60 in these four dimensions.
Each has a checked injective matching into incompatible source points, so none
of these four particular unions improves its source. Source geometry is inherited
from the pinned existing full certificates; a fresh full replay is required for
any future positive candidate. No worldwide-priority claim is made.

## Reproduce

Use the repository's existing NumPy/SciPy environment. No installation, GPU or
external input is required. Run from the repository root and adapt affinity
to leave at least one logical CPU free.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 timeout 120s \
  python research/partial_orbits_20261104/search17.py \
  --families six,dual,flat,newshort,equator --blockers 8 \
  --candidates 21000 --seconds 20 --out /tmp/new-partial17

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 timeout 180s \
  python research/partial_orbits_20261104/rational_tail_search.py 30 \
  --trials 4 --out /tmp/new-rational30

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 timeout 180s \
  python research/partial_orbits_20261104/rational_tail_search.py 27 \
  --mixed --trials 4 --out /tmp/new-rational-mixed27

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  python research/partial_orbits_20261104/antipodal_search.py 31 \
  --out /tmp/new-antipodal31
```

For a sampled 17D frontier, use `--candidates 2048 --seconds 5 --pool-seed 172701`.
Scripts refuse to replace completed experiment reports. Numerical trajectories
and timing can differ; exact coordinate validity is checked independently.
The manifest pins the committed source and evidence. Historical development
reports retain hashes of their then-loaded engine; the final sampled runs use
the committed implementation. Tests cover actual symmetry breaking, a necessary
two-for-three exchange, direct signed kernel and blocker checks, stabilizers,
the high unsigned phase bit, arbitrary precision, sampling and resource limits.
