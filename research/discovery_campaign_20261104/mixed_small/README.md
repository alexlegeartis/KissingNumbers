# Mixed-orbit 11D/12D numerical search wave, 2026-11-04

**Outcome: no new SOTA and no admissible target-size result.** The best overall
maximum cosine remains the published input, 0.5009010186022461.
All 24 saved point sets were independently replayed and all fail cosine <=1/2.

This directory contains actual numerical search runs targeting at least 605
points in 11D (comparison 604) and 842 in 12D (comparison 841). A target point
count is not an admissible kissing configuration. Numerical arrays and objective
values are not exact certificates.

The prior `verifications/closed/dim09-19-record-maximality/README.md` and
KNOWLEDGE sections 63–64 were read first. Those notes close simple append-to-841
and fixed-record-equator replacements, and price the two-block 840 family.
The first phase imports no record coordinates and uses none of those templates.
The 841-point file was not read as an optimization input.

## Fresh families actually searched

Each real cyclic representation is a direct sum of independent rotation planes
and a trivial subspace. Frequencies may repeat, so multiplicity spaces and every
seed's amplitude/phase are independently variable. Composite cyclic groups permit
mixed orbit lengths by putting each seed in the fixed space of the corresponding
power of the generator. Periods divide the group order; a rotation plane of
frequency f supports period L exactly when f*L is divisible by the group order.

| Dimension | Group | Orbit forecast | Target points |
|---|---|---|---:|
| 12 | C15 | 55 full 15-orbits, 3 short 5-orbits, 2 fixed points | 842 |
| 11 | C21 | 28 full 21-orbits, 2 short 7-orbits, 1 short 3-orbit | 605 |
| 12 | C7 | 120 full 7-orbits, 2 fixed points | 842 |
| 11 | C5 | 121 full orbits, repeated irreducibles | 605 |
| 12 | C11 | 76 full 11-orbits, 6 variable fixed-plane points | 842 |
| 12 | C13 | 65 full orbits | 845 |
| 11 | C11 | 55 full orbits | 605 |
| 12 | two independently rotated C5 actions | 168 full 5-orbits, 2 fixed points | 842 |

These forecasts only describe the number of optimization variables/points.
Closed group action makes seed-versus-orbit checks sufficient during the common
frame phase, with correct orbit-length weights. The independently rotated case
uses all pairs even during optimization. Every completed trial is separately
screened against every pair. Quadratic positive-violation and quartic objectives
alternate; analytic gradients were checked by finite differences.

Three selected outputs then have their group symmetry explicitly broken by
small random point perturbations. Every individual vector becomes free, with
joint normalization and full-pair objectives. These runs preserve no antipodality
or group equations. A further augmented-Lagrangian/smooth-maximum phase was
interrupted after about 54 seconds when the parent redirected the budget to a
competitive published near-miss. Its partial log is retained; no unexported
intermediate point set is claimed as a candidate.

## Competitive published 842-point near-miss branch

The later branch downloads `beyond841/Gram842.txt` from the primary repository
https://github.com/k-nic/841_in_12D . Source URL, byte count, SHA256 and Git blob
identity are recorded in `near_miss_source.json`, and the original source bytes
are preserved as `Gram842_source.txt`. The Git blob hash was independently
recomputed and matches the API metadata.

The 12-decimal source matrix has actual maximum off-diagonal entry
**0.500901018602**; the historical shorthand 0.500900 is not used as a precise
baseline. Its omitted eigenvalues reach approximately ±1.63e-11 from decimal
rounding. Taking the largest twelve eigenvectors and renormalizing rows gives
842 vectors in R^12 with maximum cosine **0.5009010186022461** and Gram error
6.54e-13. These reconstruction errors are explicitly numerical. No exact rank
or positive-semidefinite certificate is inferred from rounded source decimals.

All 842 point vectors are optimized jointly. The schedules include direct
pressure, a coupled anisotropic deformation followed by constraint homotopy,
low-degree cluster reseeding, high-pressure contact cluster reseeding, coupled
opposite clusters, and global symmetry breaking. A relaxed threshold is reduced
toward 1/2, followed by a smooth-maximum objective. This is not an append test on
the 841-point record and does not preserve any immutable base subset.

The six near-miss modes took 60.18 seconds. Four final modes added high-beta
minimax polishing and coupled external-force drives through relaxed constraints,
then removed those forces and compressed every point; these took 43.06 seconds.
Neither phase improved the source maximum cosine. The best retained input has
squared positive-violation sum 0.007805746737110224 both before and after.
The fresh group phase ran 44 trials in 84.08 seconds; the free releases ran
11 trials in 167.70 seconds. Fresh best maxima were 0.5608616583 for 605 points
in 11D and 0.5589977855 for 842 points in 12D.

`near_miss_results.json` and `driven_results.json` record the actual initial/final maximum, violation
objectives, phase limits and per-trial results. The global best maximum is tracked
at every objective evaluation; a worse final objective minimizer cannot silently
replace a better maximum-cosine point set. All saved best point sets receive a
separate full-pair replay in `replay_results.json`.

## Reproduction and scope

Uses existing Python, NumPy, SciPy; no installation, push, or registry edits.
Commands run from this directory:

    OPENBLAS_NUM_THREADS=1 python search.py --phase group --seconds 35
    OPENBLAS_NUM_THREADS=1 python search.py --phase release --seconds 65 --release-count 3
    OPENBLAS_NUM_THREADS=1 python near_miss_search.py
    OPENBLAS_NUM_THREADS=1 python driven_homotopy.py
    OPENBLAS_NUM_THREADS=1 python replay.py

The source Gram is fetched separately from the pinned primary URL in its metadata.
The scripts retain seeds, configurations, frames, orbit periods, phase limits,
wall times and full-pair numerical screens. Source/reconstruction hashes support
replay. No negative numerical run is a nonexistence theorem. No new SOTA is
claimed unless a feasible result is independently certified.
