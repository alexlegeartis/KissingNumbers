# Dimensions 26 and 27: finite exchange exclusions, 14 November 2026

**This campaign did not outperform PackingStar. No new lower bound is claimed.**
The checked targets remain 201,180 points in dimension 26 and 209,862 in
dimension 27. The repository's registered bounds are unchanged.

The positive result is an exact closure of three explicit candidate pools.
The 2,004-point 26D octad-transplant pool left unresolved on November 13
cannot improve the source, even when every source point may be deleted.
The same holds for larger merged pools with 4,662 points in 26D and 5,230
points in 27D. This excludes these finite unions, not other coordinates,
larger pools, deformations, or the global kissing-number problem.

| Pool | New candidate points | Exact rational gain bound | Integer gain bound |
|---|---:|---:|---:|
| Original November 13 octad pool, 26D | 2,004 | 33 / 8,388,608 | **0** |
| Merged pool, 26D | 4,662 | 47 / 4,194,304 | **0** |
| Merged pool, 27D | 5,230 | 347 / 16,777,216 | **0** |

Here gain means the number of admitted new points minus the number of
deleted source points. These are exact rational bounds, rounded down only
because the gain is an integer. They are independently checked certificates,
not numerical MILP optimality messages.

## Replay

Only Python and NumPy are needed for the checker. From the repository root,
adjust the affinity to the host and leave at least one logical CPU free:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python research/26_27_exchange_exclusions_20261114/verify.py

taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python research/26_27_exchange_exclusions_20261114/test_verify.py
```

The supplied `verification.json` is a receipt from a fresh replay. The
checker hashes the actual inputs and recalculates every inequality used in
the proof. Corruption controls reject wrong norms, negative weights, false
blockers, compatible members of a claimed conflict clique, missing coverage,
and an incorrect declared bound.

To construct another dual certificate from a saved pool, SciPy is also needed:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python research/26_27_exchange_exclusions_20261114/build_dual.py \
  research/26_27_exchange_exclusions_20261114/data/joint26.npz \
  --out /tmp/joint26-rebuilt-dual.json
```

The LP may choose different dual weights. Its proposed weights are rounded
to multiples of 2^-24; exact nonnegative repairs using indicator bounds
restore coefficient domination. The checker uses `Fraction`, without an
optimization library or a geometric tolerance.

## Why the proof works

Scale the coordinates to squared norm 128, with product limit 64. Let `x_i`
indicate inclusion of a new point and `y_j` indicate deletion of source
point `j`. If a group of new points is a conflict clique and every member
conflicts with source point `j`, then

```text
sum(x_i over the group) <= y_j.
```

At most one member can be admitted, and admitting it forces deletion of the
owner. A clique without an owner instead has sum at most one. Each stored
proof is a nonnegative rational combination of these inequalities and the
individual bounds `x_i <= 1`, `y_j <= 1`. Its coefficients dominate
`sum(x_i) - sum(y_j)` and its right side is strictly below one.

`verify.py` checks all candidate norms, distinctness, absence from the
source, every owner product, every within-clique product, all rational
weights, and coefficient domination. For the sixteen quadratic source
points in 26D, signs in Q(sqrt(3)) are decided with integer squares.
Unmentioned source deletions only decrease the gain. The proof therefore
allows deletion of **any** source point. It does not require complete
blocker lists or the geometric validity of any search oracle.

This checker does not perform a complete all-pairs verification of either
PackingStar source. That is unnecessary for an upper bound on the cardinality
of an admissible subset of these literal finite unions. It must not be
reported as a fresh proof of their lower bounds.

## Candidate pools and completed searches

The merged pools combine individually selectable Golay flat midpoints,
signed six-core/two-tail points, and integer E8-body candidates in the normal
octad sector. The 26D pool also includes the November 13 octad transplants.
All admitted new points have at most four source blockers. Duplicate points
are removed, and deletion costs are shared across families.

The archived search receipts record the following additional experiments:

| Experiment | Completed scope | Outcome |
|---|---|---|
| Golay flat-midpoint scan | All 1,130,496 two-error bodies, with 4 tail directions in 26D and 12 in 27D; blocker cutoff 4 | 1,792 and 3,584 new candidates; no improving numerical selection |
| 27D sextet images | All 35 frames compatible with the normal octad; cutoff 32 | No gain |
| Signed tail permutations | All 8 in 26D and all 48 in 27D; cutoff 8 | No gain |
| Integer normal-sector catalogue | 21,960 and 57,046 new exterior-compatible points; only candidates costing at most four deletions enter the exact merged-pool proof | No gain in the bounded searches; no optimality claim for the full catalogues |
| Half-integral mixed-parent midpoints | All 169,856 mixed-parity contact pairs with both parents absent in 26D, and all 1,394,432 in 27D; 4 and 12 tail directions | No candidates at cutoff 4; none at cutoff 8 in 26D |
| Release the 26D norm-16 head layer | Same 169,856 half-integral bodies; all norm-16 heads and sixteen radicals may be deleted, while the other layers stay fixed | No admitted candidates |
| Continuous 26D tail-circle pilot | All 169,856 bodies; numerical polygon-vertex screening | No free insertion found; this numerical receipt is not an exact exclusion |

The exact result is restricted to the explicit coordinate arrays in `data/`.
The other receipts document completed finite searches and do not certify
unrestricted optima or complete candidate-generation coverage. The pool
files preserve blocker lists and conflict edges for continuation, but the
independent checker does not trust those lists.

## Source identities and attribution

The source is PackingStar commit
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111),
whose release is labelled November 11, 2026:

* `D26_201180_coordinates.npz`: SHA-256
  `36732203081894af188e75e2dae69327961b79c7b0781385b1589025b7e0b309`.
* `D27_209862_coordinates.npz`: SHA-256
  `79f8c596a3bdffac1005032af0edaa8790d3c724f146115898237905982e757e`.

They are bundled unchanged as `data/source26.npz` and `data/source27.npz`.
The constructions are due to Ma and collaborators; see their
[PackingStar paper](https://arxiv.org/abs/2511.13391). They are inputs, not
new constructions from this campaign. `manifest.json` pins all proof inputs.

The original transplant pool comes from the repository's
[`continuous26/transplant759.npz`](../reclaim_26_31_20261113/continuous26/transplant759.npz).
This continuation began at repository commit
`921b1b653f15b599b78953139448be71d8f3812e` and used four of five available
logical CPUs at most, reserving CPU 4. No search jobs remain running.
