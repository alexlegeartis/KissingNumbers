# 248,020 points in dimension 30 from a Barnes–Wall module

13 November 2026. This complete explicit certificate proves
`K(30) >= 248020`, **28 above PackingStar's 247,992** and 48 above the
preceding registered repository value 247,972. The comparison is pinned to
PackingStar commit
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
The retained points are credited to PackingStar. This establishes a strict
improvement over that release; it does not assert global optimality or an
exhaustive worldwide priority survey.

## Complete verification

From the repository root, using Python 3.11+ and NumPy:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 \
  python verifications/improved/dim30-bw14-module/verify.py
```

Adapt the CPU affinity to the execution host and leave at least one logical
CPU free, as required by `AGENTS.md`. The proof uses no optimizer, network
access, external coordinate file, or assumed validity of the source.

`verify.py` checks the complete package manifest, then runs two complementary
checks. `verify_geometry.py` checks every norm, distinctness and all
**30,756,836,190 unordered point pairs** in the final configuration. The
123,394 actual antipodal pairs reduce repeated work; all 1,232 unpaired
points are included explicitly. Squared norms are 128 and distinct-pair
products are at most 64. Bounded integer products and all partial sums are
exactly representable in binary32; there is no numerical tolerance. A full
replay on the discovery host took about 20 seconds on two logical CPUs.

`verify_module.py` independently rebuilds the 1,260-point Barnes–Wall
section, checks its full integer geometry and exact orientation, reconstructs
the pinned source from the final coordinates and delta, and checks all
**238,098,240 inserted-point comparisons** using direct int64 coordinate
sums. This second check is a complete changed-pair verification, not a second
independent enumeration of all unchanged pairs. `data/full_geometry.json`
and `data/module_verification.json` record the fresh successful replays.

`final.py` derives the registered dimension and count from the coordinate
array. `negative_controls.py`, also replayed by `verify.py`, checks rejection of
malformed norms, duplicate points, an invalid same-norm angle and inconsistent
normalization metadata.
`data/manifest.json` pins every other package file. The certificate can be
copied and verified outside this repository.

## Construction

Use squared-norm-32 coordinates for the explanation. Split the first 24
coordinates into the Golay octad
`I = {0,2,4,5,6,10,11,23}` and its complement; append the six tail coordinates.
Inside the resulting 14-dimensional subspace, the incumbent has 1,232 points:
240 equatorial points, 60 pure-tail points and **932 mixed heads**.

Start from the 16-dimensional Barnes–Wall coordinate construction and set
the coordinates labelled 14 and 15 to zero. In the remaining 14 coordinates,
include every vector with two entries ±4 and all other entries zero. Also
include the 128 even-sign vectors with entries ±2 on each affine hyperplane
of `F_2^4` that avoids labels 14 and 15. There are seven such eight-element
supports. The total is

```math
2(14)(13)+7(128)=1260=240+60+960.
```

All norms equal 32. Distinct affine supports intersect in four positions;
within a support, even parity forces at least two sign changes. These facts
bound the relevant products by 16, and the checker also enumerates all
module pairs directly.

Orient the 960 mixed heads by the exact orthogonal matrix `R/4` in
`data/module.npz`, acting on the eight octad and six tail coordinates.
The checker reconstructs the original heads, verifies `R R^T = 16 I`,
and checks the transformed coordinates by exact integer division. This
orientation has **zero external blockers**: every incumbent point outside
the old 932-head set is compatible with every inserted point. Therefore

```math
247992-932+960=248020.
```

This is a coordinated module replacement: all 932 removed points and all
960 inserted points are distinct, and there is no overlap between the old
and new head sets. The complete proof does not rely on the search's
restriction to 152 potentially obstructing exterior points.

## Files and discovery

`data/coordinates30_248020.npz` contains the complete final certificate in
twice the explanatory coordinates, with squared norm 128 and product limit
64. `data/exchange.npz` stores the 932 deletions and 960 insertions in the
original norm-32 coordinates. `data/module.npz` contains the oriented heads,
their original Barnes–Wall heads and the exact transformation. Source hashes,
counts and attribution are in `data/provenance.json`.

The bounded search sampled E8 Weyl transformations and signed tail
permutations, finding the zero-exterior-cost orientation at iteration 11,025.
`data/search_report.json` preserves the original run. Reproducible discovery
code and subsequent experiments are in
[`../../../research/reclaim_26_31_20261113/packingstar_audit/`](../../../research/reclaim_26_31_20261113/packingstar_audit/).
Zero exterior cost is optimal for this fixed full-module objective; it is
not an upper bound on more general configurations or partial-module exchanges.
