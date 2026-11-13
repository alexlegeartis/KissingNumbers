# 248,084 points in dimension 30 from a Leech–Sloane section

13 November 2026. This complete explicit certificate proves **K(30) ≥ 248,084**,
**92 above PackingStar's 247,992**, 112 above the preceding published
repository certificate, and 64 above the earlier 248,020-point Barnes–Wall
construction from this campaign. The comparison is pinned to
[PackingStar commit 2f3453b](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
The retained source points are credited to PackingStar. This comparison does
not establish an exhaustive worldwide priority survey or global optimality.

## Verification

From the repository root, with Python 3.11+ and NumPy:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 \
  python verifications/improved/dim30-leech-sloane-section/verify.py
```

Adapt CPU affinity to the host and reserve one logical CPU under `AGENTS.md`.
The package needs no optimizer, network access, external coordinate file, or
assumption that the PackingStar source is valid.

`verify_geometry.py` checks the complete 248,084-point configuration, covering
**30,772,711,486 unordered pairs**. Squared norms equal 128 and distinct-pair
products are at most 64. There are 123,434 actual antipodal pairs and 1,216
unpaired points. All bounded integer products and intermediate sums are
exactly representable in binary32; the check uses no numerical tolerance.

`verify_module.py` independently rebuilds the 2,564-point parent from its
Steiner quadruples and affine binary glue codes, checks every parent pair,
takes the exact 1,356-point coordinate section, checks every section pair,
verifies its rational orientation, reconstructs the pinned source from the
minimal exchange, and checks every inserted point against all other final
points using direct int64 sums. The latter covers **261,975,648 comparisons**.
This independently verifies every changed pair; the complete geometry checker
establishes every unchanged pair as well.

The combined `verify.py` checks every package hash, both fresh proof receipts,
and four malformed-input controls. `final.py` reads the count from the actual
coordinate array. The package is self-contained and portable.

## Construction

The parent is the 2,564-point Leech–Sloane recipe reassembled to expose an
eight-dimensional E8 equator. Its exact reconstruction is described in the
[31D companion package](../dim31-leech-sloane-module/), and the complete recipe
and constructor are also included here. No new 15D kissing bound is claimed.

Set its seventh tail coordinate to zero and drop that coordinate. The
fourteen-dimensional section contains

**240 equatorial points + 1,056 heads + 60 pure-tail points = 1,356.**

Apply the exact orientation R/4 stored in `data/module.npz`, with RRᵀ = 16I.
Embed its first eight coordinates on the Golay octad
`{0,2,4,5,6,10,11,23}` and its last six on ambient coordinates 24 through 29.

PackingStar's corresponding sector contains 1,232 points: the same 240-point
equator, the same 60-point tail set, and 932 heads. The new oriented section
has 32 external blockers. Cancelling the 300 unchanged equator and tail points
leaves the exact exchange

**247,992 − (932 + 32) + 1,056 = 248,084.**

The complete verifier checks compatibility against every retained point; it
does not rely on the reduced exterior used during discovery.

## Why this exceeds the full-fiber model

The preceding Barnes–Wall module has 960 heads and reaches a restricted bound
for complete 16-point fibers over the fixed integer tail alphabet. The new
section uses partial fibers:

| Points in a fiber | Number of tails | Head contribution |
|---:|---:|---:|
| 2 | 192 | 384 |
| 10 | 48 | 480 |
| 16 | 12 | 192 |
| **Total** | **252** | **1,056** |

This is outside the hypotheses of that complete-fiber bound. The 96 additional
heads incur 32 new exterior deletions, improving the Barnes–Wall result by 64.

## Files and discovery

`data/coordinates30_248084.npz` contains the entire norm-128 configuration.
`data/exchange.npz` stores 964 deletions and 1,056 insertions at that scale.
`data/module.npz` stores the norm-32 section and its exact transformation.
`data/construction.json` stores the parent recipe; `data/provenance.json`
identifies the omitted tail, normal octad, source and exchange counts.
The exact proof receipts and complete manifest are alongside them.

Seven coordinate sections and their initial orientations were priced before
a bounded 120-second orientation search found the 32-blocker configuration.
A separate matching certificate shows that the literal union with the earlier
Barnes–Wall witness cannot improve this count. That exclusion applies only to
those two fixed configurations. Code, receipts and the matching certificate
are in the [research directory](../../../research/reclaim_26_31_20261113/continuous26/).

The malformed-input checks are implemented in `negative_controls.py`.
