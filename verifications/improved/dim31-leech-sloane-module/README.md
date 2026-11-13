# 251,395 points in dimension 31 from a Leech–Sloane module

13 November 2026. This complete explicit certificate proves **K(31) ≥ 251,395**,
**141 above PackingStar's 251,254** and 195 above the preceding published
repository certificate, 251,200. The comparison is pinned to
[PackingStar commit 2f3453b](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
PackingStar supplies the retained source configuration. The local construction
uses the Leech–Sloane 2,564-point recipe, reassembled to expose the required
E8 equator. This is a new 31D configuration; no new 15D bound is claimed.
The comparison does not assert an exhaustive worldwide priority survey or
global optimality.

## Verification

From the repository root, with Python 3.11+ and NumPy:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 taskset -c 0-1 \
  python verifications/improved/dim31-leech-sloane-module/verify.py
```

Adapt affinity to the host and leave at least one logical CPU free, as required
by `AGENTS.md`. No optimizer, network access, external source coordinate file,
or assumed validity of the PackingStar source is needed for this proof.

`verify_geometry.py` checks all norms, distinctness, and every one of the
**31,599,597,315 unordered point pairs**. The stored integer coordinates have
squared norm 128 and distinct-pair products at most 64. There are 125,667
actual antipodal pairs and 61 unpaired points; the latter are checked explicitly.
The bounded products and every intermediate sum are integers exactly
representable in binary32. Acceptance uses no floating-point tolerance.

`verify_module.py` is a separate direct-integer reconstruction and changed-pair
proof. It builds the complete 2,564-point module from the saved combinatorial
recipe, checks all **3,285,766 module pairs**, verifies the exact rational
orientation, reconstructs and hashes the pinned source after undoing the
exchange, and checks **557,843,286 added-point comparisons** using int64 sums.
This is an independent check of all changed pairs, rather than a second
enumeration of every unchanged source pair. The full geometry proof above
establishes those unchanged pairs without an inherited assumption.

The combined `verify.py` also checks the complete manifest and four malformed
inputs: a wrong norm, duplicate points, a same-norm forbidden angle and
inconsistent normalization metadata. It compares fresh receipts with the
saved receipts, ignoring only runtime. The complete package is portable.

## Reconstructing the local module

Work initially in 15 coordinates with quadratic form diag(1,…,1,2).
Seven pairs partition the first fourteen coordinates. The saved exact cover
contains 70 transversal quadruples, each using four different pairs. Together
with the 21 unions of two pairs, they form 91 quadruples covering every
three-element subset of the fourteen coordinates exactly once.

The weighted construction contains:

| Class | Count |
|---|---:|
| Signed coordinate axes ±4eᵢ, for i among the first fourteen coordinates | 28 |
| All sixteen sign patterns ±2 on each of the 91 quadruples | 1,456 |
| Two affine binary glue codes, each with 512 sign words | 1,024 |
| Three-coordinate points on one coordinate pair and the last coordinate | 56 |
| **Total** | **2,564** |

The recipe stores each affine code as one origin and nine independent binary
generators. The verifier expands the codes directly; it does not trust a
solver's feasibility or population report.

Take pair sums as seven Euclidean coordinates, twice the final weighted
coordinate as the eighth, and pair differences as the last seven. This map
doubles the weighted norm and gives squared Euclidean norm 32. The resulting
module has **240 E8 equator points, 2,240 heads, and 84 D7 tail points**.
Both the equator and pure-tail set are checked against explicit root formulas.

## Embedding and exchange

The eight equatorial coordinates occupy the normal Golay octad
`{0,2,4,5,6,10,11,23}` in the first 24 ambient coordinates. The seven tail
coordinates occupy positions 24 through 30, using zero-based indexing.
The stored orientation is the exact isometry R/2, with RRᵀ = 4I.

The final construction retains **2,559 of the 2,564 module points**. Keeping
individual points permits a better exchange than requiring the whole module
or antipodal pairs. There are 340 coincidences with the source, 2,219 actual
insertions, and 2,078 actual deletions. Of the deletions, 1,942 belong to the
old module and 136 lie outside it. Thus

**251,254 − 2,078 + 2,219 = 251,395.**

Equivalently, replace the old 2,282-point sector by the admitted 2,559-point
module subset and delete 136 external blockers. Five otherwise new module
heads are omitted. The independent checker enforces the subset interpretation
and does not silently assume the whole module survives.

## Files and provenance

`data/coordinates31_251395.npz` is the complete norm-128 configuration.
`data/exchange.npz` contains the minimal delta at the same normalization.
`data/module.npz` contains the full norm-32 oriented module and R/2.
`data/construction.json` is the self-contained combinatorial recipe, including
the original 15D recipe attribution and source hash. `data/provenance.json`
pins the PackingStar source, coordinate hash, counts and subset convention.
The two verification receipts and malformed-input receipt record successful
replays. `final.py` reads the registered count from the coordinate array.
`data/manifest.json` hashes every other package file.

The first full-module transplant gave 251,358 points. A bounded 90-second
orientation and exact bipartite-union search tested 3,530 orientations and
produced this 251,395-point witness. Search scope and reproducible code are in
the [research campaign](../../../research/reclaim_26_31_20261113/algebraic/).
No search optimum beyond the explicitly evaluated finite unions is claimed.

The malformed-input checks are implemented in `negative_controls.py`.
