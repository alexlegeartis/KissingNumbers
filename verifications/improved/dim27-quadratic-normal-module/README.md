# A verified 209,600-point kissing configuration in 27 dimensions

**Retained construction.** The [October 29 continuation](../dim27-31-joint-refills/) updates 27–31D against PackingStar’s later release; the earlier comparisons and certificates below are historical. The 28D sextet continuation now also exceeds the later PackingStar release.

28 October 2026. The explicit configuration has **209,600 points**, improving
PackingStar's latest checked **209,586** by **14**. It uses coordinates in
Q(sqrt(2)). This is a verified lower bound, not an optimality claim.

The comparison is pinned to PackingStar commit
[`b67832ad9970895a22ab667bb06d1d43a45c6b18`](https://github.com/CDM1619/PackingStar/tree/b67832ad9970895a22ab667bb06d1d43a45c6b18/25D-31D/new_bounds_20261026).
PackingStar supplies the retained core of 209,216 points. Our contribution
is a new 384-point normal module replacing their 370-point normal sector.
The preceding repository construction has 209,532 points.

## Complete independent verification

```bash
python final.py
python verify.py
```

Python 3.9+ and NumPy suffice. Run without `-O`. The verifier checks the
explicit coordinates without assuming validity of the baseline, the search,
or a geometric reduction. It checks all norms, distinctness, antipodal
partners, and all **5,491,467,600** distinct pairs of the **104,800**
antipodal lines. The complete replay took approximately 108 seconds on the
discovery machine; its report is in `data/verification27_209600.json`.

The NPZ stores integer arrays A and B representing X = A + B sqrt(2).
Every squared norm is exactly 32, and every distinct-line absolute inner
product is at most 16. After division by sqrt(32), these are unit vectors
with pairwise inner products at most 1/2, giving a kissing configuration.
There are 144 irrational points and 209,456 integer points.

Small integer dot products are evaluated exactly by binary64 BLAS; the
verifier bounds every coefficient so all products and partial sums are
exactly representable. For a quadratic comparison a + b sqrt(2) <= L,
it first checks the signs of b and L-a, then squares only the appropriate
nonnegative quantities. No floating-point tolerance or approximate square
root is used. All integer squaring is within the checked int64 bounds.

## Construction

Let I = {0,2,4,5,6,10,11,23}. The normal sector consists of vectors supported
on I and the three tail coordinates. All points outside this subspace are
retained from PackingStar's pinned file.

The new module consists of:

* The 240 E8 roots of squared norm 32 on I: two entries +/-4, or eight
  entries +/-2 with even minus parity.
* 144 caps with head squared norm 16 and tail squared norm 16. Each cap
  belongs to a block (+/-h, +/-t), with independent signs. The 36 blocks
  use E8 norm-16 heads and tails t = 2 sqrt(2)(e_i +/- e_j).

This changes the tail directions and radii from the incumbent's integer
module. For each of the six tail lines, all 1,080 E8 norm-16 head lines
were checked against the fixed core using exact quadratic comparisons.
The resulting 76 admissible head/tail blocks were optimized jointly;
36 fit. Optimality of this small block model is not a ceiling on other
modules or on the kissing number.

An integer normal-sector search also produced an independently verified
209,598-point configuration. Combining its sector with the quadratic
module by exact bipartite matching gave the same 209,600 count. The
packaged certificate is independently checked regardless of this search
history.

## Rebuild and provenance

`data/coordinates27_209600.npz` contains the full certificate.
`data/normal_module.npz` contains the explicit 384-point replacement.
`data/provenance.json` pins the public input and coordinate hashes.
With PackingStar's pinned file from the linked commit:

```bash
python rebuild.py D27_209586_integer_coordinates.npz --output rebuilt27.npz
python verify.py rebuilt27.npz
```

The rebuild checks the public file's SHA-256 and the canonical coefficient
hash of the resulting configuration. No optimizer is needed to rebuild
or verify the construction. The separate [asymmetric exchanges](../dim26-30-asymmetric-exchanges/) now also exceed the
pinned PackingStar files in 26D and 28D-30D.
