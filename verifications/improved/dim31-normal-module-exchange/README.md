# A verified 251,110-point kissing configuration in 31 dimensions

**Retained construction.** The [October 29 continuation](../dim27-31-joint-refills/) updates 27D and 29–31D against PackingStar’s later release; the earlier comparisons and certificates below are historical. The 28D repository result remains surpassed.

28 October 2026. The explicit integer configuration has **251,110 points**,
exceeding PackingStar's latest checked **251,002** by **108** and the preceding
repository construction, 246,638, by 4,472. This is a lower bound, not an
optimality claim.

The comparison is pinned to PackingStar commit
[`b67832ad9970895a22ab667bb06d1d43a45c6b18`](https://github.com/CDM1619/PackingStar/tree/b67832ad9970895a22ab667bb06d1d43a45c6b18/25D-31D/new_bounds_20261026),
whose 31D coordinate file was checked on 28 October. The construction extends
PackingStar's coordinates; their retained points are credited as the baseline.
Our contribution is a coordinated replacement by a larger explicit module and
exact selection of compatible points from the combined configurations.

## Complete verification

```bash
python final.py
python verify.py
```

Python 3.9+ and NumPy suffice. Run without `-O`. The geometric verifier uses no search software,
external baseline, network access, or inherited validity assumption. It checks
all norms, distinctness and antipodal partners, then all **7,881,966,235**
distinct pairs of the 125,555 antipodal lines. All squared norms equal 32;
all distinct-line absolute inner products are at most 16. Scaling by
1/sqrt(32) gives the required unit vectors with inner products at most 1/2.

The arithmetic has no tolerance: coordinates are small integers, and every
product and partial sum is exactly representable in binary64. The supplied
verifier explicitly checks the bounds that justify this. The recorded complete
replay took approximately 111 seconds on the discovery machine. A report is
in `data/verification31_251110.json`; `verify.py` recomputes it independently.

## Construction

Split the 24 Leech coordinates into the octad
`I = {0,2,4,5,6,10,11,23}` and its 16-coordinate complement. The new module
uses only I and the seven tail coordinates, a 15-dimensional subspace.

Label the eight octad positions by the binary vectors 8 through 15 and the
seven tail positions by 1 through 7. Each nonzero binary linear functional
on these labels has an eight-position support. Exclude the functional whose
support is the octad alone. Each of the remaining 14 supports contains four
octad and four tail positions. Put all 128 even-sign vectors with entries
plus or minus 2 on each support. This gives **1,792 points**. Add all
**224** vectors with one octad entry plus or minus 4 and one tail entry plus
or minus 4. The resulting module contains **2,016 points**.

Every module vector has squared norm 32. Distinct supports intersect in four
positions, giving product at most 16. Within one support, even sign parity
forces at least two changed signs, again giving product at most 16. The
224 additional vectors and their products with the signed supports obey the
same bound directly.

Most of the incumbent clears this entire module. Its regular midpoint heads
restrict to scaled E8 roots on I, and its residual heads are supported on the
complement. This makes it possible to replace a whole layer profitably,
although individual-point insertion searches found no additions.

For each union of two valid antipodal configurations, cross-conflicts form a
bipartite graph. Maximum matching and an alternating-path minimum vertex
cover select a maximum compatible subset **of that particular union**.
The first module exchange gave 251,106 points; union with PackingStar's
28 October update gave four more. The complete verifier above independently
checks the final coordinates, so it does not rely on this reduction.

## Rebuild and provenance

`data/coordinates31_251110.npz` is the complete coordinate certificate.
`data/exchange_delta.npz` records the exact points removed from and added to
PackingStar's pinned 251,002-point file. Obtain that file from the linked
commit and run:

```bash
python rebuild.py D31_251002_integer_coordinates.npz --output rebuilt31.npz
python verify.py rebuilt31.npz
```

The rebuild checks the public input's SHA-256 and the canonical hash of the
resulting point set. `data/provenance.json` records the source commit, blob
hash, certificate hashes, counts and current comparison values in 26D–31D.
This package addresses 31D. The separate [27D quadratic module](../dim27-quadratic-normal-module/)
gives 209,600; the [asymmetric exchanges](../dim26-30-asymmetric-exchanges/) also exceed
the pinned PackingStar files in 26D and 28D–30D.
