# Exact improvements over PackingStar in all dimensions 27–31

29 October 2026. Five explicit configurations exceed PackingStar’s **28 October release**:

| Dimension | PackingStar | This construction | Gain | Removed | Added |
|---|---:|---:|---:|---:|---:|
| 27 | 209,760 | **209,807** | **47** | 244 | 291 |
| 28 | 221,864 | **221,988** | **124** | 312 | 436 |
| 29 | 231,676 | **231,697** | **21** | 219 | 240 |
| 30 | 247,967 | **247,972** | **5** | 512 | 517 |
| 31 | 251,169 | **251,188** | **19** | 463 | 482 |

The sextet continuation reclaims 28D by **124 points** and further improves
27D, 29D and 31D. These are explicit lower bounds; no unrestricted optimum
is claimed.

PackingStar supplies the retained points. Its public inputs are pinned to
commit [`477df34`](https://github.com/CDM1619/PackingStar/tree/477df34cf2a09bb9476b735a685b538a82da8c55/25D-31D/new_bounds_20261028).
Our contributions are the explicit replacements above. `data/provenance.json`
records input SHA-256 hashes, canonical coordinate hashes, deletions,
insertions and the preceding repository values. Each complete configuration
is included; verification does not assume the public input is valid.

## What changed

Write I = {0,2,4,5,6,10,11,23}. In the initial 27D and 29D checkpoint,
new points occupied I and the tail coordinates. The sextet continuation
also changes points outside this subspace. The search enumerates even/odd E8-coordinate shells and
additional even-coordinate shells, with integer tails of complementary
squared norm. It allows the deletion of conflicting points outside this
subspace as well as points inside it. A conflict graph records every signed
angle violation, and individual points are selected without imposing
antipodal symmetry. Shared deletions and larger perturbations produce the
explicit improvements; the search itself is not used as a proof.

For the initial 30D and 31D checkpoint, complete reflected copies of
PackingStar's configuration provide additional candidates. The three commuting generators negate I,
negate all tail coordinates, and negate the whole vector. Their eight
images are valid individually. The union need not be valid: the selected
subset excludes all cross-image conflicts. This finite selection gives
the stated gains. The exact coordinate replay independently checks the
result without relying on the selection algorithm or its optimization status.

### Corrected sextet symmetries

A Golay sextet partitions the 24 Leech coordinates into six tetrads. On
each tetrad apply S = J/2 - I_4, where J is the four-by-four all-ones matrix,
and then negate the first tetrad. The sign correction is essential for
preserving the standard Leech shell. S is exactly orthogonal because
J^2 = 4J. The tail coordinates are fixed.

For a valid configuration X, its image T(X) is valid too. Separate their
common points and make the bipartite graph of forbidden products between
the two remaining parts. Maximum matching and its alternating paths give
a largest compatible subset of this two-copy union by Konig's theorem.
Iterating different sextets allows coordinated changes to Leech points
and midpoint layers. The final certificates below are verified without
assuming this discovery argument or the input configuration is valid.

`search_sextets.py` implements the search with exact bounded-integer products
and SciPy bipartite matching. It verifies both the normalized input and final
output completely. For example, the first frame reproduces 221,896 points
from the pinned 28D input, with exactly the same canonical coordinate hash
as our discovery run:

```bash
python search_sextets.py D28_221864_integer_coordinates.npz --frames 1 --output trial28.npz
```

All 35 frames are used by default. `--rounds` enables further sweeps;
`--max-blockers` optionally restricts candidates by their conflicts with
the current configuration. A cutoff gives a restricted search, not an
unrestricted ceiling. Fractional images are skipped rather than rounded.
The exact final configurations are reproduced from the deltas below.

## Complete verification

```bash
python final.py
python verify.py
python negative_controls.py
```

Python 3.9+ and NumPy suffice. Run without `-O`. Coordinates are integer
numerators with squared norm **128**. Every distinct-pair product is at
most **64**, so normalization by sqrt(128) gives a kissing configuration.
The denominator-two field records the common scale relative to the public
files; it does not enter the angle proof.

`verify_geometry.py` checks every norm and distinctness condition, identifies
the opposite pairs that actually occur, and checks all remaining signed
or absolute products. It does not assume antipodality. With P opposite
pairs and U unpaired points, its comparisons cover exactly

    4*binom(P,2) + P + 2*P*U + binom(U,2) = binom(2*P+U,2).

Every coordinate is bounded by 16 and the dimension by 31. Each product
and partial dot-product sum is therefore an integer of magnitude at most
31*16^2 = 7,936, below 2^24. Binary32 matrix multiplication represents all
these values exactly. There is no geometric tolerance. `--precision 64`
replays the complete comparison with binary64; both formats were checked
on the final certificates. They share the same enumeration logic.

`verify_delta.py` separately reconstructs the pinned baseline coordinate
hash and checks every inserted point against every final point using direct
int64 coordinate sums, without BLAS. This is an independent cross-check
of all changed inequalities; the complete geometry replay establishes
the unchanged inequalities too. The five configurations cover
**135,782,755,339 unordered point pairs** in total.

The default verifier checks the complete file manifest, regenerates each
geometric report and compares every deterministic field, and reruns the
direct-integer checks. Reports are written only when explicitly requested
from `verify_geometry.py`; ordinary verification leaves the package unchanged.
Negative controls require both arithmetic formats to reject a duplicate,
a wrong norm, and two distinct norm-128 points with product 96 > 64.

## Rebuild

Each `data/exchange<D>.npz` contains exact removed and added numerator
coordinates. With the pinned public input, for example:

```bash
python rebuild.py 27 D27_209760_integer_coordinates.npz --output rebuilt27.npz
python verify_geometry.py rebuilt27.npz
```

The rebuilder checks the downloaded file's SHA-256, removed membership,
insertion distinctness and the final canonical coordinate hash. Rebuilding
and verification require no optimizer, network retrieval or search state.
