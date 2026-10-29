# Further exact improvements over PackingStar in 27D and 29–31D

29 October 2026. Four explicit configurations exceed PackingStar's
**28 October release**:

| Dimension | PackingStar | This construction | Gain | Removed | Added |
|---|---:|---:|---:|---:|---:|
| 27 | 209,760 | **209,775** | **15** | 104 | 119 |
| 29 | 231,676 | **231,684** | **8** | 97 | 105 |
| 30 | 247,967 | **247,972** | **5** | 512 | 517 |
| 31 | 251,169 | **251,180** | **11** | 159 | 170 |

These are modest further lower-bound improvements. The broader search has
not beaten PackingStar's **221,864** in 28D. No unrestricted optimum or
larger structural jump is claimed.

PackingStar supplies the retained points. Its public inputs are pinned to
commit [`477df34`](https://github.com/CDM1619/PackingStar/tree/477df34cf2a09bb9476b735a685b538a82da8c55/25D-31D/new_bounds_20261028).
Our contributions are the explicit replacements above. `data/provenance.json`
records input SHA-256 hashes, canonical coordinate hashes, deletions,
insertions and the preceding repository values. Each complete configuration
is included; verification does not assume the public input is valid.

## What changed

Write I = {0,2,4,5,6,10,11,23}. In 27D and 29D, new points occupy I and the
tail coordinates. The search enumerates even/odd E8-coordinate shells and
additional even-coordinate shells, with integer tails of complementary
squared norm. It allows the deletion of conflicting points outside this
subspace as well as points inside it. A conflict graph records every signed
angle violation, and individual points are selected without imposing
antipodal symmetry. Shared deletions and larger perturbations produce the
explicit improvements; the search itself is not used as a proof.

In 30D and 31D, complete reflected copies of PackingStar's configuration
provide additional candidates. The three commuting generators negate I,
negate all tail coordinates, and negate the whole vector. Their eight
images are valid individually. The union need not be valid: the selected
subset excludes all cross-image conflicts. This finite selection gives
the stated gains. The exact coordinate replay independently checks the
result without relying on the selection algorithm or its optimization status.

The 28D searches included broader shells, core deletions, quadratic tail
rotations, reflected unions and finite integer programs. Their failure to
improve the current bound is not an optimality certificate.

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
the unchanged inequalities too. The four configurations cover
**111,131,795,527 unordered point pairs** in total.

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
