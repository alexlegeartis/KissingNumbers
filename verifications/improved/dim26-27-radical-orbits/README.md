# Exact 26D and 27D improvements over PackingStar

Complete coordinate certificates establish the following kissing-number lower
bounds, compared with PackingStar's November 11 release at commit
`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`:

| Dimension | PackingStar | This certificate | Gain | Unordered pairs checked |
| --- | ---: | ---: | ---: | ---: |
| 26 | 201,180 | **201,188** | **8** | 20,238,205,078 |
| 27 | 209,862 | **209,884** | **22** | 22,025,541,786 |

These are lower bounds, without an optimality or exhaustive worldwide-priority
claim. Both improvements replace points only in the normal sector supported
on core octad `[0,2,4,5,6,10,11,23]` and the additional tail coordinates.

In 26D the normal sector grows from 380 to 388 points: 160 Leech roots, 128
midpoint points, 64 heads, four pure tails and 32 radical points. The source
sector has 168 roots, 128 midpoint points, 64 heads, four pure tails and 16
radical points. All heads are retained. The 160-for-152 exchange inserts 120
midpoint points, 32 radical points and eight restored Leech points; it removes
120 midpoint points, 16 radical points and 16 Leech points. The final code
has 201,156 rational and 32 radical points.

In 27D a complete orthogonal E8 head frame on one tail axis is replaced.
Jointly selecting compatible midpoint orbits grows the normal sector from
422 to **444** points: 144 Leech roots, 192 midpoint points, 96 heads and 12
pure tails. The source has 146 roots and 168 midpoint points, with the same
head and pure-tail counts. The 238-for-216 exchange inserts 32 heads, 176
midpoint points and 30 restored Leech points, and removes 32 heads, 152
midpoint points and 32 Leech points. All **209,884 points are rational**.
The other 209,440 points are retained exactly.

Coordinates are `(A + sqrt(3)*B)/2`. Every point has squared norm 32 and
every distinct inner product is at most 16. Equivalently, stored numerators
have squared norm 128 and inner product at most 64.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
```

Only Python and NumPy are required. The default proof checks both latest
records. It pins every supplied file, reconstructs each exact source exchange,
verifies all inserted-point interactions with separate coordinate-at-a-time
integer arithmetic, and checks **every unordered pair of each final code**.
Final validity assumes no unverified source geometry. Use `--count 201188`
or `--count 209884` to replay one dimension.

The full geometry checker is independent of discovery. Rational matrix
products are exact binary32 operations: integral inputs and every absolute
partial sum are bounded by 128, below the 2^24 exact-integer limit. Radical
pairs use integer lower and upper bounds for sqrt(3), obtained by an integer
square root; no floating tolerance is used. Six malformed-input controls
exercise norms, distinctness, both pair types and metadata.

`data/verification26_201188.json` and `data/verification27_209884.json`
record complete verification and coordinate hashes. `final.py` reads the
actual retained counts for the registry. `verify_geometry.py` is the
standalone full-pair checker invoked by `verify.py`.

`discovery_snapshot.py`, `discovery_normal_snapshot.py` and
`discovery_frames_snapshot.py` preserve the CP-SAT research scripts, with
their original workspace paths and research dependencies. They are unnecessary
for certificate replay. `data/model_orbits27.npz` preserves the 27D model
input; `data/discovery27.log` records the successful one-worker search with
seed 641, axis 0 and selected frame 33 among 275 frames. It reached a feasible
444-point sector after approximately 632 seconds. No solver optimality claim
is used for the 27D lower bound.

The earlier 201,182 and 201,184 certificates are retained; replay with
`python verify.py --count 201182` or `--count 201184`.
`data/discovery26.log`, `data/discovery26_201184.log` and
`data/discovery26_201188.log` preserve those searches.
`data/optimization26.log` records the optimum of the restricted fixed-head
E8 midpoint/root-orbit model. That restricted optimum is not a global
kissing-number upper bound and is unnecessary for the geometric proof.

The pinned source configurations are credited to
[PackingStar](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
