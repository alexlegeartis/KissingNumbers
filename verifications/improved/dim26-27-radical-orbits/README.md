# 26D radical-orbit continuation: 201,188

The complete coordinate certificate proves **K(26) >= 201,188**, eight above
PackingStar's 201,180-point November 11 release at commit
`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`. This is a lower bound; no
optimality or exhaustive worldwide-priority claim is made. The 27D search
continues separately.

The normal sector grows from 380 to 388 points: 160 Leech roots, 128
midpoint points, 64 heads, four pure tails and 32 radical points. The pinned
source sector has 168 roots, 128 midpoint points, 64 heads, four pure tails
and 16 radical points. All head points are retained. The full 160-for-152
exchange inserts 120 new midpoint points, 32 radical points and eight Leech
points; it removes 120 midpoint points, 16 radical points and 16 Leech points.
The full code has 201,156 rational and 32 radical points.

Coordinates are `(A + sqrt(3)*B)/2`. Every point has squared norm 32 and
every distinct inner product is at most 16. Equivalently the stored
numerators have squared norm 128 and inner product at most 64.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
```

Only Python and NumPy are required. The default proof pins every supplied
file, reconstructs the 160-for-152 exchange, verifies all inserted-point
interactions with separate coordinate-at-a-time integer arithmetic, and
checks **all 20,238,205,078 unordered pairs of the final configuration**.
Final validity assumes no unverified source geometry.

The full geometry checker is independent of discovery. Rational matrix
products are exact binary32 operations: integral inputs and every absolute
partial sum are bounded by 128, below the 2^24 exact-integer limit. Radical
pairs use integer lower and upper bounds for sqrt(3), obtained by an integer
square root; no floating tolerance is used. Six malformed-input controls
exercise norms, distinctness, both pair types and metadata.

`data/verification26_201188.json` records the complete verification and coordinate
hash. `data/discovery26.log` and `discovery_snapshot.py` preserve the actual
CP-SAT discovery. The latter is a research snapshot using the original
workspace paths; it is unnecessary for certificate replay. CP-SAT's
`OPTIMAL` status here means the feasibility model found a solution, not an
upper bound on kissing numbers. The candidate files alone and the two
verifiers provide the proof.

The pinned source is credited to
[PackingStar](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).

`final.py` reads the actual retained coordinate count for the registry.
`verify_geometry.py` is the standalone full-pair checker invoked by `verify.py`.

The earlier 201,182-point certificate is retained. Replay it with
`python verify.py --count 201182`. `discovery_normal_snapshot.py` and
`data/discovery26_201184.log` preserve the stronger normal-sector search.

The 201,184 checkpoint is also retained (`python verify.py --count 201184`).
`data/discovery26_201188.log` records the final feasibility search;
`data/optimization26.log` records the reported optimum of the fixed-head
E8 midpoint/root-orbit model. That restricted upper bound is not used as a
premise for this kissing-number lower bound.
