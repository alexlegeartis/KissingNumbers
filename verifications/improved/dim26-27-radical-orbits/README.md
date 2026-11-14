# 26D radical-orbit continuation: 201,182

The complete coordinate certificate proves **K(26) >= 201,182**, two above
PackingStar's 201,180-point November 11 release at commit
`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`. This is a lower bound; no
optimality or exhaustive worldwide-priority claim is made. The 27D search
continues separately.

The construction reroutes sixteen four-point midpoint orbits inside the
ten-dimensional normal sector. It removes 64 old midpoint points, inserts
64 new ones, then replaces two Leech points with a four-point orbit over
Q(sqrt(3)). All head points are retained. The normal sector grows from 380
to 382 points. The full code has 201,162 rational and 20 radical points.

Coordinates are `(A + sqrt(3)*B)/2`. Every point has squared norm 32 and
every distinct inner product is at most 16. Equivalently the stored
numerators have squared norm 128 and inner product at most 64.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python negative_controls.py
```

Only Python and NumPy are required. The default proof pins every supplied
file, reconstructs the 68-for-66 exchange, verifies all inserted-point
interactions with separate coordinate-at-a-time integer arithmetic, and
checks **all 20,236,997,971 unordered pairs of the final configuration**.
Final validity assumes no unverified source geometry.

The full geometry checker is independent of discovery. Rational matrix
products are exact binary32 operations: integral inputs and every absolute
partial sum are bounded by 128, below the 2^24 exact-integer limit. Radical
pairs use integer lower and upper bounds for sqrt(3), obtained by an integer
square root; no floating tolerance is used. Six malformed-input controls
exercise norms, distinctness, both pair types and metadata.

`data/verification26.json` records the complete verification and coordinate
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
