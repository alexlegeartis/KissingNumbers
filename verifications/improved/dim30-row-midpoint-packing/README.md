# 248,506 points in dimension 30

An exchange adds 340 points and removes 250 from the [248,416-point parent](../dim30-reflected-certificate-unions/) at [revision fde77c3](https://github.com/alexlegeartis/KissingNumbers/tree/fde77c3e8b9124640ec65f8915337960029f30d4/verifications/improved/dim30-reflected-certificate-unions/). The resulting lower bound is **τ(30) ≥ 248,506**, a gain of 90.

## Construction

Split the coordinates into the octad I = {0, 2, 4, 5, 6, 10, 11, 23}, the remaining sixteen Leech coordinates J, and the six tail coordinates T. Use integer vectors of squared norm 128 throughout. The Leech roots whose I and J projections both have squared norm 64 form 135 classes. Within one class all combinations of 16 signed octad frame vectors and 32 signed complement frame vectors occur.

Take two roots (a, w₁, 0) and (a, w₂, 0) with w₁ and w₂ orthogonal. Their midpoint has squared norm 96. Append a tail t with two entries independently equal to ±4 and four zero entries. The resulting point (a, (w₁ + w₂)/2, t) has squared norm 128. Exchanging the roles of I and J gives a second midpoint family. Any Leech root other than the two endpoints has product at most 64 with each endpoint, hence also with their midpoint. This reduces the equatorial conflict calculation to at most two deletions per candidate.

The first family alone has 2,160 choices for a, 480 unordered orthogonal pairs per choice, and 60 tails. Among its 62,208,000 points, 80,434 were absent from the parent and conflicted with at most six of its points. The search also admitted other existing layer types. It solved exchanges in sets sharing a projection, recomputed conflicts against improved configurations, and considered signed permutations of those configurations.

To combine trades, give each candidate an insertion variable and each affected parent point a deletion variable. An insertion forces every conflicting parent deletion. Two candidates with product above 64 cannot both be inserted. Maximize insertions minus deletions, counting each distinct deleted parent once. At 248,496, several solutions achieved the same count with different points. Searching those alternatives exposed ten more net insertions. The supplied finite instance contains 1,403 distinct candidate points.

| Squared norm on I | Squared norm on J | Squared norm on T | Inserted | Removed |
| ---: | ---: | ---: | ---: | ---: |
| 32 | 64 | 32 | 132 | 132 |
| 64 | 0 | 64 | 36 | 48 |
| 64 | 32 | 32 | 144 | 30 |
| 64 | 64 | 0 | 28 | 40 |

## Verification and search replay

The coordinate verifier pins the parent and exchange by SHA-256, reconstructs the child, checks uniqueness and squared norm, and covers all **30,877,491,765** unordered pairs. The maximum inner product between distinct points is **64**, giving angular separation at least 60 degrees. A second calculation checks all 84,434,070 pairs involving an insertion directly in signed 64-bit integer arithmetic. See the [manifest](manifest.json) and the reports for [exact bounded-integer matrix products](verification.json) and [signed 64-bit arithmetic throughout](verification_int64.json).

The registry driver `final.py` reads the archive shape and checks it against the manifest. Run the geometric verifier from this directory with NumPy installed.

```sh
OPENBLAS_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 python verify.py --direct
```

To solve the frozen finite instance again, install OR-Tools as well as NumPy. The solver may produce a different optimum of the same size. The geometric verification above is independent of OR-Tools.

```sh
OPENBLAS_NUM_THREADS=1 python optimize.py --workers 4 --output /tmp/rebuilt30.npz
```

To enumerate fresh candidates on one of the 135 frames, the following two commands generate its row midpoints and optimize them against the original parent. Frames are numbered from zero by lexicographic order of the sorted octad vectors. Each frame can be searched independently for candidate generation. Trades on different frames can share deletions in the joint optimization.

```sh
OPENBLAS_NUM_THREADS=1 python enumerate_midpoints.py 23 /tmp/frame23.npz --family row
OPENBLAS_NUM_THREADS=1 python optimize.py --candidates /tmp/frame23.npz
```
