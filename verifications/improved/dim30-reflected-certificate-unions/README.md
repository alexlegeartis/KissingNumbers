# 248,416 points in dimension 30

The integral coordinate witness proves τ(30) ≥ 248,416. The construction adds
65 points to the [248,351-point contact-sum witness at revision 73981a8](https://github.com/alexlegeartis/KissingNumbers/tree/73981a8a0dc3310db3c01ff8bb592af9a74ef67d/verifications/improved/dim30-contact-sum-exchanges).

Let S be the pinned witness. Number coordinates from zero and let
I = {0, 2, 4, 5, 6, 10, 11, 23}, J = {0, …, 23} excluding I, and
T = {24, …, 29}. Start with S. At each step, make a new copy of S and
negate its coordinates in the indicated block. Form a bipartite graph
on the points unique to each of the two configurations, joining points
whose inner product exceeds 64. Their shared points remain compatible
with both sides. A maximum matching yields a minimum vertex cover and a
maximum compatible union. The reachable sets of alternating paths give
the exchange with fewest deletions from the current configuration.

| Reflected coordinates | Insertions | Deletions | New total |
| --- | ---: | ---: | ---: |
| J | 386 | 335 | 248,402 |
| I | 31 | 29 | 248,404 |
| T | 37 | 32 | 248,409 |
| All | 6 | 3 | 248,412 |

For the final four points, use either sign of
(−2, 2, 2, −6, −2, 2, −2, −2) on I and independently set coordinate
28 to ±8. Set the other coordinates to zero.

The frozen exchange inserts 456 points and removes 391. The checker
reconstructs the exchange, computes 113,173,500 unordered pairs involving an
inserted point in signed integer arithmetic, and checks all
30,855,130,320 unordered pairs of final points. Each vector has squared
norm 128 and the maximum distinct-point inner product is 64.
The replay is independent of the matching implementation.

From the repository root, use NumPy for verification and add SciPy
for reconstruction.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim30-reflected-certificate-unions/verify.py
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim30-reflected-certificate-unions/reconstruct.py \
  verifications/improved/dim30-contact-sum-exchanges/data/coordinates30_248351.npz /tmp/rebuilt30.npz
```

[Manifest](manifest.json) · [Exact verifier report](verification.json)

The repository registry reads the coordinate count through `final.py`.
