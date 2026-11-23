# 36D: multicell Coxeter–Todd fibre reassignment

The construction gives 506,376 points of squared norm 12, with distinct
pair inner products at most 6. From the pinned 502,929-point stock,
insert 13,677 triple-root points, remove 6,630 old triple-root points and
6,057 two-block points, then add 2,457 reassigned two-block points.
This improves the previous 504,365-point witness by 2,011. The coordinates
are determined by the witness and the existing sibling certificates.

Run the exact verifier with Python 3 and NumPy (without `-O`):
```sh
OPENBLAS_NUM_THREADS=1 python3 verify.py
```

An independent verifier expands the literal 36-coordinate rows. It first
certifies 505,989 points against the pinned stock, then checks every
changed row against the entire final set at each step. Unchanged rows
inherit the preceding exact verification. The last step checks
2,106,017,784 inner products.

```sh
OPENBLAS_NUM_THREADS=1 python3 verify_independent.py independent505989.json
OPENBLAS_NUM_THREADS=1 python3 verify_independent_delta.py independent505989.json independent506079.json
OPENBLAS_NUM_THREADS=1 python3 verify_independent_delta.py independent506079.json witness.json
```

SHA-256 of `witness.json`:
`c14ba2da244f718d2fc329899e657a9e30e5c3c4abdc6205a6d9f5f3727c3109`.
