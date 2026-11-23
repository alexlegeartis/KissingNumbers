# Tilted half-sign tails and coupled Γ72 contact profiles, 77–85D

| Dimension | Verified lower bound | Increase over the preceding repository bound |
|---:|---:|---:|
| 77 | 6,223,377,992 | 143,664 |
| 78 | 6,227,143,800 | 572,328 |
| 79 | 6,233,438,398 | 1,290,344 |
| 80 | 6,246,118,720 | 2,811,480 |
| 81 | 6,248,779,696 | 14,120,752 |
| 82 | 6,252,353,036 | 13,941,544 |
| 83 | 6,256,290,188 | 13,753,664 |
| 84 | 6,281,841,864 | 34,677,952 |
| 85 | 6,286,710,352 | 34,860,968 |

The construction uses Nebe’s [Γ72](https://arxiv.org/abs/1008.2862v3) as its
72-dimensional core. It retains its full norm-eight shell and the
previous pure and B-only norm-four tail layers. Some raw triple-root
classes are reallocated rather than all previous coordinates being fixed.

## Half-sign tails in 81–85D

Start with the [complete D_m rays](../dim81-85-gamma72-apn-complete-dm/)
and the [published norm-four tail layer](../dim81-85-gamma72-norm4-layer/).
New antipodal rays use entries ±1/2 on eight sites and zero elsewhere.
Use the first eight sites in 81–83D. In 84–85D, split the first twelve
sites into three four-sets and use each union of two sets. On each
support, fix the first sign positive and take the 64 even sign patterns.
Distinct rays have absolute dot at most one, also with every D_m ray.
Each uses a distinct free Γ72 phase in the three-colour palette. The
literal insertion order and phase assignments are in `phases.json`.
For the lower four dimensions the existing D/E rays and phases are unchanged.

A pure old norm-eight tail has dot at most four with a new tail by
Cauchy–Schwarz. The B-only norm-four tail vectors have L1 norm at most
four, hence dot at most two with a new tail, and their horizontal
difference has squared distance at least six. The old Γ72 core is
separated from each tilted packet by squared distance at least 6+2.

## Exact horizontal contact graph

Different ray colours already give horizontal squared separation at
least six. For equal colours and nonorthogonal tails, the only dangerous
horizontal separation is a singleton raw norm-eight difference in one
Leech block, with B-only shift the XOR of the two phase labels. Local
coordinates in the other two blocks must coincide. The full horizontal
raw norm is twelve, so the local norm in the changed block is also
unchanged. Thus 4+4+4, 4+8, and 6+6 can be optimized independently.

For **4+4+4** we first admit the signs certified by directed clipping,
then fill any class whose offending singleton images are absent from the
actual opposing packets. We further reallocate triples within independent
components. A component is contained in a fixed colour and A-label and
fixed local B-labels modulo the common shift subspace. Within each such
group there is an edge exactly when the tail directions are nonorthogonal
and exactly two local B fields agree. Those edges are checked again on the
selected set. Each selected row contributes sixteen points. The bounded
graph search is only a way to find a feasible set; no claim of solved
maximum independence is needed.

For **4+8** and **6+6**, the two occupied block positions and local A
labels are fixed across any conflict. In the changed block a horizontal
norm-eight difference occurs exactly when the inner product of the signed
local representatives is four, so when both signs are included the test
is absolute Gram product four. For a shifted norm-four class this is
automatic separation failure against the matching unchanged local axis:
the two distinct root classes have raw separation eight, for all signs.
The replay checks the resulting actual intersections between packets and
adds each surviving two-block row as a group of eight points.

## Replay

Run with ordinary assertions enabled, Python 3, NumPy, SymPy and a
little-endian C++17 compiler (`g++`):

```sh
OPENBLAS_NUM_THREADS=1 python3 verify77.py  # 77–80D
OPENBLAS_NUM_THREADS=1 python3 verify.py    # 81–85D
```

Existing proofs are pinned by SHA-256 and replayed. BLAS handles only
integer Gram products with a checked exact-in-binary64 bound. These are
constructive lower bounds, not optimality claims or alternatives to Nebe.
