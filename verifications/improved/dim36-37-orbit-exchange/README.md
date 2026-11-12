# Exact orbit exchanges in dimensions 36 and 37

The compact certificate gives **501,822 points in dimension 36** and
**505,638 points in dimension 37**, increasing each preceding repository
construction by **nine points**.

| Dimension | Previous repository bound | New bound | Gain |
|---|---:|---:|---:|
| 36 | 501,813 | **501,822** | 9 |
| 37 | 505,629 | **505,638** | 9 |

Three three-point phase orbits are removed and six are added. The 403,704-point
base, all 90 existing triangle switches, and all 636 three-point latitude-body
orbits are retained. The preservation of the latitude layer is checked directly;
the 37D improvement does not rely on a monotonicity argument.

## Replay

With Python 3.9+, NumPy, and SymPy, from this directory run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python rejection_controls.py
python final.py
```

The default commands are read-only. `verify.py --report /path/to/report.json`
and `rejection_controls.py --report /path/to/controls.json` optionally write
fresh receipts. Do not use Python `-O`; the checker explicitly rejects it.
`final.py` derives the two counts from the orbit lists and pinned preceding
certificate; it is a count driver, not the geometric proof.

The complete replay uses the adjacent
[`dim36-37-coxeter-todd`](../dim36-37-coxeter-todd/) package. Seven exact source
files are pinned by SHA-256 in both the certificate and verifier. The proof
requires no network, optimizer, external coordinates, candidate pool or GPU.

## Exact construction

Let G be the retained integral Coxeter–Todd Gram matrix. Let Phi be its stated
order-three isometry, and let x_a be the phase-zero root for label a. An orbit
(a,b,c,u,v) represents the three points

```text
(x_a Phi^t, x_b Phi^(u+t), x_c Phi^(v+t)),  t = 0,1,2,
```

with phases taken modulo three. Every point has squared norm 12 in
diag(G,G,G). The compact exchange removes

```text
(23,77,23,u,1),  u = 0,1,2,
```

and adds

```text
(23,105,23,u,1),  u = 0,1,2,
(38,77,12,u,2),   u = 0,1,2.
```

The final equator therefore has

```text
403704 + 3*(32703 - 3 + 6) = 501822
```

points. The unchanged 1,908 latitude bodies have squared norm eight and are
used at both heights +2 and −2, adding 3,816 points in dimension 37.

## What the verifier proves

First, the verifier copies the pinned preceding triangle-switch proof into a
temporary directory and freshly replays it, including its exact completion
proof. This establishes every pair wholly inside the retained baseline.

The new checker then independently reconstructs the full current base from a
fresh rational enumeration of the norm-four, norm-six and norm-eight short
shells. It independently decodes each of the 90 triangle IDs and checks the
explicit switch record against that decoding. All final phase points are
distinct, and their norms are computed exactly.

Every pair involving an added point is checked by direct integer-coordinate
products. The new proof covers:

| Pair category | Ordered products | Maximum | Allowed maximum |
|---|---:|---:|---:|
| Added points against the entire current base | 7,266,672 | 6 | 6 |
| Added points against every final phase point, excluding identical points | 1,766,106 | 6 | 6 |
| Added points against every latitude body | 34,344 | 5 | 6 |

The verifier proves an explicit int64 overflow bound. No floating-point angle
tolerance or sampled pair test is used. Distinct norm-twelve vectors with
product at most six become unit vectors with product at most one half after
normalization, giving the stated kissing-number lower bounds.

The corruption controls include restoring a deleted orbit while correctly
updating the claimed counts. Its exact product with an added point exceeds six,
so the geometric checker rejects it. Duplicate additions, inflated counts,
modified source hashes and disabled assertions are also rejected.

## Scope and files

This is an explicit improvement over the pinned repository witnesses. It does
not assert global optimality or an exhaustive worldwide-priority survey. All
unchanged geometry is proved by the freshly replayed preceding certificate.
There is no inherited unenumerated shell population or spherical-design premise:
the short lattice shells used here are enumerated exactly.

- `certificate.json`: the compact nine-orbit exchange and pinned dependencies.
- `verify.py`: the complete baseline-plus-delta proof.
- `verification.json`: a passing full replay receipt.
- `rejection_controls.py`: the positive delta replay and five rejection controls.
- `rejection_controls.json`: the control receipt.
- `final.py`: the lightweight count driver.
- `SHA256.json`: hashes of the package files, excluding the manifest itself.

The bounded discovery code and negative search results are recorded separately
under `research/dim32plus_continue_20261113/lowdim/`; none is required for replay.
