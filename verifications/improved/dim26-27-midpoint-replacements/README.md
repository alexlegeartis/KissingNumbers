# Shared midpoint replacements in dimensions 26 and 27

The current-bound column follows the [29 October joint-refill continuation](../dim27-31-joint-refills/).

Retained construction: the dated benchmark, construction and gain columns below are historical; the current-bound column uses the [26D asymmetric exchange](../dim26-30-asymmetric-exchanges/) and [27D quadratic module](../dim27-quadratic-normal-module/).

Checkpoint: **25 October 2026**.

| Dimension | Previous benchmark | This construction | Gain | Current repository bound (29 October 2026) |
|---|---:|---:|---:|---:|
| 26 | 199,868 | **200,918** | 1,050 | 201,136 |
| 27 | 201,589 | **204,636** | 3,047 | 209,775 |

The comparison benchmarks and their methods are documented by
[Takhanov and Yun](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf).
This construction starts from the standard Leech shell. The earlier coset-triangle
constructions remain in [`../dim26-27-iota-triangles/`](../dim26-27-iota-triangles/).

## Construction

All points have squared norm 32 and distinct points have inner product at most 16.
For Leech minimal vectors u,v with inner product 16, their midpoint has squared
norm 24. Adding a tail of squared norm 8 gives a new kissing point. Removing
both endpoints clears its products with every remaining equator point, since
each such product is the average of two products at most 16.

Two matchings share four deleted endpoints and provide eight inserted points
on orthogonal tail directions. Compatible groups of these replacements give
the following counts:

| Component | 26D | 27D |
|---|---:|---:|
| Leech shell | 196,560 | 196,560 |
| Deleted equator points | -4,352 | -8,064 |
| Inserted midpoint points | 8,704 | 16,128 |
| Pure tail points | 6 | 12 |
| Total | **200,918** | **204,636** |

The pure tails form a hexagon in 26D and a cuboctahedron in 27D. Cauchy-Schwarz
bounds all pure-tail/inserted products by 16. Every inserted/inserted product
is checked directly with integer arithmetic. The family is open; these achieved
lower bounds do not establish optimality of the configurations or the search.

## Reproduce

```bash
python final.py
python verify.py
```

Python 3.9+ and NumPy suffice. The wrapper checks a complete SHA-256 manifest,
validates coordinate ranges before narrowing integer types, and replays both
independent source checks in a temporary directory. Each check regenerates
the Golay Leech shell, verifies deleted-coordinate membership and midpoint
identities, and checks every inserted/inserted pair exactly. Retained equator
pairs use the standard kissing property of the Leech minimal shell. No external
cache, optimizer, coordinate download or search replay is needed.

`source/README.md` gives the full geometric argument and search provenance.
`source/verification.json` preserves the original recorded runs. `source/record26.json`
and `source/record27.json` describe the two coordinate NPZ files. `final.py` derives
the counts, `verify.py` runs the proofs, and `data/SHA256SUMS.txt` covers every
package file except itself.
