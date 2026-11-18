# Native Leech fibre obstructions in dimensions 26 and 27

The alternate 26D arrangement ties Takhanov at **201,192**. The 27D
comparison bound is **210,108**. For each of the two published configurations,
and for the alternate 26D arrangement, all **4,320** native Barnes–Wall
projections are exhausted for insertion while the configuration is fixed.
The result applies to arbitrary real remaining coordinates for those
projections. It leaves exchanges and points at other projections open.
[Source certificates](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1/certificates)

## The 26D alternate

Remove the 512 lifted Barnes–Wall midpoints and restore their 256 Leech
endpoints. In 64 other fibres, each with its own projection onto the
16 coordinates complementary to the octad, replace four signed Leech
vectors with eight signed midpoints. The two changes contribute −256
and +256. The source's 52 rational points still fit. The resulting
768-for-768 exchange and full coordinates are in `data/`. The exact
verification checks **20,239,009,836** unordered pairs.

## Fixed-projection proof

The octad is I = (0, 2, 4, 5, 6, 10, 11, 23) with zero-based indices.
Write an old point as (b, h), where b is its projection on the
complementary 16 coordinates J. For each candidate projection w of
squared norm 64, the proposed point (w, y) needs squared norm 64 in y.
Its compatibility constraints are

```text
h · y ≤ 64 − max { w · b : (b, h) is an old point }
```

Within the same fibre, the remaining coordinates have pairwise
nonpositive inner products, so the orthoplex bound is 2(d − 16).
Elsewhere, opposite inequalities with zero right sides force linear
equalities. The checker verifies a rational nullspace and certifies that
the norm on the reduced polytope stays below 64, using rational
endpoints, polygon vertices or exact dyadic subdivisions. The 26D
certificates follow using the integer core alone.

| Configuration | Orthoplex | Reduced dimension 1 or 2 | Dyadic subdivision | Norm impossible |
| --- | ---: | ---: | ---: | ---: |
| Published 26D | 1,024 | 3,196 | 0 | 100 |
| Alternate 26D | 1,088 | 3,232 | 0 | 0 |
| Published 27D | 0 | 4,078 | 210 | 32 |

## A triangle replacement in one 27D fibre

Choose three orthogonal body axes u₁, u₂, u₃ of norm 8 from one Leech
fibre and three orthogonal norm-8 tail axes v₁, v₂, v₃. On each directed
edge i → j of 1 → 2 → 3 → 1 take the four vectors with remaining
coordinates

```text
((e u_i + f u_j)/2, (e v_i − f v_j)/2),   e, f independently ±1.
```

Each edge gives two orthogonal lines. At the shared index of two edges,
the body and tail products cancel. The twelve new points consequently
occupy six orthogonal lines with both signs. The five untouched body
axes add ten points, giving 22 points in the fibre. The search over
this coordinate-axis family found no profitable single-fibre
replacement in the published 27D configuration. Coordinated changes
across multiple fibres remain possible.

## Replay

Use Python 3.10 or later with NumPy and SciPy. To replay the complete
alternate witness and its fixed-projection obstruction, run these
commands from this directory.

```sh
OPENBLAS_NUM_THREADS=1 python verify26.py data/coordinates26_201192_alternate.npz
OPENBLAS_NUM_THREADS=1 python certify_fibres.py data/coordinates26_201192_alternate.npz
```

For the published arrangements, download the two coordinate archives
at the pinned revision and check their SHA-256 digests against
`reports/`. Then run the following.

```sh
revision=dadcb2940f7c3c3bd478b0543931d3eef9f4bd65
url=https://raw.githubusercontent.com/k-nic/KissingNumbers/$revision/version1/certificates/data
curl -fL "$url/d26_coordinates.npz" -o /tmp/d26_coordinates.npz
curl -fL "$url/d27_coordinates.npz" -o /tmp/d27_coordinates.npz
OPENBLAS_NUM_THREADS=1 python certify_fibres.py /tmp/d26_coordinates.npz
OPENBLAS_NUM_THREADS=1 python certify_fibres.py /tmp/d27_coordinates.npz
python reconstruct26.py /tmp/d26_coordinates.npz data/exchange26.npz /tmp/rebuilt26.npz
```
