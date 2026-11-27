# 232,236 kissing points in dimension 29

**Six above [Rustem Takhanov's 232,230-point v0.2 configuration](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2), and 79 above [our previous 232,157-point witness](../dim29-even-grid-exchange/).** This is an extension of Takhanov's construction, not an independent replacement of its scaffold. It retains the entire equator and every original midpoint body.

## A two-colour switch

The raw squared norm is 128 and all distinct points have raw inner product at most 64. Write an old midpoint lift as `(b,t)`, with a body `b` of squared norm 96 and a 5-coordinate tail `t` of squared norm 32. A paired body carries both `(b,t)` and `(b,-t)`. The colours are unoriented tail lines.

Take the two orthogonal colours represented by

```
A = (0,4,0,-4,0)
B = (4,0,0,0,-4).
```

On paired bodies of these two colours, join opposite colours whenever the body dot product exceeds 32. Such a pair cannot become the same colour, since the maximum tail dot would be 32. It follows that flipping the colour of one endpoint requires flipping the other. Conversely, swapping the colours on a union of connected components preserves every constraint inside the two-colour bank: equal old colours had body dot at most 32, and edges above 32 between two colours are not split.

Take the singleton specified by `TARGET` in [verify.py](verify.py). Its own tail is `(4,0,0,0,4)`. The proposed new pair uses colour B instead. In the B bank, use all bodies blocking that pair as seeds, and swap the two colours throughout their forced components. This switches **3,792 existing antipodal fibres**, 1,896 of each colour, with no point loss. The switch frees the first singleton and seven more individually compatible antipodal replacements. Five of those seven can coexist; the remaining conflicts among the four competing B-colour bodies form a 4-cycle. Replacing the first singleton and the five compatible additional ones gives a net gain of six. Compatibility with all other colours and all other layers is checked directly. The earlier [single-fibre obstruction](../../../research/dim29_tail_recombination_20261127/) holds with every other fibre fixed; this coordinated change escapes that hypothesis.

## Exact replay

Download Takhanov's pinned [`d29_coordinates.npz`](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d29_coordinates.npz). The verifier checks both the distributed file SHA-256 and the coordinate-set SHA-256 before reconstruction. Takhanov's [published replay](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/) verifies all pairs of that stock; unchanged pairs are inherited from it.

```sh
OPENBLAS_NUM_THREADS=1 python verify.py d29_coordinates.npz --output coordinates29_232231.npz
```

Only NumPy is required. The optional `--full` additionally checks all retained–retained pairs and was replayed for this release. The verifier computes the forced components from the coordinates, keeps 224,640 original points and verifies 7,596 replacement points. It checks every new–retained pair and every new–new pair; both maxima are exactly 64. The binary32 matrix products are exact here: products are integers and Cauchy–Schwarz bounds each partial absolute sum by 128. The final sorted little-endian int16 coordinate-set SHA-256 is `81d5f3043438cbf3f8f3521033ecaa5de44a50e774b4610caee8ce54d2b234c6`.
