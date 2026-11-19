# 37D: 544,474 by packing six-coordinate radical supports

Keep the dense and signed octad layers of the [533,864-point QR
construction](../dim37-38-qr-merge/). Remove only the 144 signed roots
incident to the special QR coordinate. Keep all 2,520 ordinary roots.
Unlike the [536,440-point intermediate](../dim37-qr-shared-stars/),
do not use any of the star points. This permits a much larger layer
of signed six-coordinate supports.

For each of the 168 supports in `six_supports.txt`, use all 32
even-parity signings with entries ±2 on its six ordinary coordinates.
Set the special coordinate independently to ±2√2. Each of the
10,752 resulting points has squared norm 32. Two distinct ordinary
supports meet in at most two sites and each intersects every retained
octad in at most four sites. Thus if the special signs agree, the
ordinary product is at most eight (also true on identical support by
even-parity sign distance). If the special signs differ, the raw
product is at most −8 + 24 = 16. The same norm and product threshold
hold for the root and octad layers, whose supports avoid the special
coordinate. Against dense points of special magnitude m=1 or 4,
the raw product is at most 24 + 2m√2, strictly below half the norm
product; the replays check its square with exact integer bounds.

Both unit poles on the special coordinate also fit. They have
product −1 with each other, zero with the roots and octads, absolute
product exactly 1/2 with the normalized radical points, and squared
product at most 16/160 with the normalized dense points.

The set totals 180,224 dense points, 350,976 signed octad points,
2,520 ordinary roots, 10,752 radical points and 2 poles: **544,474**. This is
**10,610** above the original QR result and **8,034** above the
intermediate star construction.

The supports are explicit masks on the 36 ordinary columns. They
were selected from all 86,840 six-sets meeting the incumbent octads
in at most four coordinates, using a disjoint-triple packing search
with local rebuilding. A search restricted to orbits of an order-three
projective permutation produced a seed; dropping the symmetry
restriction and restarting from the improved packing gave 168 supports. No optimum is claimed.

Run from this directory:

```sh
python verify.py    # exact standard-library layer proof
python audit.py     # exact NumPy expansion of all 57,797,376 radical pairs
```

The full expanded check covers the radical-radical products.
Sparse/dense and cross-layer angles follow from the separate exact
bounds, and all unchanged interactions inherit the original proof.
