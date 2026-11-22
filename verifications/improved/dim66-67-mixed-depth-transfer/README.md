# Mixed-amplitude clones of equatorial owners in 66D and 67D

| Dimension | Prior lower bound | Mixed-clone lower bound | Gain |
|---:|---:|---:|---:|
| 66 | 342,390,470 | **342,390,951** | +481 |
| 67 | 355,196,300 | **355,203,873** | +7,573 |

The same ordinary-coordinate relabeling used in the 65D transfer is used
here to match BCH field labels to the physical affine labels of the cyclic
code. Relabel all BCH-based dense equator, Kerdock cap and GF(8) layers
together. Keep the old sparse stock unchanged: every earlier comparison to
it depends only on amplitudes, weighted support size, and tail occupancy.

All selected owners belong to both the BCH radical and the cyclic code.
Removing an equatorial owner permits cloning it to several signed marked
ordinary axes. A clone has one magnitude 11 (shallow) or 14 (deep), all other
ordinary magnitudes 2, and zero tail; normalize it to squared norm 64. The
sign at its own mark reverses. Each mark uses one amplitude, and every owner
has at most one shallow clone. The verifier gives exact rational bounds
for all comparisons on distinct marks and checks each same-mark pair.

For an unmodified ordinary W16 support, reversal on the mark guarantees
one mismatch by even cyclic parity. For modified or added supports,
the verifier reads the literal pinned sign checks and affine shifts for
each dimension. A W16 support with a nonzero tail has fourteen ordinary
positions, so even with the heavy coordinate and complete sign agreement
its product is at most 2(h+26) times the clone normalization. The bound
also covers all ordinary supports with the mark and at least one mismatch.
Against a Kerdock cap, membership of its ordinary owner in the exterior
BCH coset supplies distance at least fourteen before the heavy sign flip,
and its ordinary scale is sqrt(37)/8. The same distance also applies to
the retained GF(8) spikes; those are on different axes than the new marked
clones.

In 66D the replacement inserts 2,875 clones and retires 1,418 equatorial
owners and 976 GF(8) spikes. In 67D it inserts 14,777 clones and retires 6,180
equatorial owners and all 1,024 GF(8) spikes. No other points are removed.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with Python 3 and NumPy.
