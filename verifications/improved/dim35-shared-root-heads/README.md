# 415,336 points in dimension 35

This construction increases the preceding independently verified
415,334-point repository bound by two. It uses the 2,202 unchanged
octad supports of the inherited Lin-based construction, reassigning their
even/odd signs where required. All vector data are exact.

* Dense layer: the 131,072 binary signs of the given affine ([32,17,8])
  code, embedded in all coordinates other than 0, 3 and 18
  (zero-based). Norm squared is 32.
* Signed octads: all 128 sign words on each inherited octad of the
  prescribed parity, with nonzero magnitudes 2. Norm squared is 32.
* Root and star layer: retain 2,246 of the 2,380 usual D35 roots with
  nonzero magnitudes 4, and add the 150 explicitly listed signed stars
  with magnitudes ((4,2,2,2,2)), all norm squared 32.
* Heads: 12 explicit signed vectors with two magnitudes 3 and 33 unit
  magnitudes, norm squared 51. A head has its two heavy coordinates
  among 0, 3 and 18.

The stars block 132 signed roots. The heads block six more signed roots,
four of which are already blocked by the stars. Thus count the union of
blocked roots once: 132 + 6 - 4 = 134, and
131072 + 2202*128 + 2246 + 150 + 12 = 415336.
Use normalized products: an integer raw inner product d is permitted iff
d <= 0 or 4*d*d <= N1*N2, where N1, N2 are the raw squared norms.

Replay:
    python verify_shared.py
    python audit_coordinates.py

The first checker verifies all point categories and derives valid octad
parities by integer maxima. The second independently builds ordinary
coordinates for every inherited vector and checks all products involving
the changed stars and heads, without importing either the first verifier
or the discovery code. Its small-integer floating-point dot kernels are
exact since every individual product and partial sum is representable.
