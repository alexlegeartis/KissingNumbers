# Barnes–Wall-type equality sections as Γ72 tails

| Dimension | Verified lower bound | Gain over preceding half-sign package |
|---:|---:|---:|
| 82 | 6,254,779,152 | 2,426,116 |
| 83 | 6,265,203,518 | 8,913,330 |
| 85 | 6,292,777,568 | 6,067,216 |

A new tail is cut from a sixteen-coordinate integer stock, without
projecting or approximately rescaling points. Label the coordinates by
the binary four-space. Take the antipodal roots with two entries of
absolute value two, together with the antipodal vectors of eight
entries of absolute value one on each affine hyperplane, with even
sign parity. Their squared norm is eight and distinct antipodal rays
have absolute inner product at most four. The verifier checks this
directly for each selected equality section. Divide its vectors by two
to get tails of squared norm two.

Each section sets designated coordinates to zero and requires
coordinates in each nonzero block to agree. The blocks give ambient
tail rank at most 10, 11 or 13 exactly as required. The sections yield
168, 219 and 453 antipodal rays, respectively. The construction uses
the first two sections in full; for 85D it selects 384 of 453 rays
and fills the entire three-colour phase palette. Literal rays and
phases are in `phases.json`. This is a replacement of the prior tail
stock, not an additive union with it.

Attach the tail rays to the 72-dimensional Nebe core. Reuse the exact
contact-graph method from the [previous package](../dim81-85-gamma72-half-sign/).
For each ray include both antipodal tails, with raw norm-twelve
horizontal contacts clipped and reallocated in the three separate
profiles. Include its pair of pure tails (twice the tilted tail alone)
and the entire Γ72 norm-eight core. By the Gram check these pure tails
are mutually compatible; they are compatible with all mixed tails
because their inner products are twice those of the norm-two rays.
The depth-six phase proof separates the core and all differently
coloured rays. The same-colour singleton conflicts are exactly replayed.

Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` with the same dependencies
as the preceding package. This still uses Nebe’s lattice.
