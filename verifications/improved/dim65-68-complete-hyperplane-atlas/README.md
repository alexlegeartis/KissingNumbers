# All 126 affine hyperplanes in 65–68 dimensions

The earlier [24-cell construction](../dim63-68-affine-hyperplane-heads/)
used 94 of the 126 affine hyperplanes of the six-bit body space.
Here all 126 carry heads simultaneously. Every hyperplane still carries
the complete [32,11,12] trace-quadratic fibre at each retained head; we do
not thin any of these fibres.

| Dimension | Before | New lower bound | Gain |
|---:|---:|---:|---:|
| 65 | 418,874,240 | 418,995,074 | 120,834 |
| 66 | 429,011,844 | 429,274,500 | 262,656 |
| 67 | 448,201,606 | 448,596,358 | 394,752 |
| 68 | 479,699,848 | 481,013,656 | 1,313,808 |

The algebraic rule has three ingredients.

1. The hyperplane label fixes one of 32 possible five-bit sign
   signatures on its intersection with another hyperplane. An
   intersection of nonparallel affine hyperplanes is an affine
   4-flat. The five parities are evaluated on the entire flat
   and on four independent coordinate cuts. The literal
   assignment in `hyperplane_heads.txt` makes the two signatures
   different on every such pair. Quadratic words do not alter
   these parities. Use 24 distinct quadratic phases modulo the
   same trace code for the vertices of the 24-cell on each support.
   The offset and exact cross-family checks are as in the earlier proof.
2. For each old sparse support, solve its finite conjunction of
   required shared eight-site parities independently. The checker
   validates all 7,812 head-to-sparse overlap guards and all 31,512
   sparse shared-parity edges. We keep the 64 extra pure-body
   supports and leave out the 128 four-tail supports that blocked
   eight of the added hyperplanes. In 65D this head family takes the
   place of the previous polar cells, and both poles are retained.
3. Add norm-squared-four points: for every two-tail mask use 32
   unordered body pairs from one nonzero XOR direction, different
   for different masks. Take all sixteen signs on each resulting
   four-site support. In 68D also take all sixteen signs on the
   four tail coordinates alone. The four-site supports intersect
   each other and all original affine two-flats in at most two
   positions. Their other cross-dot inequalities follow directly
   from their support sizes and the tail coordinates of the 24-cell.

The 68D head population is exactly `126*24*2048=6,193,152`;
in 65–67D retain two, four or six axial heads per hyperplane.
The four-site points contribute 512, 1,536 and 3,088 in 66–68D.

Run `python3 verify.py` on a BMI2 x86-64 machine with a C++20
compiler. The old witness is loaded from the adjoining pinned package;
all additional witnesses and exact counts are checked locally. No
numerical tolerances or external searches are used. See
[local fibre caps](LOCAL_BOUNDS.md) for the limitations of two ways
of making these individual fibres larger.
