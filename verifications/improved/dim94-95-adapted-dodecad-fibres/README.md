# Adapted dodecad sign fibres in 94D and 95D

| Dimension | Previous lower bound | Appended points | New lower bound |
|---:|---:|---:|---:|
| 94 | 74,142,198,284 | 1,026,880 | **74,143,225,164** |
| 95 | 117,091,871,620 | 1,026,880 | **117,092,898,500** |

Every old point is retained. The new points all live on physical rows
2 through 23 of the common four-coordinate-row body, and therefore
extend by zero in the exceptional head directions in either dimension.
Each new mode uses a monochromatic Golay dodecad and one D4 matching
frame. All twelve rows have squared norm 8, and the total squared norm
is 96. A row's frame has four orthogonal projective axes; different
frames give absolute product 4 on a common row.

The 1,064 incumbent shifted patterns from the mask-3 sector are read
by Git blob hash. For each inserted mode, only three to five of those
can exceed 48 by the raw overlap bound: they each share eight rows,
all in the same frame, and are monochromatic with no old axis phase
shift. On those eight slots the verifier checks the zero and all
one-axis differences against the old projected dodecacode. When the
axes differ in two or more slots the product is already at most 48.
When they differ in zero or one, the old sign projection on the
equal-axis rows must not contain the new sign word's projection.

As in the [96D adapted proof](../dim96-adapted-dodecad-fibres/), the
near-axis old sign exclusions leave an affine sign fibre. Group its
coordinates by equality of their homogeneous parity-check columns.
Assign distinct small labels within each group, with nonzero labels in
the zero-column group, and retain one coset of the XOR of negative
signs' labels. The same syndrome cannot hold after a one- or two-sign
change. The verifier reconstructs the groups by membership in the
actual sign fibre and checks every active new-new projection directly.

The 85 existing short D4 matching modes are also protected. For such
a mode of weight w, with q common rows and r common matching frames,
the raw product is at most 2√2(q+r). The replay checks
(q+r)²≤12w, which is the exact squared cap for norms squared 96
and 4w. Cross-products with the remaining dense, grouped roots,
RS, BCH, and six-row triality layers inherit the same per-row
bounds as the shifted dodecad proof.

To replay, run `python3 verify.py` with assertions enabled.
The base64/zlib witness stores a little-endian 16-bit geometry count
and records `(uint32 mask, uint8 matching tag, uint16 axis count,
uint32 axis words[])`. Axis words store the twelve two-bit labels in
increasing physical row order, with the adapted sign syndrome in bits
24–27. Every sign is reconstructed exactly; the search optimizer is
not part of the certificate.
