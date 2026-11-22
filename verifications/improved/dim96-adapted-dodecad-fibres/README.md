# Adapted affine sign fibres for 96D Golay dodecads

This construction proves **τ(96) ≥ 209,295,249,216**, an increase of
**827,308** over the preceding **209,294,421,908** bound. It uses the same
post–Reed–Solomon base **209,294,032,480**, but **replaces** the former
joint dodecad append of 389,428 by a larger append of **1,216,736**. The
two append layers must not be added together. The post–Reed–Solomon
base stock is retained in full.

As before, the new modes are monochromatic Golay dodecads using one D4
matching frame and twelve occupied rows of squared norm 8. Every
incumbent shifted mode with potential raw dot over 48 shares exactly
eight same-frame rows, and is itself monochromatic. The old axis and
sign codes are independent: the 4,096-word dodecacode and a
256-word shortened Hamming code. For each new axis word, the verifier
checks the projected old dodecacode at distance zero and one on these
eight rows. Every such nearby old axis word imposes a sign exclusion:
the new signs on the equal-axis rows must avoid the projected old sign
set. Otherwise a raw dot of 56 or 64 could occur. This test is applied
to the actual 2,970 older outer patterns.

## Adapt the sign code to each affine fibre

For each surviving new axis word, the old sign exclusions cut out an
affine binary set S. A nontrivial inherited sign projection on seven or
eight rows, when not full, is a single hyperplane. Avoiding it is an
affine equality; avoiding a full projection is impossible. The
verifier first recovers S exactly, with no added sign-code constraint.

In the homogeneous kernel of S, a change in a single coordinate is
possible precisely when its constraint column is zero. A change in two
coordinates is possible precisely when their constraint columns agree.
Thus an extra parity check needs distinct labels **only within a group
of equal columns**. The zero-column group also needs nonzero labels.
Give the zero group labels 1 through its size, and every other group
labels 0 through its size minus one, in increasing physical slot
order. Choose one coset of the XOR of labels on negative signs. Within
this coset, no one- or two-sign change is possible. The verifier
derives the groups by direct membership in S, independently of the
linear elimination used in discovery. This avoids the needless
constraints of prescribing twelve globally distinct syndrome labels.

For distinct new axis words in the same frame, the raw product is
8 times the shared equal-axis count minus 16 times the number of
opposite signs on those axes; every active projected sign pair is
verified exactly. Different D4 matching frames cost at most 4 per
common row. All cross-products with the dense, six-row triality, RS,
BCH and root stock use the same axis-independent bounds as the
preceding shifted dodecad geometry.

Run `python3 verify.py` from the repository with assertions enabled.
The witness is base64 of a zlib stream. Its little-endian schema is a
16-bit count of geometries and then records `(uint32 mask, uint8 tag,
uint16 count, uint32[count] words)`. The low 24 bits of each word
are the twelve two-bit axes in increasing physical row order; bits
24–27 give the selected adapted sign syndrome. All signs are
reconstructed by finite projections and parity checks, without a
numerical solver or a heuristic step in the verifier.
