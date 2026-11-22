# 96-dimensional partial triality fibres

This extension proves **τ(96) ≥ 209,292,462,656**, keeping all **209,292,294,976**
points in the affine-phase construction and adding **167,680**. Its sixty
six-row geometries carry 5,240 explicitly listed quaternary axis words in
total, each with all 32 even flips of its six signs. Unlike the preceding
append, the axis sets may be different sizes and may independently permute
the four labels at each row.

In each of the 24 physical rows of four coordinates, use one of the three
orthogonal line-frames in the 24-cell: the ordinary coordinate axes, the
even-parity sign-flat axes, and the odd-parity sign-flat axes. In the integer
normalization each vector has norm squared 4 per occupied row. Different
lines of one frame are orthogonal; across distinct frames the product has
absolute value 2. When two geometries have q common rows and r matched
frames among these rows, a Hamming distance d between axis words projected
to the matched frames bounds their raw dot by 2(q+r)−4d. The target is 12.
The verifier enumerates the two projected sets directly for every potentially
active pair; it does not assume that a partial fibre is a linear code.

Within a geometry the explicit axis words have quaternary distance at least
three. When the axes agree and the signs differ, their even sign parity
forces at least two changes. The raw product is at most 12 in both cases.
All compatible configurations are automatically without duplicate vectors
since their squared norm is 24.

The mixed products with older dense, shifted Golay and sparse RS/BCH layers
use the same geometric envelopes as in the preceding proof. In particular
the new per-row dot against a shifted D4 row (two entries of magnitude 2)
is at most 4, so the six rows cost at most 24, exactly the mixed cap with
a squared-norm-96 point. The remaining support constraints are replayed against
the actual ordinary/paired BCH masks and physical RS groups. All 1,711 earlier
six-row modes (1,170 mixed, 438 pure hexary, and 103 affine) are retained.
The predecessor input, witness and count bytes are pinned by cryptographic digest.

Every line of witness.txt contains disjoint axial/even/odd frame masks, the
number of axis words and those words. Each word is packed in two-bit digits
in increasing physical row order. Run `python3 verify.py` (assertions enabled).
No solver, approximate distances or heuristic search assumptions are needed
to reproduce the bound.
