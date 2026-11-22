# Joint axis and sign rescue in dimension 96

The lower bound is **209,294,421,908**, an increase of **389,428** over the
post–Reed–Solomon bound **209,294,032,480**. All previous points are retained.
The append uses 1,619 new monochromatic Golay dodecad modes, containing
12,688 choices of the twelve projective axes. Their sign sets are explicit
subsets of binary Hamming-code cosets determined by the predecessor
intersections.

In each of the twelve occupied physical rows, the direction lies on one
of the four orthogonal lines of one D4 matching frame. The row has two
entries of magnitude 2, hence squared norm 8; the total squared norm is
96. If two twelve-row patterns share q rows, of which r have the same
frame, the raw product is at most 4(q+r) before using any phases.
Different frames contribute absolute value 4 per common row. In the
same frame, different projective axes are orthogonal; identical axes give
±8 according to their signs.

For every new geometry, the only incumbent shifted dodecad patterns
with q+r>12 share exactly eight rows and have the same frame throughout.
These nine incumbent patterns are all monochromatic, so the prior
dodecacode phase is zero and their slots are sorted in physical row
order. The new six-row triality and all other inherited layers use the
same axis-independent interfaces as the shifted dodecad proof. The
append does not assume the prior dodecacode data are merely an estimate:
the actual 2,970 older patterns are loaded by Git blob hash.

For a fixed candidate axis word, compare its eight shared slots with
the literal projection of the 4,096-word predecessor dodecacode. If the
axis words differ in at least two of those slots, the mixed cap 48
already holds. Otherwise, remove from the candidate sign set every
word that has a matching predecessor sign assignment on all equal axes.
After this removal at least one matched row has the opposite sign,
reducing the raw product 56 or 64 by 16. This condition admits more
than pure axis-code separation: the sign and axis labels have to be
checked jointly. The verifier applies this test to the zero error and
each possible one-axis error, against all relevant predecessors.

Among new axis modes of the same frame, the product is 8 times the
number of shared equal axes minus 16 times their sign disagreement
count; the verifier checks the exact necessary sign separation.
Within one identical axis word, every chosen sign set is a subset
of an affine [12,8,3] binary Hamming-code coset, so its cardinality
and minimum separation are also exact.

Run `python3 verify.py` from the repository with assertions enabled.
The one compressed witness `witness.b64` is base64 of a zlib stream.
The uncompressed little-endian schema is the 16-bit number of geometries,
then records `(uint32 mask, uint8 tag, uint16 axis_count, uint32 words[])`.
In each word the lower 24 bits hold the twelve two-bit axes in increasing
physical row order; the next four bits give the sign coset's syndrome
(XOR of one-based positions of negative signs). The explicit sign words
are recovered by the deterministic projection rule above, without a
search program, floating-point rounding, or optimizer.
