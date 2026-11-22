# Coordinate-aligned dodecad fibres in 94D and 95D

| Dimension | Previous lower bound | Net gain | New lower bound |
|---:|---:|---:|---:|
| 94 | 74,143,250,508 | 3,663,872 | **74,146,914,380** |
| 95 | 117,092,923,844 | 3,663,872 | **117,096,587,716** |

This construction retains the entire post–Reed–Solomon stock of
74,142,198,284 and 117,091,871,620 points, respectively. It **replaces** the
previous 1,052,224-point adapted dodecad append by a larger 4,716,096-point
append. The two dodecad layers are not added together.

The geometry and old interfaces are the same as in the
[preceding adapted dodecad proof](../dim94-95-adapted-dodecad-fibres/).
In each D4 matching mode we occupy twelve of the common physical rows
2 through 23, with four possible projective axes in each row. For each
new axis word and every potentially active old mode, the replay checks
the old dodecacode projected to the eight shared same-frame rows. Only
zero- and one-axis differences can exceed the cap; for each such possible
difference the new sign on the equal-axis rows must avoid the actual
old sign projection. This leaves an affine set of signs. As before,
coordinates with equal homogeneous kernel columns need distinct added
XOR labels; the zero-column group uses nonzero labels. One syndrome
class then gives sign distance at least three.

The new idea is on the **axis** side of the same argument. Instead of
cutting all the locally admissible axis words by one fixed ordering of
the twelve-slot dodecacode, assign the dodecacode coordinate labels to
the physical rows separately in each geometry. In particular, copying
the old slot labels on one eight-row interface makes that projected
code identical; now translating it can move an entire fibre into a
favorable projection coset. Different assignments are scored by the
number of admissible sign words, not merely by the number of axes.
Finally select collectively compatible fibres in the exact sign-overlap
graph. The optimizer does not enter the verification.

Run `python3 verify.py` with assertions enabled. The base64/zlib witness
stores the actual geometry, packed axis words and adapted sign syndromes.
Its little-endian schema is a 16-bit geometry count followed by records
`(uint32 mask, uint8 tag, uint16 axis_count, uint32 words[])`. The lower
24 word bits are the twelve two-bit axis labels in increasing physical
row order; bits 24–27 specify the adapted sign syndrome. The verifier
decodes all signs, checks every active old and new overlap, preserves the
85 earlier short matching modes by the squared geometric cap, and
recomputes both lower bounds with integer arithmetic.
