# Reed–Solomon axis repacking in seven dimensions

| Dimension | Previous bound | Repacked six-row modes | Gain | New bound |
|---:|---:|---:|---:|---:|
| 89 | 15,184,866,834 | 859 | 989,568 | **15,185,856,402** |
| 91 | 25,922,303,504 | 850 | 979,200 | **25,923,282,704** |
| 92 | 34,512,238,096 | 850 | 979,200 | **34,513,217,296** |
| 93 | 47,397,139,984 | 850 | 979,200 | **47,398,119,184** |
| 94 | 74,141,118,156 | 832 | 958,464 | **74,142,076,620** |
| 95 | 117,090,791,492 | 832 | 958,464 | **117,091,749,956** |
| 96 | 209,292,462,656 | 1,334 | 1,536,768 | **209,293,999,424** |

The axis word in each six-row mode specifies one of four mutually orthogonal
lines in each of the six rows. The 32 even sign choices on those axes are
independent. The previous inner code had 128 axis words. Replace it in each
indicated mode by 164 axis words, keeping every other mode and layer unchanged.
This is a repacking, not a claim that the extra 36 axes can be appended without
replacing the old axis code.

## The inner code

Evaluate every cubic over `F_5` at zero through four, with the leading
coefficient as the sixth coordinate. These 625 words are the extended
`[6,4,3]_5` Reed–Solomon code: if two words have the same last coordinate,
the difference polynomial has degree at most two, and otherwise it has
at most three finite roots. Omit one symbol in each position, specifically
`(0,0,0,0,4,4)`, and relabel each remaining four-symbol alphabet to
`{0,1,2,3}`. Exactly 164 words survive, with minimum Hamming distance 3.
The size is also suggested by averaging over omitted-symbol tuples:
`625*(4/5)^6 = 163.84`. Thus the code uses 5,248 rays per mode, up from
4,096, a gain of 1,152 each. A quaternary code of length 6, size 164 and
distance 3 is a known bound; see [Bogdanova–Brouwer–Kapralov–Östergård
(2001)](https://doi.org/10.1023/A:1011275112159) and
[Andries Brouwer’s table](https://aeb.win.tue.nl/codes/quaternary-1.html).
The novelty here is the application to the seven kissing configurations.

## The interfaces

In one mode, distinct axes differ in at least three positions, so at most
three squared-norm-4 row vectors can contribute 4 to the dot product.
With equal axes, at least two sign positions differ, so the dot is at most
8. Either bound is at most 12, the kissing threshold at squared norm 24.

For the D4 matching modes in 89D, 91–93D and 94–95D, the pinned preceding
proofs bound every cross-mode product by the absolute maxima of the
occupied row frames regardless of actual axis labels. This includes
interfaces with the old Golay, dense, grouped roots and flat layers, and
with D4 matching modes of other row counts. In 94D and 95D, the preceding
mixed-triality and complementary hexary geometry also has an
axis-independent cross-mode bound and can be repacked: 473 three-frame
patterns plus 304 complementary even/odd patterns and 55 D4 matching
patterns.

In 96D, the preceding affine and partial fibres do need genuine
axis-projection checks. Across two triality modes with q shared rows and r
of them in the same frame, the axis-independent bound is 2(q+r). The
verifier checks all 1,771 mode labels, identifies the 317 pairs with
q+r>6, and uses the 164-word replacement only in the 1,334 full modes
incident to none of these pairs. The other 437 modes, including all modes
with affine or partial axes, are unchanged. The cross products with
all earlier dense, D4, grouped root, and RS/BCH layers already use row
maxima, not the old quaternary code.

Run `python3 verify.py` with assertions enabled in the repository.
It checks the Reed–Solomon list and all pairwise distances, pins the exact
preceding geometry and witness files, recounts the eligible modes and
reconstructs every bound. No search or solver results are needed.
