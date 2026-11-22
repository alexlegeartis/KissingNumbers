# Dimension 89: mixed five- through eight-row matching frames

**15,184,866,834 kissing points**, up **783,360** over the first mixed selection and **7,476,752** over the published
variable-puncture Golay configuration.

Keep the [signed row design](../dim88-89-signed-row-designs/) through
its four-row triality and single-row D4 layers (15,178,573,330
points). Replace its 885,760 matching-layer points by 6,293,504
points on the tagged patterns in `witness.txt`. The matching
coordinates and zero height are as in that construction.

| Active rows | Patterns | Rays per pattern | Rays |
|---:|---:|---:|---:|
| 5 | 2 | 1,024 | 2,048 |
| 6 | 859 | 4,096 | 3,518,464 |
| 7 | 3 | 8,192 | 24,576 |
| 8 | 61 | 45,056 | 2,748,416 |

Each tag specifies a matching of the four coordinates in its row.
The four orthogonal axes consist of two nonzero magnitudes √2 on
a matched pair. On a (w)-row pattern the squared norm is (4w).
The axis codes of lengths 5, 6, 7 and 8 have sizes 64, 128, 128,
352 and minimum distances 3, 3, 4, 4. The verifier generates the
projective five-coordinate code, a seven-coordinate additive
code and its puncture, and [Laaksonen and Östergård’s published 352-word code](https://www.cs.helsinki.fi/u/ahslaaks/codes/) (section 2, `codes2a.zip`). The verifier checks all distances in the supplied numerical codeword list.
Each takes all even row-sign words. An axis disagreement contributes
zero, while a sign reversal changes the dot by 8.

For two patterns let (a) count common rows and (b) their matching
tags on those rows. The dot is at most (2(a+b)) between new
patterns, and (2\sqrt{2}(a+b)) against the old Golay layer of
squared norm (8r). In either case it suffices to check
((a+b)^2\le wr). Each pattern with (w<8) avoids containing a
Steiner four-row support, giving at most three occupied rows against
the triality stock; for (w=8) four shared rows are allowed at equality.
Each new row contributes at most (2\sqrt{2}) against the dense and
triality alphabets, and the pure-height points are orthogonal.
All dimensions and sparse input identities are checked against
pinned siblings. Run `python verify.py` with assertions enabled.
