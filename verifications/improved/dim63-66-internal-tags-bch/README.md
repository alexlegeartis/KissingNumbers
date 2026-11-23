# Internal weight-12 tags and BCH polar cells in dimensions 63–66

Start with the [cubic self-derived shells](../dim62-67-rm2-selfsparse/).
The 66D calculation also retains the [external axial weight-12
additions](../dim66-67-rm2-axial-weight12/). These older certificates
supply the counts and validity of the unchanged points; the new replay
reconstructs every old support relevant to the added coordinates.

## Amplitude-four tags on body coordinates

Each new vector uses a twelve-site body of amplitude 2 and one
amplitude-4 tag at a different site *among the original 64 coordinates*.
The squared norm is 64. On the twelve body positions, ordered by index,
let the sign bits range over the kernel of the parity columns `16 | i`
for `i=0,...,11`. For a fixed body with more than one tag choice, use
distinct translates of that same shortened extended-Hamming code by
masks 0, 3, 5, 9 in the listed order. Each has 128 words; add either tag
sign. On one body the within-tag distances are at least four and
different tag palettes have distance at least two. Distinct bodies are
tested using worst signed overlap. The verifier also checks the
weight-16, weight-14, and external-tag weight-12 old stocks as applicable.
An old four-site body is an affine 2-flat. A possible violation with
the heavy site would require that site to be the XOR of three body
positions; all these equalities are excluded literally.

All 256 listed tags work in 64D; 64 of them work against the stronger
retained 66D stock. For 63D restrict to the previous section, coordinate
0 equal to coordinate 20. Sign choices surviving this equality are
counted directly. The exact additions are 39,424 in 63D, 65,536 in 64D,
and 16,384 in 66D.

## The 65D polar replacement

Use the 64D construction above as the equator and keep its two poles.
A coordinate-permuted `[64,18,22]` extended binary BCH code supplies
three of its four cubic cells, each a translate of its 16-dimensional
quadratic kernel. The witness gives the 16 truth-table generators and
the three shift words in the actual 64-coordinate order. On the positive
side append (+sqrt{17}); on the negative side shift all three cells
by the fourth literal word in the witness and append (-sqrt{17}).
There are 196,608 words per side. This replaces, and is not added to,
the earlier alternating-checksum polar layer.

The exact checker regenerates all 65,536 quadratic-kernel words. It
checks minimum 22 on each side and minimum 8 across sides, and
minimum 14 against every old dense quadratic coset. All polar
vectors have squared norm 81, so these distance thresholds imply
each required normalized inner-product bound. Every equatorial
sparse vector, including the new internal tags, has L1 norm at most
32, which also suffices.

| Dimension | Verified lower bound |
| --- | ---: |
| 63 | 214,125,316 |
| 64 | 413,592,448 |
| 65 | 413,985,666 |
| 66 | 421,604,228 |

Replay: `python3 verify.py` (requires a C++17 compiler).
