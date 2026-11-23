# Axial weight-12 additions in 66D and 67D

These configurations keep every vector in the [preceding self-sparse
construction](../dim62-67-rm2-selfsparse/). The additions come from the
128 weight-12 words of its binary code \(C\) after translation by
cubic label 262160. Use amplitude 2 on each of the twelve body sites
and amplitude 4, with both signs, on one of the added axes. The
squared norm is \(12\cdot4+16=64\).

The four old support bodies are the weight-16 stock and three
weight-14 stocks, in cubic-shift order 0, 822310, 213640, 779738.
The maximum size of the intersection with each whole stock, for
each new twelve-set, has one of two exact profiles:
\((8,6,7,7)\) for 96 supports, \((8,6,6,6)\) for 32.
The cross product with any old sparse point is at most
four times the body overlap, plus 8 if the head axis lies in its
two-site head. This permits all 128 supports on both axes in 66D
and the latter 32 on all three axes in 67D. Against a dense body
the bound is 24, against a four-site body at most 32, and against
an axial pole at most 32.

Within a support, take the shortened extended-Hamming \([12,7,4]\)
kernel with columns \((1,j)\), \(0\le j<12\), in ascending order of
the body sites. On the three respective axes translate its bit
patterns by the XOR masks 0, 3 and 5. The translates have even weight,
are disjoint, and are at cross distance at least two. The verifier
also tests every same-support pair of signs and head choices
directly. Distinct twelve-sets meet in at most four sites, so their
body product is at most 16 and the head product is at most 16.

- 66D: 421,522,308 + 65,536 = **421,587,844**.
- 67D: 437,513,094 + 24,576 = **437,537,670**.

The verification is independent of a stored support catalogue. It
regenerates every relevant old and new stock by Gray traversal of
the quadratic polynomials, checks the complete intersections, and
replays the sign inequalities. The parent certificate supplies
compatibility within the unchanged old configuration. This package
does not claim compatibility with the separate 65D polar construction.

    python3 verify.py
