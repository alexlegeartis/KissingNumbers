# 38D: coordinated octad switches and signed stars

This certificate constructs **652,072** unit directions in dimension 38
with pairwise inner products at most (1/2). It adds 28 over the
[preceding 652,044-point certificate](../dim38-qr-star-trades/).

| Layer | Points |
| --- | ---: |
| Dense QR fibres (unchanged) | 253,952 |
| 3,088 full signed octads | 395,264 |
| Retained signed roots | 2,584 |
| Signed stars | 272 |
| **Total** | **652,072** |

The [witness](witness.json) specifies eight old/new octad replacements
at zero-based indices of the
[original support file](../dim37-38-qr-merge/supports_3088.txt). Each
eight-coordinate support has all 128 even-parity signings. The 17 star
templates have a center of magnitude 2 and four leaves of magnitude 1;
each uses both signs of the center and all eight odd-parity leaf
signings. Every two-coordinate root of magnitudes 2 is retained unless
its raw product with a star exceeds 4. Exactly 228 roots are removed.
The octads, stars and roots have squared norm 8.

The octad substitutions leave the octad count unchanged but create new
reciprocal incidences between star centers and leaves. Opposite
orientations of the same center pair can eliminate the same signed
roots, letting more stars share the deletion cost. In this witness the
stars exceed the removed roots by 44, versus 16 in the previous
certificate.

Run `python verify.py` here (Python 3.10+, standard library only). The
checker reconstructs the [unchanged QR dense layer](../dim38-qr-free-offsets/)
and checks all 465 pairs of dense fibres by exact code distances. It
checks every octad support pair, star pair and cross interaction using
integers. For the dense-to-sparse pairs a uniform magnitude bound
uses only the dense coordinate magnitudes and the one marked coordinate.
The computation neither uses a tolerance nor asserts optimality.

The independent replay, `python audit.py`, reconstructs the binary
syndrome from the full code span and checks every changed sparse layer
using explicit signed coordinates.
