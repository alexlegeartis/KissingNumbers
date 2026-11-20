# 38D: coordinated octad switches and signed stars

This certificate constructs **652,076** unit directions in dimension 38
with pairwise inner products at most (1/2). It adds 4 over the
[preceding 652,072-point certificate](https://github.com/alexlegeartis/KissingNumbers/tree/31a99695ad9ee55d4f63f9e0f478fdd9c11edc17/verifications/improved/dim38-qr-octad-switches/).

| Layer | Points |
| --- | ---: |
| Dense QR fibres (unchanged) | 253,952 |
| 3,088 full signed octads | 395,264 |
| Retained signed roots | 2,572 |
| Signed stars | 288 |
| **Total** | **652,076** |

The [witness](witness.json) specifies ten old/new octad replacements
at zero-based indices of the
[original support file](../dim37-38-qr-merge/supports_3088.txt). Each
eight-coordinate support has all 128 even-parity signings. The 18 star
templates have a center of magnitude 2 and four leaves of magnitude 1;
each uses both signs of the center and all eight odd-parity leaf
signings. Every two-coordinate root of magnitudes 2 is retained unless
its raw product with a star exceeds 4. Exactly 240 roots are removed.
The octads, stars and roots have squared norm 8.

The octad substitutions leave the octad count unchanged but create new
reciprocal incidences between star centers and leaves. Opposite
orientations of the same center pair can eliminate the same signed
roots, letting more stars share the deletion cost. In this witness the
stars exceed the removed roots by 48, versus 44 in the previous
certificate. The two new octads cannot be inserted separately: each
conflicts with both removed octads. Their joint replacement lets 16
additional stars share four already-deleted roots, for a net gain of four.

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
