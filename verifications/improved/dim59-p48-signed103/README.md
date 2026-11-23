# 59D: eight new tails in place of six

The explicit lower bound is **56,969,166**, two above the
[inherited construction](../dim59-p48-signed101/). Remove signed pure
tail indices 19, 33, 79, 82, 88 and 89 (zero-based); add the eight
directions in witness.json at pure amplitude 1. All equatorial and
mixed points remain unchanged. The two-point change is relative to the
previous repository result, not the cumulative gain over
[Kamenetsky's P48p construction (2026)](https://arxiv.org/abs/2611.01906).

Run `python verify.py` with Python 3 and assertions enabled (standard
library only). The checker pins the old literal tails and exact
tetrahedral and weighted interfaces. Each inserted direction is
normalized and compared with the retained poles and the other donors.
The exact interface inequalities are the same as in the pinned
parent, whose verification covers the stock left unchanged.
