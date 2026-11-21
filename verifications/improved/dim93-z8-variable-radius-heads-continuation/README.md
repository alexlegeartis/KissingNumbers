# 44 variable-radius Z8 fibers in dimension 93

**47,389,968,384 kissing points**, up **1,073,741,824** from the preceding variable-radius 93D witness, with the hexacode layer retained. This is a lower bound, not an optimality claim.

## Variable-radius dense fibers

Use the [pinned Z8 parent](../dim88-89-z8-prefix-anticode/), with its exact length-24, rank-12 generator and homogeneous minimum 24. Drop the first two four-coordinate blocks as in 88D. Every body has squared norm 88, and each chosen two-symbol prefix leaves (8^{10}) words. The selected prefixes and new heads in `heads.json` assign 44 distinct pairs (color, prefix). The three body colors are the same two Walsh and one axial crosspolytope alphabet as in the parent. Every new head is an integer five-vector of squared norm at most 8.

Normalize **individually** by (sqrt{88+|h|^2}). If two distinct words share a fiber, their body dot is at most 40 by subtracting the deleted prefix dot 8 from the inherited cap 48; thus (40+|h|^2le(88+|h|^2)/2). In different same-color fibers, subtract the actual deleted prefix dot from 48 and add the new head dot. In different colors, use 22·2=44 for the body before adding the head dot. For every one of the 946 fiber pairs these signed raw upper bounds are at most half the product of their radii; the exact checker tests by signed squaring with integers. Thus no duplicates appear, and the dense count is (30cdot8^{10}=47,244,640,256).

## Different sparse layer

Take the [88/89D sparse variable-puncture layer](../dim88-89-golay-variable-punctures/) and place it entirely in the unchanged 22 body blocks, with zero on all five head coordinates. It consists of 982 weight-11 pattern fibers with 1,024 axis and 144 sign choices and 11 weight-10 fibers with 256 axis and 72 sign choices. Its own pair proof is unchanged. Each occupied row has dot at most 4 against every dense body row. For (w=10,11), the squared sparse radius is (8w), and (4(4w)^2le8w(88+|h|^2)) guarantees normalized compatibility even for the smallest dense radii. Dense bodies have all 22 rows occupied, so sparse points do not duplicate them. This adds precisely 145,004,544 points.

Run `python3 verify.py` in the installed sibling layout, or give the common parent directory as a positional argument (Node.js is needed by the inherited sparse replay). It checks all integer head comparisons, SHA-bound Z8 source, rank and prefix counts, sparse replay and all new cross-layer bounds. The pinned Z8 minimum theorem remains an inherited proof premise.

The published table includes the additional flat hexacode layer below. The raw head entries in heads.json are integer numerators with denominator 1.

## Flat hexacode appendix

An additional 323,584 points use 79 six-row supports in flat_supports.json, pairwise meeting in at most three rows, none contained in the support of any retained sparse pattern. On each support take two complementary parity patterns (three even and three odd flat rows). The four projective axes in each parity class are the sign vectors in four coordinates modulo overall sign, with raw entries of magnitude one. Choose the axis sequence with the length-six quaternary hexacode and choose an even-parity overall sign on the six axes. The 64 axis words and 32 sign words are checked directly.

The squared radius is 24. New points from different supports have dot at most 12; opposite parity assignments on the same support also have dot at most 12. Within a parity assignment distinct axis words share at most two axes, while different signs on the same axes disagree at least twice. Dense rows contribute at most 18 across the six slots by the balanced parity assignment; retained sparse patterns meet at most five rows, giving dot at most 20. These satisfy the normalized cap for all old radii. Occupied-row counts distinguish all three structures.

Run python3 verify_flat.py (Node.js for the inherited sparse proof) for the current published count in claims.json. It replays verify.py and the count in claims_base.json before checking all flat supports, code distances and layer geometry.
