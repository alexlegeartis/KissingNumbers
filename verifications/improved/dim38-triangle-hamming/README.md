# 38D: 611,798 by parity-split triangular Hamming latitudes

The construction keeps the equator, axes, and Fano grid of the
[611,638-point witness](../dim38-joint-latitudes/). It replaces the
512-point orthogonal Hamming layer with 768 points and makes 96 net
additional cap deletions, gaining **160**. For comparison, Qihong Lin's
[external preprint](https://arxiv.org/abs/2611.03682) gives 605,322 in 38D.
The older comparisons elsewhere in the repository are historical.

## Construction and angle check

Let b be a signed norm-32 minimal Leech vector and let y be a seven-sign
tail of squared norm 7 in one of the two orthogonal seven-coordinate
tail halves. The inserted unit vectors are
`(b/sqrt(96), sqrt(2/21) y)`.

The original assigned an antipodal pair of body vectors to each tail.
Here every tail receives three oriented minimal bodies forming an
A2 triangle. Hamming parity chooses between two disjoint five-coordinate
body blocks, occupying Leech coordinates 14–18 and 19–23 (zero-based).
In either block four explicit triangles and their negatives supply eight
mutually disjoint oriented packets. The binary Hamming syndrome of the
tail selects a packet. The 16 tail labels at each syndrome have mutual
Hamming distance at least three.

For one tail, the three bodies have mutual normalized product −1/2,
so the pairs are exactly at the kissing threshold. Tails at Hamming
distance one have opposite parity, hence orthogonal bodies. At distance
two the syndromes differ and no oriented body is repeated, so the body
product is at most 1/2. For distance at least three even equal bodies
fit. Tails in different halves are orthogonal and need no further
restriction. The same triangle palette works for both halves.

There are 24 distinct new cap-body lines, disjoint from the 14
Fano-grid body lines: 38 cap lines removed in total versus 22 before.
Of the eight former Hamming owners, five are restored, and 21 other
lines are newly removed. The final counts are:

| Component | Points |
|---|---:|
| Retained caps | 589,452 |
| Equator | 13,710 |
| Axes | 1,596 |
| Fano grid | 6,272 |
| Triangular Hamming | 768 |
| **Total** | **611,798** |

Both point generators check all 294,528 new unordered pairs with
integer dot products. Independent cross-layer checks include the exact
cap and old-grid geometry of restored cap lines. The original cap-cap
proof uses the inherited Golay/Leech minimum-norm argument, not a new
enumeration of all unchanged pairs. The existing cap/axis and cap/equator
bounds apply to new cap subsets and restored caps. The manuscript and
previous proof package are unchanged.

## Replay

From the repository root, with Python and NumPy:

```sh
python verifications/improved/dim38-triangle-hamming/verify.py
python verifications/improved/dim38-triangle-hamming/audit.py
```

The first verifier checks pinned data hashes and rebuilds the triangles
from signed coordinate cycles. The second independently enumerates
Leech minimal lines and constructs the 768 new points using explicit
root palettes and a parity-check matrix.
