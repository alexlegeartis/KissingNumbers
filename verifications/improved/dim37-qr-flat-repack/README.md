# Dimension 37: 544,666 by paired-coordinate parity refill

Run `python verify.py`. The standard-library-only replay verifies the
544,474-point predecessor and all the new interactions exactly.

Use the ordinary coordinate flat U = 12..23 (zero-indexed). Its three original
octads carry 384 points. Every remaining octad and radical six-support
meets U in at most four coordinates. Any norm²32 ordinary vector on U with
eight coordinates of magnitude two is therefore compatible with the exterior
octads and radical supports, the roots, the marked poles, and the dense layer:
raw dot ≤16 with each norm²32 point, zero with marked poles, and at most
32 with dense points of squared norm at least 145.

## The new local code

Split U into six consecutive pairs and split the pairs into triples
A=(0,1,2) and B=(3,4,5). Use each of the nine 8-supports obtained by
omitting one pair from A and one from B. Each support contains 64 phases:
- take independent simultaneous sign flips on each occupied pair (16);
- take an affine copy of the projection of the binary [6,2,4] code whose
  two repeated length-three blocks are equal and have even parity (4).

A parity word constrains the product of the two signs in each coordinate
pair. Distinct phases on one support are at sign distance at least two:
pair flips change both signs, and the restricted parity code has minimum
distance two. Supports sharing only two pairs are automatically compatible.
For supports sharing three pairs, require that their restricted affine
parity codes be the two distinct cosets of the common projected code. This
is a linear condition on affine offsets and ensures a sign mismatch in a
shared coordinate, so their signed overlap is at most four.

`win.txt` lists the sign group (basis on line one), then support masks,
initial negative masks and orbit sizes. `search.py` optionally derives
suitable offsets using sparse GF(2) elimination. The complete integer
verifier also directly checks all 165,600 new pairs. The resulting count
is 544,474 − 384 + 576 = **544,666**.

## Alternate Golay construction

`python golay_alternate.py` independently generates 480 projected Golay
points and 96 extra points from three shortened orbits of the involution's
fixed subcode. It checks all 165,600 pairs. These are a separate 576-point
refill on the same flat, not extra points added to the first refill.
