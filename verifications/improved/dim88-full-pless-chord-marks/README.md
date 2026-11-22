# Full Pless census and new chord marks in dimension 88

The explicit lower bound is **13,358,363,692**, up **569,772** from
13,357,793,920. The dense and small-layer counts are unchanged;
the sparse layer grows from 959,965,380 to 960,535,152 rays.
This is a replacement of the previous 88-dimensional realization,
not an append to it. Its other layers use the same recipes in the
re-labelled coordinates.

Keep fixed ternary sites 0 and 3. Change the four chord sites
from {1,2,24,28} to {9,30,39,44}. Projection onto these six coordinates
has full rank six in the same QR48 code. Hence the dense count,
norm and distance proof of the pinned 88D construction still apply.
All old small layers transfer by assigning the 42 ordinary planes in
the same order to the new 42 ordinary slots. Their comparisons
depend only on the number and arrangement of ordinary planes and
the per-plane alphabets, not on the old coordinate labels.

Instead of selecting a Pless automorphism orbit, enumerate all
weight-15 words in the pinned systematic [48,24,15] ternary parent
with zeros at sites 0 and 3. A weight-15 word has exactly one
information half of weight ≤7. If the generator is [I|A], the
self-orthogonality identity AA^T=−I gives A^(-1)=−A^T. Thus either
half determines the other, and the two enumeration branches never
overlap. The exact counts by axial occupancy are
37,392 / 80,576 / 57,696 / 17,088 / 1,552. Apply the same phase
counts 7,812 / 6,561 / 2,187 / 729 / 729 as in the predecessor.
The parent distances and the geometry of the three basic shifted
letters are unchanged; the original per-word and mixed distance
arguments therefore apply. The old marks are also enumerated as
a check, reproducing their earlier orbit count exactly.

Run `python3 verify.py` with a C++17 compiler available. The source
certificates and predecessor proofs are hash-pinned; their ternary
parent minimum distance is inherited. No large projected point list
or optimizer is needed for replay.
