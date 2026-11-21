# BCH syndrome gluing doubles the flat layers in dimensions 94–96

| Dim. | Pinned incumbent | Added net points | Combined lower bound |
|---:|---:|---:|---:|
| 94 | 74,134,063,820 | 2,097,152 | 74,136,160,972 |
| 95 | 117,083,737,156 | 2,097,152 | 117,085,834,308 |
| 96 | 209,278,995,264 | 4,194,304 | 209,283,189,568 |

Retain the dense, shifted Golay–dodecacode, sparse Reed–Solomon and
root layers of the linked predecessor configurations. **Replace**, do not
append to, their previous flat Reed–Solomon layer. The old settings
retained in 94D and 96D are the [flat Reed–Solomon configuration](../dim92-96-flat-rs-layers/).
The enlarged [108-head E7 95D configuration](../dim95-e7-108-heads/) is retained.
The result is a lower bound and not an optimality claim.

As before, divide the ordinary four-coordinate blocks into six groups:
sizes (4,4,4,4,4,4) in 96D and (2,4,4,4,4,4) in 94–95D. All new
heads in 94–95D are zero. A support word is a six-tuple of block
indices with binary XOR zero. It selects one block per group. Inside
the six selected blocks place signs of magnitude one. The new points
have squared norm 24.

To choose their 24 signs, use the length-24 shortened extended binary
BCH code in GF(32) defined by the 11-bit columns
(H_x=(1,x,x^3)), for field elements indexed 0 through 23. Use
the polynomial (x^5+x^2+1) to realize the field. Partition these
columns in their displayed order into groups of four, independent
of which physical block a support word prescribes. The full check
matrix has rank 11, and every subset of one to five columns is
independent in the sense that its sum is nonzero. Therefore every
kernel coset has exactly 8,192 sign patterns with pairwise Hamming
distance at least six. Within the same physical support their dot
product is at most (24-2cdot6=12).

For a block index (ain{0,1,2,3}) in group (i), define the map
(phi_i) by binary-linearly combining the corresponding two phase
vectors with entries (0,4,913,1759,417,772) and
(0,1025,1987,327,1632,1281) respectively. These entries encode elements
of the 11-bit syndrome space, and combinations use XOR. The signs
allowed on support (s) form the full affine coset of syndrome
(igoplus_iphi_i(s_i)).

Distinct support words differ in at least two groups. If they differ
in at least three, they share at most 12 coordinates and clear the
dot threshold regardless of signs. If they differ in exactly two
groups (i,j), the XOR difference of their block indices is the
same nonzero value (d) in both missing groups. Their sign syndromes
differ by (phi_i(d)oplusphi_j(d)). In the quotient by the columns
of the two erased groups, this difference is neither zero nor any
single surviving column. Thus the signs differ on at least two of
the 16 shared coordinates, also giving dot at most 12. The exact
checker tests all 45 quotient exclusions directly.

New rows are the same flat rows as in the preceding geometry: every
active row has cross dot at most 4 against any pinned dense body
alphabet or shifted Golay–dodecacode body row. The sparse RS points
have six raw entries of magnitude four, and roots have two entries
of magnitude sqrt(48), so cross bounds are also at most 24. Every
inherited point has squared norm 96, making 24 the required cross
cap. New and inherited directions are different by their active-row
and nonzero-coordinate profiles, even after normalization. Different
flat supports are distinct, and their signs are uniquely decoded.

Run `python3 verify.py` (Python 3.10+, standard library only). It
replays local field arithmetic, the kernel distance, all quotient
exclusions, row compatibility and counts.
