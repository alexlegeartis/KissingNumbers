# Delsarte–Goethals donors on (D_4) tail roots

| Dimension | Preceding lower bound | New lower bound | Net gain |
|---:|---:|---:|---:|
| 66 | 434,139,558 | **434,387,366** | 247,808 |
| 67 | 453,138,696 | **454,068,488** | 929,792 |
| 68 | 483,338,010 | **485,390,106** | 2,052,096 |

The [preceding three-site construction](../dim65-68-dg-diagonal-triples/)
uses full 64-coordinate signed bodies from 80 cubic translates of
RM(2,6), plus weighted variants and sparse sign banks. The only
predecessor deletions here are the explicitly counted sparse banks
and hyperplane-head fibres; every old dense body and weighted variant
is preserved.

Select new cubic cosets at Hamming distance at least 14 from all 80
old cosets. In each chosen coset take a quadratic translate of the
same 131,072 binary Delsarte–Goethals words, with minimum distance 24.
Use signs of magnitude one for the 64 coordinates and the tail
((31/10)r), where (r) is a root with two signed unit coordinates.
The 66D, 67D and 68D tails use respectively 2, 3 and 4 coordinates.
Each row of selection.json specifies a two-coordinate tag, tail
sign parity, cubic label, and quadratic shifts for both orientations.
Quadratic shifts assigned to the same cubic label produce disjoint
DG stocks. Distinct body points in such a label are separated by
at least 16. For each other pair, compute the actual cubic-coset
distance. The same-tail, mixed-tail, full old-body and old weighted
body inequalities are all exact. The latter crucially use cancellation
between opposite signs on the first two tail axes.

A donor needs a sign mismatch only against a norm-16 sparse support
with exactly its own pair of tail coordinates. Since the 14 body
coordinates are orthogonal to RM(2,6), compare total parity on the
support, including the two signed tail coordinates. An opposite
parity forces a mismatch of donor magnitude at least one. Delete
only the bank if this guarantee is missing. For an atlas fibre,
use the exact minimum body Hamming distance after restriction to its
affine hyperplane and the literal tail dot. Count each dropped fibre
once, both across orientations and across new labels. Each bank or
fibre costs 2,048 points. In 68D, the 22 tails add 2,883,584 points,
at a cost of 312 banks and 94 atlas fibres.

Run python3 -B verify.py with a C++20 compiler on BMI2-capable x86-64.
This replays the hash-pinned predecessor, recomputes the inherited
bank parities, then runs clearance.cpp to check all required
finite Hamming distances. Only integer and rational inequalities
are used. No numerical search is needed for replay; this proves
lower bounds, not global optimality or worldwide priority.
