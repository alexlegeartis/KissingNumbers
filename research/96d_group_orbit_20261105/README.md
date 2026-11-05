# 96D group-orbit exchange and an exact matching obstruction

This is a **restricted search result**, not a new kissing bound and not a
no-go theorem for other transformations or partial-orbit constructions.

Let C be the free Z8 Golay code of size 2^36. Its subcode H=2C has size
2^24, and C is the disjoint union of 4,096 H-cosets. Choose a diagonal
coefficient vector α with entries 1 or 5 and a shift t in Z8^24, and put
D=αC+t. Multiplication by an odd unit and ring translation preserve the
homogeneous metric, so D is another distance-24 code. Since α≡1 modulo 4,
αH=H, so D is also a union of 4,096 H-cosets. On the binary Gray images,
the group H acts isometrically on both families.

Construct a bipartite graph between the C-cosets and D-cosets. An edge means
that some pair of words from the two cosets has binary distance below 24.
Any collection of whole cosets whose union is a code is an independent set in
this graph. A matching with 4,096 edges therefore bounds its size by 4,096
cosets: this particular C∪D whole-orbit exchange cannot add a single H-orbit
beyond C. This does not exclude breaking H-orbits, modifying other layers, or
using a different transformation.

`matching_certificate.json` records one explicit α and t and 64 exact
conflict witnesses. The α=5 positions form one binary Golay octad A. The
subgroup K={h∈G: A∩h∈G} of the binary Golay code has size 64. Translating
each witness by the corresponding common symmetry produces 64 disjoint
conflict edges; the 64 witness families cover every C- and D-coset exactly
once. `python research/96d_group_orbit_20261105/verify_matching.py` derives K,
reconstructs all 4,096 word pairs and checks their distances and both sets of
orbit labels exactly. Thus the result does not rely on the numerical decoder
or a claim that unobserved edges are absent.

The same matching rules out a larger set of equivalent shifts for this α.
Since α−1 equals 4 on A and zero elsewhere, C∩αC consists exactly of the
2C-cosets whose residue h lies in K. Consequently |C∩αC|=2^30 and
|C+αC|=2^42. Replacing t by t+c+d, with c∈C and d∈αC, gives an isomorphic
conflict graph: d leaves the affine code unchanged and translation by c
preserves C. The certificate therefore covers 2^42 shifts, representing 64
distinct affine αC-cosets. It also transports under the verified cyclic shift
of the first 23 coordinates; this is symmetry reduction, not a conclusion
about every diagonal transformation or every shift.
