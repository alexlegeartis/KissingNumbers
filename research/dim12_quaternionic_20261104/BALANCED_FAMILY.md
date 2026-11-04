# A three-parameter deformation of the quaternionic 840-point arrangement

This is a derivation made during the 4 November 2026 campaign. It gives a family
of **840-point** kissing arrangements, not a new kissing-number lower bound.
The notation refers to Sections III–IV of
[Takhanov, arXiv:2609.09179](https://arxiv.org/abs/2609.09179).

Keep the three axial 24-cells and the first mixed family M1 unchanged. For each
of the three 128-point pieces M2,r, replace its quaternionic component weights by

\[
 (1/2,\sqrt{a_r},\sqrt{3/4-a_r}),\qquad 1/4\leq a_r\leq1/2.
\]

The ternary relation and all component signs remain unchanged. Every point has
norm one. The canonical construction is the corner a0=a1=a2=1/2.

**Claim.** All pairwise inner products stay at most 1/2 throughout this cube of
parameters.

For the axial constraints, a C direction has weight 1/2 and a D direction has
weight at most 1/sqrt(2). The maximum absolute inner product of a C direction
with a D direction is 1/sqrt(2). Thus every axial-to-mixed inner product is at
most 1/2.

Within one M2,r piece, different unsigned label triples differ in at least two
positions, by the xor parity relation. At most one component contributes, and
its squared weight is at most 1/2. If the unsigned labels agree but signs differ,
at least one component is flipped. Each squared component weight is at least
1/4, so the resulting inner product is at most 1-2(1/4)=1/2.

Between different M2 pieces, the component frames have coherence 1/2. The
inner product is at most one half the scalar product of their weight vectors,
which is at most 1/2 by Cauchy–Schwarz. The unchanged M1-to-M1 constraints are
already valid.

For M1,r versus M2,s, the affine parity choices have a crucial property: the
three quaternionic component scalar products cannot all be nonzero. The exact
integer check in `exact_support_check()` enumerates all 9*16*16=2304 unsigned
label pairs and verifies this property. A nonzero directional scalar product
has absolute value 1/sqrt(2). The maximum contributions of the three blocks are

\[
 1/4,\qquad \sqrt{a_s}/(2\sqrt2),\qquad
 \sqrt{3/4-a_s}/(2\sqrt2).
\]

Each is at most 1/4, and at most two occur. Their sum is at most 1/2. This
completes the proof for every real parameter triple in the stated cube.

Writing T=a0+a1+a2, the frame operator is

\[
 \operatorname{diag}\bigl(78I_4,(30+32T)I_4,(102-32T)I_4\bigr).
\]

At a0=a1=a2=3/8 it becomes diag(78 I4,66 I4,66 I4), compared with
diag(78 I4,78 I4,54 I4) at the canonical corner. These two arrangements are
not isometric, since the frame spectra differ. At the balanced point the
coordinates lie in Q(sqrt(3)): the two D components of each M2 vector are
integer multiples of sqrt(3)/4, while M1 and the axes are rational.

This suggests replacing the two last axial 24-cells by two interacting
25-point regions. The search first optimizes 50 unit vectors and the three
parameters against the remaining 792-point core. It then releases the core,
including starts formed by superposing the published 841 displacement field
with its copy under exchange of the last two quaternionic factors.

The completed 512-start experiment did not produce a valid 842-point code.
Neither that failure nor the separate whole-cell obstruction proves a ceiling
for this continuous family with arbitrary additional points or deformations.
