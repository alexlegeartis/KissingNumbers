# Forty-six Z8 fibres from bipartite root heads in dimension 93

**49,545,659,928 kissing directions**, an increase of **2,147,485,448**
over 47,398,174,480 (about **4.53%**). This is an explicit lower bound.
It replaces the entire preceding 44-fibre dense stock with 46 fibres
and keeps all 153,534,224 certified body-only rays in the current 93D
construction. It also adds 1,800 D93 roots incident to the head (40 pure and 1,760 mixed).

The body is the exact [preceding 88-coordinate Z8 body](../dim93-z8-variable-radius-heads-continuation).
A two-symbol prefix fixes two deleted blocks and leaves 8^10 codewords.
Each retained row has squared norm 4, and its alphabet is one of the
same two Walsh frames or the axial frame. The body has norm squared
88. Each new head lies in five orthogonal coordinates.

Take the 40 heads 2(±e_i±e_j) of the D5 root arrangement, together
with six signed unit heads. Give the 24 roots crossing the split
{0,1}|{2,3,4} one body colour. Their two prefix symbols record the
signed axes i and j separately. This equates their head Gram with
the original two-block Gram, so their compatibility follows directly.
The remaining 16 roots are on same-side axis pairs. Assign them and
the six unit heads to the two other colours; the explicit choices
are in the witness. This is a joint recolouring/replacement of the
previous head design, not an append to its head rays.

Normalize each fibre by its own radius sqrt(88+||h||²). For distinct
words of the same fibre the body dot is at most 40 by the pinned
homogeneous minimum 24 of the source Z8 code. The head squares
are 8 or 1, so adding the head preserves the kissing inequality.
For different prefixes of the same colour subtract the exact
two-block prefix dot from the source cap 48. For different colours
all body rows have mixed dot at most 2, and hence body dot at most
44. Add the literal head dot and check the normalized constraint
by integer squaring for all 1,035 pairs. Strict separation rules out
duplicates. The source minimum theorem is inherited, not inferred
from the rank test.

All retained points have zero on the five head coordinates and
their previous dense-interface proofs use only the same three
body frames and the body squared radius floor 88. This includes
the Golay weight-ten and weight-eleven families, the retained
seven-row flat modes, all matching modes including the six-row
Reed–Solomon repack, and both triality appends and ordinary roots.
The verifier pins those geometry descriptions and their complete
count chain, and checks their independent conservative cross bounds.
Distinct body occupancy separates the sparse rays from all new
dense directions. The previous sparse–sparse comparisons are unchanged. The additional
D93 roots at ±e_i±e_j use at least one head axis, have squared norm 2,
and distinct pair dots at most 1. Against a dense point their raw
dot is at most 4, since every dense coordinate has magnitude at most 2
and the dense squared radius is at least 88. Against a retained sparse
point a root tests at most one ordinary coordinate. The verifier
checks the corresponding mixed inequality for each retained family,
including equality for the ordinary D4 roots. No old sparse ray is
duplicated because it is head-zero, and no dense ray is duplicated
because every dense body occupies all 22 rows.

Run python3 verify.py with assertions enabled from this directory. It
uses only the standard library and the pinned sibling packages. Each
witness row contains the colour, the packed prefix 8*a+b, and the
five integer head coordinates. No optimizer is used.
