# 16D: exact real scaffold rigidity and signed cheap-trade obstruction

**No new lower bound.** The fixed comparison is the classical **4,320** of [Barnes–Wall](https://doi.org/10.1017/S1446788700025064) and the Latin families in [Rustem Takhanov, arXiv:2611.00927](https://arxiv.org/html/2611.00927v1), §§2.1–2.3. These are scoped obstructions and an exact failure of one more general weighting strategy, not a 16D upper bound. The separate E8⊕E8 equal/complement result is not used here.

## Complete 1,248-point quaternionic scaffold forces an exact alphabet

In each of four orthogonal \(\mathbb R^4\) blocks, take the unit 24-cell \(C=\{\pm e_i\}\cup\{(\pm1,\pm1,\pm1,\pm1)/2\}\). Its three levels \(C_0,C_1,C_2\) are respectively the signed axes, even half signs, and odd half signs. The alternate cell \(D=\{(\pm e_i\pm e_j)/\sqrt2:i<j\}\) has three orthoplex levels, one for each perfect pairing of the four coordinates. Keep the four axial copies of \(C\) and, for every block pair \(i,j\), all \((d_i,d_j)/\sqrt2\) with the three \(D\)-levels coupled by any fixed permutation on that pair. These are the 96 plus 1,152 scaffold points of the cited construction.

**Theorem.** The complete set of real unit rays compatible with this fixed scaffold is *exactly* \((c_0,c_1,c_2,c_3)/2\), with each \(c_i\in C\). No sign symmetry, equal weights, or alphabet is assumed in the input.

For a block of the new ray, sort its absolute coordinate magnitudes \(a\ge b\ge c\ge d\ge0\). The axial inequalities give \(a\le1/2\) and \(a+b+c+d\le1\). The largest support sums for the three matchings, before dividing by \(\sqrt2\), are \(a+b,a+c,\max(a+d,b+c)\). Put \(p=a+b\), \(r=\max(a+d,b+c)\). Write \(a=b+e,c=d+f\). If \(r=a+d\), then \(e\ge f\) and
\[
pr-(a^2+b^2+c^2+d^2)=(b+d)e+2d(b-c)-f^2\ge0.
\]
If \(r=b+c\), then \(f\ge e\) and the difference is \(a(c-e)+c(b-c)-d^2\ge0\). Thus the squared block norm is at most \(pr\).

For distinct blocks \(i,j\), the two-block walls imply \(p_i+r_j\le1\), regardless of their level permutation: the largest level sum in block \(i\) is paired to some sum at least \(r_j\). Take a block with maximum \(p=P\) and its \(r=R\). Every other block has \(p_j\le\min(P,1-R)\) and \(r_j\le1-P\), so the total squared norm is at most
\[
PR+3\min(P,1-R)(1-P).
\]
For \(R\le1-P\) this is at most \(4P(1-P)\), with equality only \(P=1/2\). For \(R\ge1-P\) it is \(3(1-P)+(4P-3)R\); if \(P\le3/4\) maximize at the lower endpoint, and otherwise bound \(R\le P\) to obtain \(3-6P+4P^2\le1\). That latter inequality is strict for \(3/4<P<1\). The endpoint \(P=1\) forces \(a=b=1/2,c=d=0\) by the axial constraints, so \(R=1/2\) and the bound is strict. Consequently for a unit ray \(P=1/2\), and all four blocks must have squared norm \(1/4\), equal to their \(p^2\). The elementary inequality \(a^2+b^2+c^2+d^2\le(a+b)^2\) is tight only for one nonzero coordinate, or all four equal, since \(c^2+d^2\le2b^2\le2ab\). These are precisely \(C/2\). Conversely each such vector meets all axial and matched-level walls.

This theorem is independent of which Latin cubes and which D-level permutations are chosen. It reduces all scaffold-preserving real replacements to the same finite signed \(C^4\) domain.

## Literal counts against every full uniform-level Latin incumbent

Take three uniform-level layers \(M_t\), each made of one arbitrary 64-word Latin cube on the four axes of each level, with all 16 independent signs per word. These add 3,072 points to the scaffold for the 4,320 comparison. An added signed ray has cost equal to the number of distinct incumbent rays of \(M_0\cup M_1\cup M_2\) with which it strictly conflicts.

For a given level \(t\), let \(k\) be the number of that level among the four coordinates of the proposed ray. Scale by the common factor, so compatibility means the sum of four inner products of unit C letters is at most two. The exact blocker counts in \(M_t\) are:
* \(k=0\): zero, since all four contributions are at most \(1/2\).
* \(k=1\): 16, one sign choice for each of the 16 old words matching that ordinary axis, with all other contributions \(+1/2\).
* \(k=2\): four, one sign choice for each word matching both axes, with both other contributions \(+1/2\).
* \(k=3\): 20. The unique old word matching all three axes blocks in two sign choices. For exactly two matching axes, each of three axis pairs has three other words, and two sign choices on the mismatching same-level axis, giving 18 more.
* \(k=4\): if the proposed word is absent, eight. Each choice of three axes has one distinct old word which matches those three and not the remaining axis, with two signs free in the remaining block.

The Latin projection properties give these counts exactly, independently of signs and of the detailed cubes. Distinct new candidates therefore have **8** blockers (uniform absent word, or 2+2 level counts) or **36** (3+1 or 2+1+1). These are per **individual signed point**, not per sign orbit.

## A pure counting proof forbids all arbitrary cost-eight trades

Fix any old anchor in \(M_t\). In each block rotate the anchor to \(e_0\), relabeling the three levels by the same cyclic shift. This is possible by the coordinatewise left multiplications in the binary tetrahedral group from Takhanov, because the anchor levels are all \(t\) and the quotient by the signed-axis group is cyclic. The 216 cost-eight neighbors of this anchor are of two types:

* 24 one-slot switches: keep three letters equal to the anchor, and put another signed ordinary axis in the changed slot.
* 192 two-slot mixed neighbors: keep two letters equal to the anchor. On the other two slots, both use the same non-anchor \(C\)-level, and the first half-sign coefficient is \(+1/2\) so its dot with the anchor letter is \(+1/2\). There are six choices of changed pair, two levels and four oriented letters in each.

We show that at most eight neighbors of these types can be mutually compatible. Write their non-anchor half signs as \((1,\sigma_1,\sigma_2,\sigma_3)/2\). In an active slot, their three-bit patterns are distinct: for two mixed rows with the same changed pair their distances on those slots must sum to at least four, and with different pairs sharing one slot the distance there must be at least two. Each active word has at most one opposite-parity partner at distance three (its bitwise complement), and at most one at distance one: a distance-one partner must use the same changed pair, and the word at the other end is then forced to be the unique complement there. Two such partners would repeat an active word, impossible. Therefore at a fixed slot there are at most four active mixed rows, either all in one parity or at most two per parity.

At most two one-slot switches can use the same slot, and if two they are opposite versions of the same alternative axis. A single switch forces one prescribed free bit on every mixed row active at that slot: its local product must be \(-1/2\), since the other three slots contribute \(2+1/2\). In that 2-dimensional face of the bit cube, distance-three partners do not exist; the preceding argument gives at most two active mixed rows. If two switches occupy the slot, no mixed row can use it, since it would need two opposite bit values. If \(m_i\) is the active mixed-row count and \(u_i\) the one-slot switch count, this proves \(m_i+2u_i\le4\). Summing over four slots counts each new candidate twice and gives at most eight. The eight one-slot switches formed by one antipodal pair per slot attain this local cap.

Each new cost-eight point requires eight old deletions. Each deleted old point can be charged by at most eight mutually compatible new ones, by the cap just proved. Thus **any simultaneous trade using only cost-eight rays cannot gain**, with the scaffold retained, even with arbitrary individual signed deletions and no orbit invariance.

## Two further exact escape requirements

**One unchanged Latin layer.** If a uniform Latin layer is kept (even if fewer than four of its points are removed), no truly new ray can use its level because the minimum blocker count in that layer is four. New rays avoiding it lie in the union of the other two \(C\)-levels in every block. After rotating levels, these are 16 ordinary binary signs; the required Hamming distance is four. Puncture to 15 bits, apply the binary radius-one sphere packing bound, and obtain at most 2,048 such points. Thus a gain under this scaffold needs at least four deletions from **each** old uniform layer.

**Saturated two-level complement.** Conversely, any 2,048-point code using only the other two levels blocks every mixed ray with exactly \(k=1,2,3\) occurrences of the absent level, even if that binary code is nonlinear and unrelated to any Latin cubes. Scale the coordinates so the candidate has one entry of magnitude two and three zero entries in each of the \(k\) axis blocks, and four entries of magnitude one in the other blocks. Its total absolute sum is \(16-2k\); dot threshold is eight. Puncture one zero coordinate; the optimal binary layer punctures to a perfect radius-one code on 15 bits. Consider the binary fill subcube matching all candidate nonzero signs, with \(3k-1\) free bits. If the candidate had no strict blocker, the perfect neighbor of each fill must make exactly one error at an occupied position whose weight drops the dot to at most eight. There are respectively zero, two, or seven possible occupied error locations for \(k=1,2,3\). Two fills with distance one or two cannot use the same error location (otherwise their perfect neighbors differ only in those fills, at forbidden distance below three). But the center and its single-bit flips give \(3k\) pairwise distance-at-most-two fills. This exceeds the number of available occupied error locations, a contradiction.

## Reciprocal weighted Hall does not close expensive trades

Charging all candidates by the reciprocal of their individual blocker count is **not** a valid universal anchored inequality. Against the same anchor there are eight cost-eight rays and one cost-36 ray that are mutually compatible. The [companion replay](verify_anchor.mjs) gives their nine literal signed rows, checks the 36 mutual products and their strict products with the anchor, and classifies costs by the level counts above. This makes the proposed sum \(8/8+1/36=37/36>1\). It is a counterexample to this charging strategy, not a profitable global trade or an upper bound on expensive trades.
