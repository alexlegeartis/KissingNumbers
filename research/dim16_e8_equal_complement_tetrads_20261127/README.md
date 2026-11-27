# 16D E8⊕E8 equal/complement tetrad obstruction (sharp)

**No new lower bound.** The comparison is **4,320**, attained by the classical [Barnes–Wall 1959 construction](https://doi.org/10.1017/S1446788700025064), also recorded with attribution by [Cohn](https://cohn.mit.edu/kissing-numbers/) and discussed by [Takhanov, arXiv:2611.00927](https://arxiv.org/abs/2611.00927). The project main at admission, [e8d76b4](https://github.com/alexlegeartis/KissingNumbers/tree/e8d76b400beb62502b83408cee22675d366b1b24), did not claim a new 16D bound. This is a sharp obstruction for a specified coordinate family, not a general 16D upper bound and not the same as a restriction on all mod-two pairings.

## Statement and geometry

Use the conventional E8 lattice
\[
L=\{x\in\mathbb Z^8 : \textstyle\sum x_i\equiv0\pmod 2\}
 \cup \{x\in(\mathbb Z+\tfrac12)^8 : \textstyle\sum x_i\equiv0\pmod 2\}.
\]
Keep the 480 axial points formed by the norm-two E8 roots in either block, normalized to unit length. Among new candidates, permit
\[
  (v,w)/\sqrt8,\qquad \|v\|^2=\|w\|^2=4,
\]
in the following three distinct classes: (A) \(v,w\in\{\pm2e_i\}_{i=1}^8\); (B) \(v,w\) are signed unit tetrads with **equal** supports in the common coordinate indexing; (C) signed unit tetrads with **complementary** supports. Permit arbitrary signs, parity twists, nonantipodal selections, and partial or whole coordinate families in these three classes. Then the total is at most **4,320**, and this is attained.

We compare all new pairs in their unnormalized coordinates; their raw dot threshold is four. A norm-four E8 vector has dot at most two with every norm-two root: integrality and Cauchy–Schwarz bound it by two since \(\sqrt8<3\). This is exactly the new-to-axial threshold. Thus nothing conceptually restricts the freedom of signs within the three families except the new pair angles.

* **A: at most 256.** There are only \(16^2\) signed-axis pairs.
* **B: at most 1,792.** Attach to each new point its four signed-three-position restrictions, recording for each selected coordinate **both** block signs. For any fixed triple and its six signs, the common raw dot is six. Two points attached to this token but with different fourth coordinate would conflict; if the fourth coordinate agrees, they coexist only if both remaining signs are opposite, giving at most two completions. There are \(\binom83 2^6\) tokens, each with capacity two, and each point uses four. Thus the cap is \(2\binom83 2^6/4=1792\).
* **C: at most 1,792.** Encode each coordinate by a letter of \(\mathbb Z_4\), with parity specifying which block is occupied and the other bit specifying its sign (use \(a+2s\) if block \(a\) and sign bit \(s\)). The norm-eight dot of two rows is \(8-d_{\rm Lee}\), because changing side lowers the local product by one and changing sign on the same side lowers it by two. Thus the Lee distance is at least four, and every word has four odd letters. The Gray isometry \(0,1,2,3\mapsto 00,01,11,10\) gives 16-bit distance-four words with exactly four unequal bit pairs. Puncture one coordinate. The remaining radius-one Hamming balls are disjoint, and each center lies in the region with three or four unequal pairs among the seven complete pairs. This region has \(2^8(\binom73+\binom74)=17,920\) words. Each ball meets it in exactly ten: the center and its flip in the unmatched bit, plus the eight flips in the four pairs that can change their parity without leaving the allowed 3/4 counts. Therefore the cap is \(17,920/10=1,792\).

These bounds add to \(480+256+1792+1792=4320\), even before requiring all cross-family constraints.

**Attainment.** Index the coordinates of each block by \(\mathbb F_2^3\). For each of the 14 affine planes \(A\), include both (balanced) octad supports \((A,A)\) and \((A,A^c)\), with all 128 even-total signings. Take all 256 signed-axis pairs in addition to the 480 axial roots. Different plane supports in \(\mathbb F_2^3\) are equal, complementary, or meet in two. Consequently the 28 distinct 16-coordinate octad supports have intersections at most four, including the case of the same left plane with different right choices. On a fixed support, even-total signings differ in at least two signs. The signed-axis pairs have raw dot at most four with a paired tetrad (at most two per block) and with each other (one equal signed axis can contribute only four; the other is then zero or negative). This proves the literal stock is a kissing code of size 4,320 and attains all three class counts. We make no new historical or isometry claim about this explicit realization of the classical count.

## Exact closure of the real block domain, but not of the support choices

Any real nonaxial unit point compatible with the **two complete fixed E8 axial shells** necessarily has the form \((v,w)/\sqrt8\) above, where \(v,w\) are arbitrary norm-four E8 vectors. This is a useful reduction, but does not require them to be among the equal/complement-support classes.

For one block set \(u=\sqrt2 x_{\rm block}\). Write its sorted absolute coordinates \(a_1\ge\cdots\ge a_8\ge0\), and let \(\epsilon\in\{0,1\}\) be the parity of its signs (irrelevant at a zero). The roots \(\pm e_i\pm e_j\) and the even half-sign roots impose \(a_1+a_2\le1\) and \(\sum a_i-2\epsilon a_8\le2\). In differences \(\delta_i=a_i-a_{i+1}\), \(a_9=0\), this is an orthant with two additional linear inequalities. A vertex has at most two positive differences, and if there are two, both inequalities are equalities. The [companion integer checker](verify_e8.mjs) enumerates the two-axis vertices from precisely these equations. The maximum squared norm is one; its magnitude profiles, multiplied by four, are \((4,0^7),(2^4,0^4)\), and \((3,1^7)\) only at odd sign parity. Strict convexity shows any equality case is among these vertices. They are exactly \(v/2\) for E8 norm-four \(v\). Both block norms \(\|u\|\le1\) must be tight for a unit point of the two-block product.

In particular novel opportunities must involve a different pairing of support partitions, mixed axis/tetrad support, half-integral phase letters, or a changed old axial stock. A separate same-residue theorem also controls some half-integral directions, but is not a premise of the equal/complement argument.
