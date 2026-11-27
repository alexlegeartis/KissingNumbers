# 34D/35D: whole QR puncture obstruction and exact 35D split-phase gate

**No lower bound changed.** At [main 944a7ed](https://github.com/alexlegeartis/KissingNumbers/tree/944a7edcab2771dc37271fba583c2b0b19ff8048), the verified bounds remain **391,004** and **415,360**, versus Qihong Lin's **390,980** and **415,308** [arXiv:2611.03682](https://arxiv.org/abs/2611.03682). The companion exact Node script replays both finite gates.

## Whole-binary-layer QR obstruction, arbitrary offsets

Let \(H\) be Yang et al.'s actual \([39,15,12]\) shortened binary QR kernel, with literal generators in [their Appendix B](https://arxiv.org/html/2609.35051v1#A2). For a marked set \(M\), puncture those coordinates and place all signs from a union of complete cosets of the punctured \(H\) in the ordinary Euclidean embedding. A 34- or 35-coordinate binary sign code needs Hamming distance at least 9. Five complete 32,768-word cosets would have been a material gain over the 131,072-point dense layer of either incumbent, with all full signed octads and roots safe (ordinary dot at most 16 and 8, respectively). This hoped-for mechanism is false:

* There are exactly **23** four-coordinate marks preserving minimum distance at least 9 (and **zero** five-coordinate marks). Enumerate the actual 970 weight-12 supports of \(H\). Since the next weight is 16, a four-mark set survives iff not contained in any of those supports; a five-mark set survives iff all its four-subsets survive. Also the punctured kernel retains rank 15 because the original minimum is 12.
* For **each** of the 23 marks the punctured length-35 kernel has covering radius **8**. Hence at most one full coset, even if the offsets are **arbitrary**, not only Yang's ten named cosets. This leaves other kernels, partial cosets and genuinely nonbinary/variable norms outside the result.

The covering calculation uses a single exact 24-bit syndrome quotient of the unpunctured \(H\). Map each unit-coordinate vector to its integer syndrome by Gaussian elimination. Distances \(D\) are the minimum number of original coordinate columns adding to the syndrome. Dynamic programming on an increasing collection of coordinate columns updates the XOR-paired cells in both directions, once per column. The full histogram for distances 0 through 11 is
\[
(1,39,741,9139,82251,575757,2836720,6944617,5468341,859047,554,9).
\]
For puncture marks \(M\), let \(L\) be the span of their four syndrome columns. The quotient distance is exactly \(\min_{\ell\in L}D(s+\ell)\): the removed coordinates are free and deleting them cannot increase the cost. A projected syndrome has distance at least nine exactly when all sixteen syndromes of its \(L\)-coset have master distance at least nine. The replay checks all members of that master tail and verifies zero such cosets for all 23 marks. It also finds at least one projected coset of distance eight in each.

## Different route: 35D two-coordinate skew and two spare axes

Use the actual [35D root/star incumbent](https://github.com/alexlegeartis/KissingNumbers/tree/944a7edcab2771dc37271fba583c2b0b19ff8048/verifications/improved/dim35-joint-star-root) with fixed 32 Cheng–Sloane dense coordinates and spare sites \(0,3,18\) (zero-based), and restore the actual signed-octad support list from its pinned predecessors. Choose two active spare coordinates \(p,q\) and leave \(r\) zero, and two distinct ordinary coordinates \(i,j\). For a parent dense word \(x\) of 32 signs, define
\[
v_i=b x_i,\quad v_j=-b x_j,\quad
v_p=\epsilon_p h,\quad v_q=\epsilon_q h,\quad h=\sqrt{16/3-b^2},\quad
b\in[1.05,1.06].
\]
Replace \(x\) by \(x+v,x-v\). Their angle is exactly \(60^\circ\), since \(x\cdot v=0\), \(\|v\|^2=32/3\), and both squared norms are \(128/3\). Every other old dense point is safe by \(16+2b<\sqrt{1024/3}\), for original signed difference distance at least eight. Every signed root is also safe: the largest pair of coordinate magnitudes is below \(2(2.06)\), giving a dot below \(8(2.06)<\sqrt{1024/3}\). Old star/head direct products can only be further constraints.

Represent each old octad by a full sign parity (may be changed independently, but keep all 128 rows and all its supports). For a proposed daughter, let \(L\) be the sum of absolute magnitudes on an octad and \(m\) the minimum. The worst old signed dot is \(2L\) for a matching parity and \(2(L-2m)\) otherwise. The angular threshold on \(L\) is \(\sqrt{256/3}\). The following **exact occupancy classification**, constant on the displayed interval, yields a small finite gate:

* Geometrically impossible if an octad contains \(p,q\) but not \(r\), and both \(i,j\); also impossible if an octad contains all spare coordinates but exactly one of \(i,j\). In the first case even the forced parity mismatch at the small coordinate does not suffice; in the second case there is a zero coordinate so parity cannot help.
* After ruling those out, the only mandatory parity changes are on octads containing \(p,q\), not \(r\) (zero or one of \(i,j\)), or octads containing exactly one active spare, not \(r\), and exactly one of \(i,j\). On the latter, only the daughter with the heavy ordinary coordinate triggers the condition. If both branches trigger on a common octad, their aligned parity is identical. All other octads are safe regardless of parity.

These inequalities can be checked without an optimizer: \(2.05<h<2.057\) and \(9.237<\sqrt{256/3}<9.238\) on the interval. For instance, with both active tails and both perturbed bodies and no unused spare, \(L-2m=6+2h>9.238\); with both active tails, exactly one perturbed body and no unused spare, the heavy branch must mismatch, and afterwards both bounds are \(4+2h+b<9.237\). If all three spares and both perturbed bodies are present, the maximum is \(3+2h+2b<9.237\); if exactly one is present then the heavy branch exceeds 9.238 with no parity help. The remaining cases follow by the same sum formula.

Exact finite conclusions on this family:

1. Only **19** of the 1,488 physical geometries pass the occupancy conditions. For each of those and each of four head sign patterns, the allowed old parity choices under all 262 actual signed stars and 12 actual full-support heads force a contradictory GF(2) system in the 17 parent coefficient bits. Thus no such parent can be split with the whole old sparse stock retained, even if octad sign parity is reassigned.
2. With the octad parities fixed to *any* shared assignment, the mandatory phase restrictions have rank **17** for each of the 19 geometries. Switching exactly one active tail sign cannot preserve all the required parities under any change of Cheng–Sloane body word. Switching both has only the all-ones body difference. Therefore any one such geometry permits at most one antipodal pair of donor parents for a fixed full octad parity assignment. This conclusion does not require the old stars/heads.
3. A separate exploratory owner-exchange reached only a tie, by splitting a single antipodal pair, retiring four old full-support heads and restoring two old pure-tail roots. This exploratory tie is not used as a premise and is not a certified positive gain. Partial phase deletion, additional shapes and mixed geometry are not ruled out.

Run `node research/dim34_35_qr_and_split_gate_20261127/verify.mjs .` from a copy of the repository (the four inputs are validated by hash).
