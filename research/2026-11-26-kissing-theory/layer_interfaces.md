# Exact interface rules for layered kissing constructions

This note isolates the transferable part of the [14D Hamming–Witt construction](https://github.com/alexlegeartis/SpherePackings/tree/main/constructions/dim14): design a compatibility rule from the actual geometry, then choose code families and shifts to realize it. A probability, average distance, or rounded Gram matrix is not a substitute.

## Affine block compatibility with unequal amplitudes and different codes

Let the common horizontal coordinates be indexed by J. A layer has the form
\[
B_i = \{(a_i\odot(-1)^\sigma, b_i):\sigma\in s_i+C_i\},\qquad
\|a_i\|^2+\|b_i\|^2=1,
\]
with nonnegative amplitude profile \(a_i\) and an affine binary code. Zero/missing amplitude coordinates can be deleted from its individual sign word. For two layers project their linear codes onto their common nonzero coordinates and put \(D_{ij}=\pi_J C_i+\pi_J C_j\) in \(\mathbb F_2^J\). The differences of the two sign words run through *exactly* the coset \(\pi_J(s_i+s_j)+D_{ij}\); this does not require independence of individual bits. For each coordinate \(\ell\), put \(w_\ell=a_{i\ell}a_{j\ell}\). Their worst possible dot product is exactly
\[
b_i\cdot b_j+\sum_{\ell\in J} w_\ell-
2\min_{e\in\pi_J(s_i+s_j)+D_{ij}}\sum_{\ell:e_\ell=1}w_\ell .
\]
Thus with a parity check H for D_ij, the syndrome Hπ_J(s_i+s_j) must not be reached by any error whose weighted cost is strictly below \((b_i\cdot b_j+\sum w_\ell-1/2)/2\). Equalities are allowed. Within a block check its nonzero differences separately; if J is empty only the tag dot remains. This recovers the usual unweighted coset-leader condition but also handles different linear codes, different supports and weighted profiles. Related weighted-distance instances already occur in the [62D archive](../../verifications/improved/dim62-coordinate-weighted-distance/).

**Cheap structural gate before searching shifts.** Every syndrome of H has a representative on any basis of the columns of H. Consequently the largest possible coset-leader cost is at most the total weight of a minimum-weight column basis (compute it greedily by matroid independence). If even that number is too small for the required cost, **no** relative shifts can make the two families coexist. In the unweighted case, if the relative tag geometry demands ordinary distance at least k, the two projected codebooks must share at least k independent parity checks: \(\dim D_{ij}^\perp=\dim((\pi_JC_i)^\perp\cap(\pi_JC_j)^\perp)\ge k\). This is only necessary; the cheap-error syndrome table gives the exact test.

## A concrete consequence: at the rank-three interface, two 128-word d4 layers must share a tetrad partition

Take two binary linear [12,7,4] codebooks before shifting. Suppose they share exactly three independent parity checks, and some choices of cosets have cross-distance at least three. Let H be their common rank-three check matrix, and s the syndrome of the difference of the shifts. There is no error of weight 0, 1, or 2 with syndrome s, so s≠0 and s is neither a column of H nor a sum of two distinct columns. The four pairs in \(\mathbb F_2^3\) under addition of s include {0,s}; from each other pair at most one distinct nonzero column can occur. By rank three there are exactly three distinct nonzero columns, forming a basis. The label s is their sum.

Complete H by two more check rows for either of the original codes. The twelve full columns must be distinct and nonzero and no three distinct columns sum to zero, by minimum distance four. Within any fixed common-prefix group, the remaining suffixes are distinct two-bit words. If a zero-prefix column exists, its suffix is a nonzero t. In each nonzero-prefix group, no two suffixes differ by t, or those two columns and the zero-prefix column would sum to zero. Each such group then has size at most two, and the zero-prefix group has at most three. This gives at most nine columns, contradiction. So the twelve columns split into **three common tetrads of four**, each with a nonzero, independent prefix, and the two private check bits run bijectively through their four values within each tetrad. This is a real restriction on choosing different large layer banks, not a new record or a claim that the case of four or five common checks is excluded.

## 26D: an exact orbit interface for the 21-row gate

The [24-donor joint construction](four_octahedra_21.md) pays three skeleton rows and needs one 21-site clipped chamber to improve Takhanov's 201,212. The chamber cutoff is invariant under the involution \(R(u,r,v)=(u,-r,v)\). Each pair {x,Rx} is valid within the chamber: its dot equals \(2(u^2+v^2)-1\le0\). Between two such pairs, the worst dot is exactly \(uu'+vv'+|r\cdot r'|\). Thus whole two-site families can be chosen and checked as single objects.

For any admissible chamber code, hold all complete mirror pairs fixed. On the k singleton representatives put an undirected edge ij iff \(x_i\cdot R x_j>1/2\). This is an exact conflict because R is an orthogonal involution; a reflected singleton can never conflict with the symmetric core or its own mate. Within the union of the original k points and their k reflected copies, the conflict graph is bipartite: edges go only across the two orientations. Let ν be the maximum matching size in this doubled graph. By König's theorem the exact optimum keeps 2k−ν of those rows. Equivalently doubling an independent set I of singleton labels and removing all its neighbors gains |I|−|N(I)|, keeping every other original label; the best such gain is k−ν. Indeed the ordinary bipartite Hall deficiency maximizes |A|−|N(A)| over subsets A. Replace A by I=A\N(A), an independent set with N(I) contained in N(A)\A and no worse deficiency.

In an odd chamber of size 21, k is odd. A single broken pair can always be doubled, producing 22. If k=3, any graph other than a triangle has positive gain (double an isolated label, or double the two endpoints of the path and delete its middle vertex). Thus a 21-row chamber that cannot be enlarged even within its mirror union either has at least five broken pairs or exactly three whose actual mirror conflicts form a triangle. The latter has a geometric restriction: with h_i=(u_i,v_i) and r_i the three regular coordinates, every edge implies \(r_i\cdot r_j<h_i\cdot h_j-1/2\). Hence the three squared regular norms and their pair products imply
\[
0\le\|\sum r_i\|^2 < \|\sum h_i\|^2-2\sum\|h_i\|^2.
\]
This is impossible if the triangle contains both opposed mixed quadrants (+,−) and (−,+). In each hot coordinate the two opposite signs then force its squared sum to be at most twice the sum of the three coordinate squares (cancel the opposite part and apply Cauchy–Schwarz to at most two remaining terms). A mirror triangle is genuinely possible elsewhere: take u=v=−1/2 and regular triples (0,1,−1)/2, (1,−1,0)/2 and (−1,0,1)/2.

This gives a structured comparison of a few complete layer orbits and genuine symmetry-breaking rows. It does not make the 21-point chamber finite in amplitudes and does not yet supply a 21-row chamber or any improved numerical bound.
