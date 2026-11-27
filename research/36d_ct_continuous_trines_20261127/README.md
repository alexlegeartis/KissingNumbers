# 36D: exact real radial rigidity of the Coxeter–Todd two-block layers

Lower79, 27 November 2026. Current 36D bound **506,628**, from [live project witness](https://github.com/alexlegeartis/KissingNumbers/tree/a304d9f04756c67f07c7df9fc846b989c1bb0826/verifications/improved/dim36-ct-neutral-reheat). Named outside comparison **484,760** is [Yang et al., arXiv:2609.35051](https://arxiv.org/abs/2609.35051), also identified as unchanged by [Lin, arXiv:2611.03682](https://arxiv.org/abs/2611.03682). No claimed lower increase, no Git by worker, no native process. The present result concerns continuous directions that were outside earlier finite triple-root and mask exchanges.

The physical code is in \(K_{12}^3\) with raw squared norm 12 and raw cutoff 6. Let \(G\) be the explicit Gram in the [pinned CT skeleton](https://github.com/alexlegeartis/KissingNumbers/blob/a304d9f04756c67f07c7df9fc846b989c1bb0826/verifications/improved/dim36-ct-nonlocal/certificate.json). Write \(M=I-R\), where \(R\) is the Eisenstein \(120^\circ\) rotation. The quotient label map of \(K_{12}/MK_{12}\) has 729 values. The 224 norm-six classes have 18 short representatives each and split as six orthogonal regular trines. All three pairs of physical blocks keep a complete one-to-one matching of these norm-six labels, with the entire \(18\times18\) Cartesian product for every pair; their 217,728 old rows and all 2,268 axial rows are unchanged by the latest record delta. A trine in a norm-six class consists of \(r,Rr,R^2r\) of raw norm six and mutual products \(-3\). This structure is verified from coordinates and was also classical for \(K_{12}\); see [Conway–Sloane](https://doi.org/10.1017/S0305004100060746) or [Elkies' lecture notes](https://people.math.harvard.edu/~elkies/M272.19/index.html).

## Fixed norm-six head: single-plane trine is exactly rigid

Fix any norm-six head \(v\) in block \(h\) and any one of the original trine planes in a matched partner block \(p\). The original trine in that plane is given by \(r_i\). Matching preserves negation. Thus the actual old \((-Rv,-r_i)\) is present in the same physical pair and is outside the modified head fibre; \(v\cdot(-Rv)=3\). For any modified residual \(t\) in the chosen plane, \(-r_i\cdot t\le3\) for \(i=1,2,3\). The three walls cut out exactly the original closed triangle: its vertices are the \(r_i\). Its sole points of squared norm six are the vertices, by strict convexity. This excludes arbitrary real reorientation and arbitrary insertion within any original plane, not just the obvious squares, even if all other rows in the modified fibre are removed. Other physical block projections and mixing planes are not covered by this argument.

## The larger theorem: any genuine radial shrink fails if only two blocks are occupied

Fix any norm-eight direction \(v\) in any physical block \(h\) and any different block \(e\). Preserve only the 6+6 Cartesian layer and the axial layer identified above. A proposed real row \((\lambda v,y,0)\) on these two blocks (up to placement), with \(q=\|y\|^2=12-8\lambda^2\), is impossible for **every \(3/4\le\lambda<1\)**. At \(\lambda=1\) there is exactly one possible oriented root trine in the other block, before checking other old layers. No claim that all three rows of that trine automatically meet the rest of the record. This result does not exclude candidates that also leak into the third block or change the direction of \(v\).

### Why fourth moments arise from the old contacts

For every \(v\) and pair \(h,e\), among the 224 norm-six body classes exactly 32 have maximum raw body product 5 against \(v\), and their negatives all have max 4. The opposite class maps to the opposite matched class. Collect all 576 rays \(w\) in the matched 32 *tail* classes. Their trines give a zero-mean random scalar \(s=w\cdot y\), with second moment \(\mu_2=q/2\), and the old points with head products 5 and 4 give the two-sided bound
\[
   -B\le s\le A,\qquad B=6-4\lambda,\quad A=6-5\lambda .
\]
Both sides use actual old rows, and the maximizers of the body and tail contributions can be chosen independently because the rows are Cartesian products.

There are 192 complex lines in the selected tail stock. With \(P_i\) their complex rank-one projectors, the exact frame
\[
   \mathcal F(H)=\sum_i (P_i-I/6)\,\mathrm{tr}((P_i-I/6)H)
\]
has traces \(\mathrm{tr}\,\mathcal F=160\) and \(\mathrm{tr}\,\mathcal F^2=2,488,320/54^2=160^2/30\). The replay finds six mutually orthogonal norm-four root complex lines unbiased to every selected line, so each line's squared projection into any of the six is \(1/6\). Thus \(\mathcal F\) is positive semidefinite and supported on the 30-dimensional real space of off-diagonal \(6\times6\) Hermitian matrices in that basis. Equality in the trace Cauchy inequality proves \(\mathcal F=(16/3)I\) there. This avoids numeric eigenvalues and any group transitivity assumption.

For completeness, the exact projector product of two CT norm-six complex lines, with representatives \(r,u\) and \(a=r\cdot u\), \(b=Rr\cdot u\), is \((a^2+b^2+(a+b)^2)/54\). This follows from trine projection: \(\sum_i (r_i\cdot x)^2=9\|\operatorname{proj}x\|^2\). The replay sums \((a^2+b^2+(a+b)^2-9)^2\) over the selected ordered pairs via an exact bank-by-bank table. Six lines are tested unbiased by the same formula with a norm-four representative, whose condition is \(a^2+b^2+(a+b)^2=6\). Across all 20,412 body directions, there are 6,804 distinct 32-label patterns. Each is checked against all six ordered physical matches (40,824 exact cases).

Write \(t_j=\|\operatorname{proj}_{L_j}y\|^2\) for these six root lines. For projections on one tail trine, the three scalar values sum to zero and their fourth powers sum to half their second-power sum squared. Applying the frame identity to the Hermitian matrix \(yy^*\) gives
\[
  \mu_4 = \frac38\left(2q^2-\sum_{j=1}^6t_j^2\right).
\]
Explicitly the sum of squared projection energies over the 192 selected lines is \(192(q/6)^2+(16/3)(q^2-\sum t_j^2)\); multiplying by \(81/(2\cdot576)\) proves the stated formula.

The actual axial layer contains the six vectors \(\pm M R^k r\) for any norm-four line. They are a regular hexagon of raw length \(\sqrt{12}\) in that physical plane. Their support is at least \(3\|\operatorname{proj}_{L_j}y\|\), so each \(t_j\le4\). For \(3/4\le\lambda\le1\) we have \(4\le q\le7.5\), and convexity with \(\sum_jt_j=q\) gives \(\sum t_j^2\le16+(q-4)^2\). Therefore 
\[
  \mu_4 \ge \frac38(q^2+8q-32).
\]
On the other hand \((s+B)(s-\mu_2/B)^2\ge0\) gives \(\mu_3\ge-B\mu_2+\mu_2^2/B\). The nonnegative polynomial \(s^2(s+B)(A-s)\) gives \(\mu_4\le (A-B)\mu_3+AB\mu_2\); because \(A-B=-\lambda\), this becomes
\[
  \mu_4 \le B^2 q/2 - \lambda q^2/(4B).
\]
After multiplying the difference of the lower and upper bounds by \(8B\), elementary expansion gives
\[
 96(1-\lambda)(28\lambda^4-80\lambda^3+4\lambda^2+124\lambda-69).
\]
The quartic is positive at \(3/4\) (value \(87/64\)) and increasing until 1: its derivative at 1 is 4, and its second derivative \(336\lambda^2-480\lambda+8\) is negative on this interval. Thus the moment bounds contradict each other for every \(3/4\le\lambda<1\).

At \(\lambda=1\) both moment bounds equal six, forcing all tail energy into just one of the six unbiased root lines. In that plane the axial hexagon, at raw radius two, allows only its polar vertices, the six signed norm-four roots. To see their orientation, \(Mr=r-Rr\) and \((Mr)\cdot r=6\); the six roots are precisely the bisector directions. Check each of the two trine orientations per one of the six root lines against the old selected bank \(s\le1\). The exact replay finds precisely one orientation in all 40,824 cases (other constraints can only reduce the possible set). It suffices to check a single signed representative of each rotation orbit, because the selected old tail bank is \(R\)-invariant.

### Consequence for a tempting constructive strategy

A possible 3-to-4 trade takes the old heavy trine with head norm squared 8, shrinks that head by \(\lambda=\sqrt{15}/4\) to squared norm \(15/2\), and hopes to replace its three norm-four tails by a regular tetrahedron of tail squared norm \(9/2\). The easy construction lifts the three original tails \(r_i\) to \((r_i,u)\) and uses a fourth \((0,-3u)\), where \(\|u\|^2=1/2\). The fourth row occupies only two physical blocks, so it is impossible, independently of \(u\)'s direction and even before touching the many triple-root rows. To beat the old full trine one must deform the direction of the head, make all new rows genuinely three-block, or edit the old 6+6/axial layers together. Merely changing constant \(\lambda\) slightly within the proved range is not a different route.

### Replay

Save the companion note as `verify.mjs` and run Node 20+: `node verify.mjs certificate.json`, with the pinned skeleton above, SHA-256 `c4598a0dfc49365ae2b082745454c38b53b76703e5eef5dafbd71e2d06c59180`. The exact rational LDL uses integer scale 16 and exact norm revalidation for every lattice leaf. The finder verifies all trines and partner matches, unbiased root axes, trace squares, the unique endpoint trine, and all contact patterns. No file from the old optimizer, floating eigenvalues, or angular tolerance is needed. Companion: `verify.mjs`.


### Outward deformation is impossible even with both other blocks occupied

For any norm-eight \(v\) consider all 81 shortest vectors in its own CT quotient class. They have sum zero: the order-three rotation \(R\) preserves the label and \(I+R+R^2=0\). Therefore the 80 squared distances from \(v\) to the other members average \(81\cdot16/80=16.2<18\). Each is the squared norm of a nonzero vector in \(MK_{12}\), an integer multiple of six and at least 12. Some \(w\) in that same class therefore has \(\|v-w\|^2=12\), whence \(v\cdot(v-w)=6\). The old axial vector \(v-w\) occupies the physical head block alone and forbids \(\lambda>1\) for any prospective new row with that head direction, independently of both off-block tails. This adds no finite search premise beyond the already enumerated classes and the verified similarity.

Separately, an [exact exclusion on the 1,134 orphan directions](ORPHAN_AXES.md) checks the actual live witness and all signed two-axis blends; use [its independent replay](verify_orphan_axes.mjs).
