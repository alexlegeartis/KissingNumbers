# 19D: exact radius barriers for masked and multilevel Golay heads

**No new lower bound.** The comparator remains **12,460 by Kevin Li**, [Improved Kissing Numbers in Nineteen and Twenty-One Dimensions, §4](https://hexagonmath.org/pdf/2610.00248v1). His 21D 31,560 also remains unchanged. Project main at admission was d53024aff3f114773d87db90a8b60740ef63876e. This is a different low-support and multilevel analysis, not a repeat of the near-flat shoulder-charge or 21D periodic-line work.

## Body and exact angular constraints

Use Li's physical $\mathbb R^{16}\oplus\mathbb R^3$, with marked tetrad $T=\{1,4,7,9\}$ in the source length-20 coordinates. Let $E$ be his dimension-12 punctured Golay code, $E_H$ its projected 16-head code with kernel $\{0,T\}$, and $Z$ the 256 projected signs having a lift with zero on $T$. The body contains 480 head roots, 12 tail roots, 3,840 odd-signed disjoint-octad points, 6,144 shoulders and 192 full-$T$ caps, all of squared norm 8, total 10,668. A candidate is scaled to squared norm 8 so every body product must be at most 4.

The 30 disjoint eight-coordinate masks have binary rank 5 and are exactly the affine hyperplanes on the head $AG(4,2)$; their orthogonal complement is $E_H$. For each of the six antipodal directions $v$ of shape $(0,\pm1,\pm1)$, the 16 six-site shoulder supports form a symmetric **2–(16,6,2)** design: every coordinate occurs 6 times and every pair twice. Each $+v$ and $-v$ occurs with the complete even head parity on each of those supports. The four full-$T$ cap head masks partition the head into tetrads; each has all odd head signs and all six axial tails $\pm2e_i$. Every pair of head tetrads is a disjoint octad. In each shoulder family, a support picks two sites from each of three head tetrads, with each omitted tetrad occurring four times. These finite claims are checked from Li's twelve source generators.

For a candidate of specified sign and individual magnitudes, the exact maximum over a complete old sign parity on a support is the sum of magnitudes if the signs agree in parity, and otherwise that sum minus twice its smallest magnitude. A zero makes either parity equally good.

## 1. Arbitrary signs, arbitrary masked-flat support and arbitrary real tail

Let a head have common magnitude $u>0$ on any mask $M$ of size $k$, zeros elsewhere, and arbitrary signs and real tail $t$. For $1\le k\le15$, **its squared raw norm is at most $64/9$ in the old body's polar feasible region**, even relaxing sign compatibility to allow each whole support independently to choose the most favorable candidate parity.

For each of the six shoulder families define $p_v$ as its maximum intersection with $M$, replacing a full intersection 6 by 4 (optimistic mismatch). Similarly $p_0$ is the max disjoint-octad intersection with 8 replaced by 6; $p_T$ is the cap intersection with 4 replaced by 2. A partial intersection gets no subtraction. If $h_k$ is 2 for $k=1$, otherwise 4, necessary constraints are
\[
0\le u,\quad h_k u\le4,\quad p_0u\le4,\quad
p_vu+|v\cdot t|\le4,\quad |v\cdot t|\le2,\quad
p_Tu+2|t_i|\le4.
\]
This gives 158 signatures on 65,535 nonempty masks, of which 65,534 have $k<16$. Put $r=|t|_{(1)}+|t|_{(2)}$ (the top two absolute tail components). Then $\|t\|^2\le r^2$, $r\le2$, and the inequalities imply $p_{\min}u+r\le4$, $p_Tu+r\le4$. This two-variable rational polygon proves $ku^2+r^2\le64/9$ except at nine signatures; for those nine enumerate independent four-wall vertices of the actual profile polytope by Cramer's rule, getting the same $ku^2+\|t\|^2\le64/9$ bound. A convex quadratic attains its maximum at some vertex. Dropped vertices with $u=0$ have tail squared norm at most 4.

For $k=16$, any sign outside $E_H$ matches the old odd parity on some disjoint octad, so $8u\le4$. Reaching squared norm 8 would force $u=1/2$ and an axial tail of length 2 (because $\|t\|^2\le r^2\le4$, with equality only on an axis). A facing full-$T$ cap then contributes 4 in the tail and optimistically at least $2u$ in the head, contradiction. Hence all full uniform heads lie in $E_H$.

A full uniform head from $Z$ is impossible: the shoulders give $6u+r\le4$. With raw squared norm 8 and $\beta=\|t\|^2\in[0,4]$ the left side is at least $(3/2)\sqrt{8-\beta}+\sqrt\beta>4$ by concavity and the endpoints $3\sqrt2,5$.

Each of the other 1,792 $E_H$ signs occurs **at most once**. Even the optimistic shoulder and cap inequalities give $r\le B=4-4u$ and $\|t\|_\infty\le A=2-u$. Supposing $\beta\ge2$ gives $u\in[1/2,\sqrt6/4]$ and $B/2\le A\le B$. In sorted absolute tail coordinates take the third equal to the second for an upper bound. Maximizing over the resulting polygon gives $\beta\le\max(3B^2/4,A^2+2(B-A)^2)$. Subtracting $8-16u^2$ from these expressions gives $4-24u+28u^2$ and $4-28u+35u^2$; both are negative at both endpoints, hence throughout by convexity. Thus $\beta<2$, body-head norm squared $>6$, so two points with the same signs have raw dot $>6-2=4$. Keeping the body and extending by the entire masked-flat class therefore gives **at most 12,460**. This is not a claim that the restricted class attains that many; Li's axis points are nonflat.

## 2. Zero-tetrad signs with arbitrary two levels or four separate tetrad levels

For any zero-tetrad sign, allow magnitudes $b+c$ on any high mask $M$, $b$ on its complement, with $b,c\ge0$, and any real tail. The exact squared polar radius is **$64/9$**, with no restriction on amplitude ratio. Let $k=|M|$ and now define $p_v$ *without subtraction* because a zero-tetrad sign aligns with every even shoulder. Use the same $p_0,p_T$ as above: these are the coefficients of $c$ after the required mismatch on the odd octads/caps. The constraints imply
\[
4b+h_k c\le4,\quad 6b+p_0c\le4,\quad
6b+p_{\min}c+r\le4,\quad 2b+p_Tc+r\le4,\quad 0\le r\le2.
\]
Across the 241 wall signatures, the three-variable bound $16b^2+2kbc+kc^2+r^2\le64/9$ succeeds except at nine. For each remaining one introduce $a=\|t\|_\infty$, with $0\le a\le r\le2a$ and $2b+p_Tc+2a\le4$, and use the sharper upper bound $a^2+2(r-a)^2$ on the tail. The four-variable vertices prove the same radius. Uniform $u=2/3,t=0$ really attains it.

Separately allow independent head magnitudes $A\ge B\ge C\ge D\ge0$ on the four whole tetrads. The pair-union octads force $2A+B\le2$; the shoulder census forces $2(A+B+C)+r\le4$; caps force $A+a\le2$. With the same tail cone, enumerate five-wall vertices and bound $4(A^2+B^2+2C^2)+a^2+2(r-a)^2$ by **$64/9$** again, attained by the uniform template. Thus taking different amplitudes on the four whole tetrads is also not enough. Neither result covers arbitrary within-tetrad amplitudes.

Any of these radius-gap statements persists for a **3% same-norm Euclidean head perturbation** at fixed tail and sign: if a raw norm-8 point were feasible, its same-norm head template at distance at most $0.03\|H\|$ would have body products at most $4+0.24$. Scaling it by $1/1.06$ would give a feasible template of squared norm $8/(1.06)^2>64/9$.

## 3. Arbitrary amplitudes on at most five head sites

This geometry does not need a common head amplitude, a Golay sign, or numerical tests. The head and tail roots enforce top-two absolute sums at most 2 in each block. Up to four head sites thus give head squared norm at most 4 (if top value $x$ exceeds 1, use $x^2+3(2-x)^2\le4$). Total squared norm 8 would then require an axial tail of length 2. The facing full-$T$ cap on any active head site contributes at least that head magnitude plus 4, contradiction. Zero head also cannot reach 8.

Five head sites give head squared norm at most 5, with equality only five ones, so the tail squared norm is at least 3. If the largest absolute tail coordinate $a\le1$, all three absolute tail coordinates are 1, and the five head amplitudes must all be 1. In a shoulder family facing tail product 2 there is a support meeting the fixed five-set in at least three positions: by the 2-(16,6,2) identities the ratio of the sum of squared intersection sizes to their sum is $(5+2)/3>2$. At least three head signs can align since there is an unused sixth site, producing product at least 5. If $a>1$, the tail inequality $a^2+2(2-a)^2\ge3$ forces $a\ge5/3$. Facing axial caps give each individual head amplitude at most $4-2a\le2/3$, and five such amplitudes plus tail squared norm at most 4 cannot reach 8.

## Replay and scope

Save the companion JavaScript as verify.mjs and run node verify.mjs (no input files, no packages). It reconstructs the exact source octads and checks all ranks, incidence and rational vertices in about five seconds hosted. The obstructions only use the unchanged base, so no dense witness is needed. They are not a general upper bound on the kissing number nor a result for an altered base.