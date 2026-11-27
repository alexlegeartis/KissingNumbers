# 13D: arbitrary-real cap bounds for PackingStar's E7 × E6 fibre product

Lower80, 27 November 2026. Named comparison **1,154**, [Zinoviev–Ericson (1999), Theorem 1](https://www.mathnet.ru/eng/ppi457); strict target 1155. No new bound. The mixed product is **already published** by Ma et al., [PackingStar arXiv:2511.13391v5, Appendix F.3–F.4](https://arxiv.org/html/2511.13391v5#A6.SS3). Its connection to [Ganzhinov's 1932 in 14D](https://arxiv.org/abs/2207.08266) motivated this independent reproduction, not a new attribution. The new result is that keeping the entire 1,008-point mixed core with its conventional E7 and E6 projections makes **any** replacement of its two pure blocks, by arbitrary real kissing codes, total at most **1,151**. This applies to any gluing edge set of that size which uses all these root directions in both projections; it is not an unrestricted 13D upper bound.

## Independent certification of the published 1,146

Take \(E_7\) in \(\mathbb R^7\) with raw norm squared 4: the 14 signed axes of amplitude 2 and the sixteen signed unit entries on each of the seven complements of the cyclic Fano lines \(\{i,i+1,i+3\}\). An \(E_6\) is its 72-root slice perpendicular to \(t=(1,1,0,1,0,0,0)\). Partition the first set into nine 14-point orthoplexes and the second into nine 8-point orthoplexes; the checker constructs both. Corresponding frame products give 1008 mixed points of squared norm 8. With neither coordinate repeated, products are bounded by \(2+2\); if one repeats, the other block is nonpositive. The 84 pure-first points are signed \(D_7\) pair roots of amplitude 2. The 54 pure-second points are the \(E_6^*\) minuscule directions: for every \(r\in E_7\) with \(r\cdot t\ne0\), take \((3r-(r\cdot t)t)/\sqrt3\), of raw norm squared 8. They pair with E6 roots by at most \(2\sqrt3\). The checker tests all 1,146 literal rows and every distinct product by exact integer norm-ratio comparisons. No new isomorphism class asserted.

## Pure E6: at most 54 without a lattice or equal-amplitude restriction

In unit length, the pure E6 directions belong to the ordinary E6 Voronoi cell: against each raw root \(s\) of squared norm 4 the condition is \(y\cdot(s/\sqrt2)\le1\). The traditional E6 roots \(s/\sqrt2\), of squared norm 2, are its relevant vectors. The Voronoi polytope has 54 vertices of squared norm \(4/3\), the two minuscule classes. See [Pook-Kolb–Allen–Agrell, *Glued lattices are better quantizers than K12*](https://arxiv.org/abs/2312.00481) and [the root Voronoi characterization](https://doi.org/10.1112/jlms.12456).

For unit \(y\) in the convex hull of the 54 vertices \(w\), take its scalar product with its convex decomposition: some \(w\) has \(y\cdot w\ge1\). Thus the normalized vertices cover the spherical slice by closed 30-degree caps. A cap containing a strict interior site cannot contain any other kissing site. Charge each such site to one distinct cap. Any remaining site is on at least two boundary caps: if exactly one vertex attains \(y\cdot w=1\), the convex decomposition uses that vertex alone, contradicting the norms. Each remaining cap contains at most two boundary sites: equality of their angular distance at 60 degrees forces mutually opposite transverse directions. Count these incidences. There are at most 54 sites, and the old pure E6* block attains 54.

## Pure E7: at most 89 by a weighted real cap inequality

For any unit kissing code and unit pole \(u\), give a site with \(\alpha=u\cdot x>1/\sqrt2\) weight \(2\alpha^2-1\), and all others zero. The total at each pole is at most one. Indeed let \(a_i=\sqrt{2\alpha_i^2-1}\), \(z_i=x_i-\alpha_i u\). Then
\[
 z_i\cdot z_j\le \tfrac12-\alpha_i\alpha_j\le -a_i a_j/2\quad(i\ne j),\qquad
 \|z_i\|^2=(1-a_i^2)/2 .
\]
The squared difference between \(2\alpha_i\alpha_j-1\) and \(a_i a_j\), in the sense of the difference of their squares, is \(2(\alpha_i-\alpha_j)^2\). Therefore \(0\le\|\sum_i a_i z_i\|^2\le(S-S^2)/2\) for \(S=\sum a_i^2\), proving \(S\le1\).

Use as poles the 56 unit vectors supported with all signs on a Fano line. For a raw norm-squared-eight pure E7 direction let \(x\) be its absolute coordinate magnitudes, so each \(x_i\le2\) and the sum on each Fano complement \(Q\) is at most 4. These are precisely the full signed E7-root tests. On a Fano line, sort its three magnitudes \(a\ge b\ge c\). Only the signed sums \(a+b+c\) and \(a+b-c\) can exceed \(\sqrt{12}\); any other signed sum is at most 2. The total pole load is
\[
 \mathcal L(x)=\sum_{\text{Fano lines}}
 \left[\left((a+b+c)^2/12-1\right)_+ + \left((a+b-c)^2/12-1\right)_+\right].
\]
The exact certificate proves \(\mathcal L(x)>5/8\) on the permitted norm-eight sphere. The verified 168 permutations preserving Fano lines reduce the 5,040 coordinate orders to 30. Each order makes the fourteen sums literal linear forms \(f_i\). If \(\mathcal L(x)\le c_0=5/8\), then for any subset \(I\), \(\sum_{i\in I} f_i(x)^2\le12(|I|+c_0)\). Consequently, for each integer vector of weights \(w_i\), Cauchy gives
\[
   \sum_{i\in I}w_i f_i(x)\le \sqrt{12(|I|+c_0)\sum_{i\in I}w_i^2}.
\]
Only **42** such tangents across the 30 orders are needed; the rounded-outward rational right sides are checked by squaring. Intersect these cuts and the seven quartet walls with the ordered cube. The latter is the simplex of points with its first \(k\) coordinates equal to 2 and all remaining zero. For each added wall the checker enumerates every vertex exactly by rational edge clipping. Two current vertices form an edge iff their common active inequalities define a face containing exactly those two vertices. Clipping retains old vertices inside and creates new vertices only on crossing edges. The largest squared norm among all final vertices of all 30 outer polytopes is
\[
 \frac{1027520673009327277}{128457728100000000}=7.998901181\ldots<8 .
\]
This proves the strict load bound. Summing the at-most-one capacity over all 56 poles gives \(N(5/8)<56\), hence \(N\le89\). Consequently cap-only rebuilding cannot reach 1155: \(1008+89+54=1151\). No claim the pure E7 ceiling 89 is attained.

## The entire natural 812-point catalogue has sharp capacity 84

The raw norm-eight integer catalogue comprises 84 D7 points and 672 Fano spikes: a signed amplitude 2 at a heavy site and four signed units on two of the three pairs from the Fano lines through it. Also include all 56 deep Fano rays. The checker verifies every integer donor is the scaled midpoint of an edge joining two deep directions of unit product \(1/3\): sum the two signed triples of raw norm 3 and dot 1. The finite lists coincide, all 756 distinct. A midpoint pair sharing one endpoint \(u\) conflicts if its two other endpoints have product \(+1/3\), because the normalized midpoint product is \(3/4\). Thus all their other endpoints have mutual product \(-1/3\) (they cannot be antipodal with equal positive projections on \(u\)). Their normalized orthogonal projections away from \(u\) have mutual dot \(-1/2\); their number is at most three by taking the squared norm of their sum. Every chosen deep ray disallows all incident midpoint rays. If \(d\) deep and \(m\) midpoint rays are used, \(2m\le3(56-d)\), giving \(m+d\le84\), attained by D7. This is exact for arbitrary signed selections of the catalogue, not a numerical optimization bound or a global spherical-code bound.

## Consequence / exact replay

A new result in this product must genuinely change the mixed geometry. Deleting an E6 root direction loses all fourteen mixed partners whenever a proposed pure point breaches its old wall. Independently deforming only the two cap shapes is insufficient; unequally weighted fibres with simultaneous changes within the full orthoplex, a projected root system other than the old E6/E7 pair, or a coupled equator/latitude trade in the external 1154 stock remain distinct mechanisms.

Save the companion note `verify.mjs` as `verify.mjs` and run `node verify.mjs` (Node 20+, no external files). It generates the two root systems, verifies both frame partitions, all 1,146 products, all 756 midpoint identities, and the 42 exact tangent walls with integer arithmetic. Hosted replay passed. No native/shared job or Git.

A separate [weighted-cap theorem for the actual 1,154-point Zinoviev–Ericson code](ZE_CAP.md) includes [an independent replay](verify_ze_cap.mjs).
