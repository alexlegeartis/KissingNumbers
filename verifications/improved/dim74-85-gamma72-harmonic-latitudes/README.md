# Harmonic-theta recovery of Γ72 shell profiles

Twelve improved kissing lower bounds: dimensions 74–85.

| Dimension | Preceding | Verified bound | Gain |
|---:|---:|---:|---:|
| 74 | 6,218,211,606 | **6,218,213,549** | +1,943 |
| 75 | 6,218,276,173 | **6,218,331,159** | +54,986 |
| 76 | 6,218,400,136 | **6,218,635,448** | +235,312 |
| 77 | 6,218,669,870 | **6,219,013,264** | +343,394 |
| 78 | 6,218,992,206 | **6,219,716,626** | +724,420 |
| 79 | 6,219,416,084 | **6,220,590,652** | +1,174,568 |
| 80 | 6,219,962,106 | **6,221,713,348** | +1,751,242 |
| 81 | 6,220,774,236 | **6,222,197,896** | +1,423,660 |
| 82 | 6,222,089,254 | **6,224,427,690** | +2,338,436 |
| 83 | 6,223,768,006 | **6,227,374,552** | +3,606,546 |
| 84 | 6,225,863,462 | **6,231,140,058** | +5,276,596 |
| 85 | 6,228,418,680 | **6,235,824,404** | +7,405,724 |

The 76D case illustrates the construction. Keep all 6,218,175,600 norm-eight
equatorial vectors of Nebe's lattice Γ72 and the 24 pure tail D4 roots scaled to squared norm 8. For every oriented
norm-ten body v and every nonzero integer ray r in rays.txt propose
((4/5)v, sqrt(8/5) r/||r||). All proposed points also have squared norm 8;
their pairwise products with both retained codes are at most 4. For two
proposed points with body product j, the strict conflict threshold on
unit-tail dot products is (25−4j)/10. Since both heads have norm ten and
Γ72 has minimum8, only j=4,5,6 and 10 matter. Exact integer tests build
the conflict matrices G_j. A graph independent set is a valid latitude.

Fix an oriented adjacent pair of bodies v,w with v·w=j. For a third
body u put a=u·v and b=u·w. We simultaneously constrain the joint
(a,b) profiles of both the norm-eight and norm-ten shells. The replay
uses the moduli of the supported inner products, positive semidefinite
Gram tests and minimum-eight tests on short linear combinations. The
exceptional vectors ±v, ±w, and ±(v−w) in the appropriate shells are
accounted for individually. Antipodal symmetry identifies (a,b) with
(−a,−b); translations by small integer combinations of v,w identify
matching profiles in or across the two shells. Neither arbitrary
swap symmetry nor a transitivity hypothesis is assumed.

The new constraints come from the classical modularity of harmonic
weighted theta series. For each homogeneous harmonic polynomial H of
even degree k, its weighted theta series has weight 36+k. Since H(0)=0
and the lattice has minimum8, it is divisible by Δ^4; the quotient is
holomorphic of weight k−12. For k=12,14,16,18,20,22 the quotient space
has basis 1, none, E4, E6, E4², E4 E6, respectively. Comparing q^4 and
q^5 gives these factors between norm-ten and norm-eight harmonic sums:
−96, 0, 144, −600, 384, −360. For k=14 both sums are zero (the replay
only needs the norm-ten identity). Together with the ordinary low-order moments,
these constraints leave just 7, 4 and 3 free variables for j=4,5,6.

For the harmonic projection of m=(u·v)^r(u·w)^s of degree k=r+s,
write H=sum_l c_l||u||^(2l) Δ_Euclidean^l m, with c_0=1 and
c_l=−c_(l−1)/(2l(2k−2l+70)). The Laplacian on a,b is
10∂_a²+2j∂_a∂_b+10∂_b². The modular relations are applied to these genuine
harmonic polynomials, not to raw high-order moments.

For a fixed first tail ray i, summing over the other two tail rays gives the
joint profile cost sum_(t,z) G_j[i,t] G_a[i,z] G_b[t,z]. The cost uses the
same body profile for every t, so sum first, then take the LP majorant.
The provided nonnegative rational multipliers give an upper bound on this
cost after exact fraction-free row reduction. Any uncovered residual is
bounded using 0≤profile count≤shell population; no floating optimizer
output is a premise. This gives at most 2,894,818,655,349,910,178
triangles through a graph vertex (same-body terms also bounded explicitly).
The validated blocked Boßecker–Rautenbach recurrence yields an independent
set of at least 459,824 vertices. Its elements are not materialized.

The same oriented construction uses 13 integer rays and a pure A2 tail in
74D, and the 60 oriented icosahedral rays and a pure D3 tail in 75D.
For 77–85D the bound uses signed permutations of the parent integer tail
shapes, with squared body scale X from the parent table and squared tail
length 8−10X. It uses one representative per **body** line, testing absolute
products, and outputs both signs of each accepted head-plus-tail pair.
The tail sets and their stabilizer orbits are explicitly regenerated.
When counting third-body candidates in this graph, symmetry of the tail
palette makes the summed profile cost invariant under (a,b)→(−a,−b);
we check this exactly and divide the oriented joint count by two.
The graph calculation and final integer ceilings respect these differing
orientation conventions. Pure D_k roots are added in dimensions 75–85;
the exact cross-tail angular constraint is checked for each shape.

Source inputs: [Nebe's minimum-eight lattice](https://arxiv.org/abs/1008.2862);
[Nebe's exposition of Venkov's weighted theta theory and shell designs](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf);
[the parent shell arithmetic and independence recurrence](../dim77-85-delta-marginal/).
Run python verify.py with NumPy, SymPy, and assertions enabled. All geometry,
linear algebra, and summation in the replay are exact. Supply dimension numbers
as arguments to restrict the requested output rows; the same common
harmonic-system proof is replayed in either mode.
