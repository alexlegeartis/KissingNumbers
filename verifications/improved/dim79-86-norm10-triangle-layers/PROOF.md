# Triangle-sparse norm-ten tail layers (79–86D verified existence bounds)

The final marginal certificates passed independent exact replay and a separate recurrence review. The82D example gives6,221,706,936 points; the dimension-specific JSON files supply all79–86D values. These are existence results, not coordinate lists or global-priority claims.

## Inherited inputs
Nebe's even unimodular Gamma72 has minimum8, minimal shell6,218,175,600, norm10shell N=15,281,788,354,560; the norm10shell is a spherical11-design. Its absolute product4,5,6 populations around each norm10vector are respectively24,994,840,000;500,966,912;2,112,250. These are inherited from the public KissingNumbers norm10 proof at `verifications/improved/dim73-88-gamma72-cosets/source/norm10/`, whose theta/moment replay and literature inputs must be retained, not replaced by this search output.

## Graph and kissing geometry
Take the M=16 binom(10,4)=3360 signed weight4 tails t in R10. For each antipodal norm10line {v,-v}, choose an orientation and each t, and form the antipodal point pair ±(sqrt(8/13)v,sqrt(6/13)t). Every point has norm8. Vertices are these pairs; changing the chosen line orientation only relabels t to -t. Edges are precisely incompatible pairs, characterized by |8(v.u)+6(t.s)|>52. Same-source vertices form a clique: even opposite tails give80-24=56>52. Distinct source lines have absolute product≤6; products≤3 never conflict. The resulting graph is regular, with N M/2 vertices and degree26,035,294,328.

All retained Gamma72 norm8 equatorial points are compatible because |v.u|≤5 for norms10and8 and 5sqrt(8/13)<4. The180 pure D10 roots 2(±e_i±e_j) have norm8 and are compatible with the signed tails, since their products have absolute value≤4sqrt(6/13)<4. Thus an independent set of A graphvertices gives6,218,175,600+2A+180 kissingpoints.

## Local triangle upper bounds
Fix a vertex with tail r=(1,1,1,1,0,...,0). For a neighboring source line choose orientation w so v.w=j>0, j=4,5,6. Neighbor tail s is classified by its support overlap c with r and dot d=r.s. Coordinate permutations inside/outside r and sign changes outside r prove common-neighbor bounds depend only on(c,d). Multiplicities are counted by enumerating all3360tails. Such a tail is a neighbor exactly when8j+6d>52. Source-line multiplicity at absolute product j is the published shell population divided by2.

For an oriented third norm10vector u, let a=u.v,b=u.w. Except four vectors ±v,±w, a,b are integers in[-6,6]. Gram positivity gives10a²−2jab+10b²≤10(100−j²). Every integer combination xv+yw+zu used in the finite screen has squared norm either0or at least8. We impose this for x,y,z∈[-2,2] with z nonzero and x or y nonzero; omission of other necessary restrictions only relaxes the LP.

All moments a^r b^s, r+s even≤10, follow from the11design. Gaussian pairing moments W satisfy W(r,s)=(r−1)10 W(r−2,s)+s j W(r−1,s−1) (and the analogous one-variable rule). The shell moment is N 10^m W(r,s)/(72·74···(72+2m−2)), m=(r+s)/2. Remove the known four exceptional vectors from each moment; the LP variables count all remaining oriented u profiles.

For each profile(a,b), the objective is one half the number of signed tails t satisfying BOTH |8a+6r.t|>52 and |8b+6s.t|>52. Division by2 quotients the involution(u,t)→(-u,-t). Exceptional vectors are added explicitly. Counting graph endpoints as potential common neighbors can only overestimate, which is harmless.

LP dual coefficients are rationalized. Any coefficientwise shortfall on the finite profile domain is repaired by adding its maximum to the constant polynomial coefficient. This produces an exact polynomial majorant of the objective. Its moment expectation plus exceptional contributions gives a rigorous common-neighbor upper bound conditional on the inherited inputs. No floating solver optimality claim is used.

Sum common-neighbor upper bounds over all neighbor classes, and add(M−1)D for same-source neighbors using the trivial upper boundD. Half this sum, rounded upward, bounds triangles through every graphvertex. The `exactN_marginal.json` files store all rational repaired duals and the result.

## Independence theorem and finite rational evaluation
Boßecker–Rautenbach's interpolation between Caro–Wei and Shearer states α(G)≥sum_v f(d(v),T) when every vertex belongs to at most integerT triangles. It defines f(d,T)=1/(d+1) when T≥binom(d,2), and otherwise

f(d,T)=[1+(d²−d−2T)f(d−1,T)]/(d²+1−2T).

Primary author abstract: https://www.imada.sdu.dk/u/btoft/GT2009/Rautenbach.pdf . Full publication: Interpolating between bounds on the independence number, DOI10.1016/j.disc.2010.05.026. The theorem statement and hypotheses were independently checked against the primary source.

Put g_d=(d+1)f(d,T)−1. Then g_d=A_d g_{d−1}+B_d, with
B_d=(1−(d+1)/(d²+1−2T))/d,
A_d=1−2(d+T)/(d(d²+1−2T)).
Caro–Wei and the recurrence give g_d≥0. On a degree block a≤d≤b lying above the threshold, let
B=(1−(a+1)/(a²+1−2T))/b,
C=2(a+T)/(a(a²+1−2T)), L=b−a+1.
Monotonicity of (d+1)/(d²+1−2T) and (1+T/d)/(d²+1−2T) gives B_d≥B and A_d≥1−C. All blocks check B≥0,0<C<1. Therefore

g_b≥max(0,1−LC) g_{a−1}+BL/(1+LC).

Indeed(1−C)^L≥max(0,1−LC), while its geometric sum is≥L/(1+LC), by (1−C)^(-L)≥1+LC. Forty contiguous blocks starting at the dimension-specific endpoint recorded in the JSON (above sqrt(2T)), with initial lower boundg=0, yield a rational lower bound on f(D,T). Multiply by the exact graph order and take the ceiling. The `blockedN_marginal.json` files record every block;82D yields1,765,578 antipodal pairs.


## Exact marginal strengthening used in the final files

In addition to the36 normalized joint monomial rows, append26 indicator rows: first a=−6,…,6, then b=−6,…,6. The exact population for product0 is the full absolute-product-zero population; for each nonzero signed product it is half its absolute-product population. Subtract the four exceptional vectors from each marginal, just as for the moment rows. These distributions are independently regenerated by the inherited univariate shell-moment proof.

The stored duals use this unscaled row basis. Row scaling was used only during numerical discovery and was algebraically removed before freezing. Exact verification evaluates every dual on every allowed profile and evaluates its objective using the exact unscaled moment and marginal right-hand sides. Thus numerical conditioning and solver tolerances supply no theorem premise. Adding these valid equalities can tighten common-neighbor bounds without changing the graph or the independence theorem.
