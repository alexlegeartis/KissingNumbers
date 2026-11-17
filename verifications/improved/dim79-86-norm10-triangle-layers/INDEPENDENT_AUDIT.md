# Independent norm-ten triangle / independence audit (79–86D;82D counts illustrated)

**PASS**, using the explicit published Gamma72 existence/minimum and shell-design theorems. The frozen 40-block result proves existence of at least **6,221,706,936 kissing points in dimension82**. This is an existence proof, not a supplied list of selected heads or a lattice packing construction. No worldwide-priority assertion is made.

## What was independently checked

* Independently enumerated all3360 signed weight-four tails and their neighbor orbit multiplicities.
* Reconstructed Gram-positive integer profile domains for all three source products4,5,6, retaining exactly the necessary small-combination lattice restrictions.
* Reimplemented mixed spherical moments via a closed combinatorial Gaussian pairing formula rather than the producer's recursive Wick function.
* Verified all nine rational LP dual majorants with36 moment rows plus26 exact marginal rows on every allowed profile, including exceptional vectors and their moment corrections. No floating solver optimality or feasibility output is trusted.
* Recomputed exact regular degree26,035,294,328, graph order25,673,404,435,660,800 and local triangle bound40,824,250,566,624,574,710.
* Recomputed all40 rational recurrence blocks, the ceiling1,765,578 selected antipodal pairs, and final kissing count.
* Verified scaling/equatorial/pure-tail compatibility using exact squared inequalities. Everywhere “norm8” means squared Euclidean norm8.

`verify_one.py` is the independently written checker. Each `audit_resultN_marginal.json` pins its candidate hashes, and `replayN_marginal.log` contains its run. Source theta and univariate moment functions were separately replayed from the existing public norm10 package; this source replay is not a new implementation (`source_replay.json`).

## Graph normalization and upper bounds

A graphvertex is an antipodal point pair indexed by a norm-ten antipodal source line and a signed tail. Oriented source vectors and signed tails represent each vertex twice, by (u,t)~(-u,-t), explaining the factor1/2 in common-neighbor counts. Neighbor source lines of absolute productj are the source's oriented absolute-product population divided by2. Changing the chosen orientation permutes the signed tail set, so no global orientation assumption is needed.

Every same-source pair conflicts, giving a clique of3360vertices. All other edges use exactly |8a+6d|>52. Products at most3 cannot conflict. The source product distribution and signed-tail transitivity make the graph regular; graph automorphism transitivity of Gamma72 is NOT assumed.

The common-neighbor LP may count the two endpoints as neighbors of themselves, but that only enlarges the common-neighbor upper bound. Half the sum over neighbors bounds triangles through a fixed vertex. Same-source contributions use the safe boundD. The final T is an INTEGER rounded upward, as the independence theorem requires. A stale noninteger-T candidate was identified and replaced before this audit; frozen hashes identify the corrected data.

## Independence theorem and recurrence

The applicable result is Boßecker–Rautenbach (2010), not the triangle-free-only Shearer theorem. Its hypothesis is at most integerT triangles through EVERY vertex, and it bounds independence by the sum of f(degree,T). The theorem and recurrence were checked against the author's [primary abstract](https://www.imada.sdu.dk/u/btoft/GT2009/Rautenbach.pdf) and [published article](https://doi.org/10.1016/j.disc.2010.05.026).

For g_d=(d+1)f(d,T)−1, the producer's formulas A_d=1−2(d+T)/(d(d²+1−2T)) and B_d=(1−(d+1)/(d²+1−2T))/d are algebraically correct. Above threshold, both are nonnegative. On each block[a,b], the rational bounds A_d≥1−C and B_d≥B follow because (d+1)/(d²+1−2T) and (1+T/d)/(d²+1−2T) decrease. The bound on g_b uses Bernoulli for the retained initial term and the finite geometric sum; both inequalities have the correct direction. The initial g≥0 follows from the recurrence and its Caro–Wei baseline. All blocks are contiguous and reach D.

## Source premises and kissing validity

[Nebe's construction](https://arxiv.org/abs/1008.2862) supplies the even unimodular72D lattice of squared minimum8. [Nebe's Venkov survey](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf), Theorem4.1 and its discussion, supplies strength11 for its shells. These published theorems are premises, not re-proved by this audit. The theta/moment calculations regenerate N8,N10 and the required source product counts.

The new points have squared norm8, and independent graphvertices ensure all their mutual inner products are at most4, including antipodal signs. Against the equator, the squared maximum product is25·8/13<16; against pure D10 roots it is16·6/13<16. Therefore all retained/additional points form a legitimate spherical kissing configuration after common rescaling. The three kinds of points have distinct head/tail structure, so none coincide.

## All eight final dimensions

The same independent checker separately passed79–86D with the frozen marginal certificates. The26 indicator populations follow from the source univariate product counts, with the four exceptional vectors removed. All final duals use the unscaled basis; the reader need not reproduce numerical row scaling. A separate checker, `check_recurrence.py`, independently rederived and checked all eight blocked recurrences, including their actual dimension-specific starting endpoints.
