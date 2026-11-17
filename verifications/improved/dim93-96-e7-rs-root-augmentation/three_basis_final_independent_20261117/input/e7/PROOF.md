# Three-basis E7-head construction in95D

This package uses the same three exact four-dimensional alphabets and inherited Z8 parent minimum24 as `three_basis_z8_93_96_20261117`. That construction has independently passed, including fresh replay of the inherited parent proof. All new head geometry is checked here directly; discovery by Weyl reflections is not a proof dependency.

## Source section and head tables

The first two ring symbols a,b determine two signed crosspolytope vectors2*(-1)^floor(a/4)*e_(a mod4) and2*(-1)^floor(b/4)*e_(4+b mod4) in R8. Retain the38 prefixpairs whose sum lies in the hyperplane x0+x4=0. Their heads are36 cross edges between axes1,2,3 and5,6,7 with all signs, plus the two heads+/-2*(e0-e4).

The certificate supplies three38x8 head tables indexed by these same source prefixes. Their complete Gram matrices are exactly equal to the original two-row Gray Gram matrix. Every head has squared norm8 and lies in x0+x4=0. They were found through exact E7 Weyl reflections, but the verifier only needs their explicit entries and pairwise products.

Their union contains100 distinct heads. Assign each repeated head to its earliest table and discard its later occurrences, leaving38,32,30 selected heads in the three colors. `candidate.json` lists all100 selected color/prefix/head records. Every pair of selected heads from different colors has inner product at most4, verified exactly. No optimality of the100-head union is asserted.

## Full points and separation

For each selected prefix, use every word in the corresponding source Z8 fiber, retaining its last22 ring coordinates. Map those coordinates through the even Gray sign alphabet, the odd sign alphabet, or the coordinate alphabet according to the head color, then append the assigned head. Each fiber has8^10 words by the source generator's invertible first-two-coordinate projection.

The retained body has squared norm88 and the head squared norm8. Within each color, exact prefix-Gram preservation proves the total distance is the source distance, at least96. Between colors, every bodyrow has cross inner product at most2, so bodyinner<=44; headinner<=4; totalinner<=48. Hence all vectors form a kissing configuration after normalization.

Bodies lie in R88 and heads in the seven-dimensional hyperplane x0+x4=0 in R8. The whole configuration therefore lies in a95-dimensional Euclidean space. An explicit orthonormal head basis is (e0-e4)/sqrt(2),e1,e2,e3,e5,e6,e7. No numerical projection is needed.

Distinct heads distinguish all selected prefixes and colors. Within a fixed prefix, deleting it is injective. Thus the exact count is

    100*8^10 = **107,374,182,400**.

No sparse layer or poles are included, and no compatibility with earlier sparse layers is assumed. No publication or worldwide-priority claim is made.

## Replay

Run `python verify.py`. It checks generator uniformity, all three alphabet tables, all three complete source/head Gram matrices, the100 selected records, the hyperplane, all cross-color products, injectivity and count. The only inherited premise is the already audited Z8 minimum24. Independent audit is requested separately.
