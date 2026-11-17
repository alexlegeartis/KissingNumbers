# Second independent audit of the 82D blocked independence recurrence

PASS, conditional on the source graph, triangle bound and Boßecker–Rautenbach theorem, which the separate full auditor checks. Reimplemented exact arithmetic and symbolic recurrence algebra in check.py; did not execute the producer's block_recurrence.py.

Write H=d²+1−2T. Substitution g_d=(d+1)f(d,T)−1 gives exactly

    g_d = [1−2(d+T)/(dH)] g_(d−1) + (d²−d−2T)/(dH).

At and below the threshold T>=binom(d,2), g_d=0 by the theorem's base branch. Above it both recurrence coefficients are nonnegative, proving g_d>=0 inductively. No unjustified inference from an independence bound to a pointwise property of f is needed.

For a block [a,b], the verifier checks a(a−1)>=2T, positive denominator, B>=0 and 0<C<1. The derivative of (d+1)/(d²+1−2T) is −(d²+2d−1+2T)/(d²+1−2T)²<0 for d>=1,T>=0. Its complementary numerator is therefore nonnegative and increasing. Dividing by d<=b proves B_d>=B_a*numerator-normalization a/b, exactly the producer B. The quantity (1+T/d)/(d²+1−2T) is decreasing as a product of nonnegative decreasing factors, proving A_d>=1−C. These comparisons preserve lower bounds because g>=0.

Unrolling over L steps gives g_b>=(1−C)^L g_(a−1)+B sum_(j=0)^(L−1)(1−C)^j. Bernoulli supplies (1−C)^L>=max(0,1−LC). Conversely (1−C)^(-L)>=(1+C)^L>=1+LC, so the geometric sum is at least L/(1+LC). Thus the recorded block lower bound is rigorous for every block, including LC>1 where the previous g contribution is discarded rather than multiplied by a negative number.

All40 contiguous blocks from16,141,882,484 through26,035,294,328, all exact B,C values, the accumulated rational factor and expected number of independent vertices match the frozen certificate. Ceiling the rational expectation is legitimate because alpha(G) is an integer. It yields1,329,536 antipodal pairs and the total6,220,834,852.

Frozen blocked82.json SHA256:0799fc0895493b0a14a9bf6a52b795acb0260936330cf885384766c19b2a5853. Input digests and arithmetic report are in RESULT.json. This report does not independently establish the source graph/populations or the external theorem; combine with research/norm10_triangle_independent/REPORT.md.
