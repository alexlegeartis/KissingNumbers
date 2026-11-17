# Nine-prefix signed-grid construction in dimension92

This extends, without modifying, the verified `latitude_sparse14_20261117` package. Keep all its92D points, including the equatorial configuration, four weight1-prefix layers, and optional six orthogonal hexagon points. The Pless parent minimum15 and audited equator validity remain explicit premises.

Add the four deleted prefixes110,120,210,220. Each supplied array contains6930 distinct Pless codeword bodies of weight13. Give them heads sqrt(3)*(f(p1),f(p2)), where f(0)=0,f(1)=1,f(2)=-1. These cornerheads have squared norm6; the bodies use thirteen signed A2 roots and have norm squared39. Total norm squared45.

Use the supplied phase13 code, obtained by fixing coordinates0,3 of the explicit7812 phase code to0,2 and deleting them. There are888 words, all pairs have minimum distance5, checked directly. Each corner layer therefore has6930*888=6,153,840 points.

General cross-layer proof: for distinct prefixes p,q of weights w,v, the full Pless code minimum15 implies body Hamming distance at least15-d_H(p,q). Thus the body squared distance is at least45-3*d_H(p,q), using the same signed A2 alphabet. The signed-grid heads satisfy exactly

    ||h_p-h_q||^2 = 3*sum_i(f(p_i)-f(q_i))^2 >= 3*d_H(p,q).

The combined squared distance is at least45. This includes the weight0 equatorial sparse15 layer and the four prior weight1 layers, regardless of their different phase codes. Within a corner layer, different outerwords have body Hamming distance at least15, and identical outerwords with different phasewords have at least five root differences of squared distance9. Therefore all within-layer squared distances are at least45.

The dense equator cross inner product has absolute value at most13*(3/2)=19.5. The old orthogonal regularhexagon has norm sqrt(45) and the largest head norm sqrt(6), so absolute cross inner product is at most sqrt(270)<22.5. Different head locations ensure layer disjointness; root alphabet injectivity ensures internal counts.

Four new corner layers add24,615,360 points, beyond the330,612,480 points in the preceding four square layers. Total direct92D gain over the90D equator is355,227,840, plus6 optional orthogonal points. The final92D bound is **22,298,868,780**.

The91D result is unchanged: **22,108,947,176**.
