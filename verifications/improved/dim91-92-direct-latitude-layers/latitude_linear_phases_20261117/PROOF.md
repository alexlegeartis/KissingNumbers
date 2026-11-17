# Stronger exact linear phase codes for direct91D/92D latitude layers

This package replaces only the phase codes in the independently audited `latitude_sparse14_20261117` and `latitude_grid92_20261117` constructions. Their outer arrays, heads, equatorial points, and every geometric argument are unchanged. The prior phase codes are not retained as additional layers.

Primary sources:
- https://www.codetables.de/BKLC/BKLC.php?q=3&n=14&k=8
- https://www.codetables.de/BKLC/BKLC.php?q=3&n=13&k=7

Grassl's tables specify a ternary constacyclic [14,8,5] code with shiftconstant2 and generator polynomial g(x)=x^6+x^5+x^3+x+1, and its shortening at the last coordinate as [13,7,5]. These claims are proved finitely here, without relying on the table's distance assertion.

Let G have rows the coefficient vectors of x^i*g(x),0<=i<8, in positions0..13. The first8 columns are upper triangular with diagonal1, proving rank8. Enumerating all6561 coefficient tuples gives the complete linear code. The only weightzero word is0; weights1..4 are absent; weight5 occurs140times. Therefore every distinct pair has distance at least5, by linearity.

Only the eighth row has a nonzero last coordinate. The shortened code therefore has generator G[0:7,0:13], dimension7, and2187 words. Enumeration gives minimum nonzero weight5, occurring90times. Both full exact weight distributions are recorded in verification.json. The verifier also checks g divides x^14+1 in F3[x].

Apply the14-symbol code to each of the31,680 outerwords in each weight1-prefix latitude layer, and the13-symbol code to each of the6,930 outerwords in each weight2-prefix corner layer. Distinct layers may use independent phase codes because the already audited cross-layer bounds depend only on their outerwords and heads.

Exact counts:
- Each weight1 layer:31,680*6,561=207,852,480.
- Each weight2 layer:6,930*2,187=15,156,660.
-91D:21,943,640,934 equatorial points +2*207,852,480 latitude points +2 poles = **22,359,345,896**.
-92D:21,943,640,934 equatorial points +4*207,852,480 weight1 latitude points +4*15,156,660 corner latitude points +6 orthogonal hexagon points = **22,835,674,500**.

The substantial new-coordinate gains are415,704,960 in91D and892,033,560 in92D. Inherited premises remain the Pless parent minimum15 and validity of the audited equator. No publication is performed by this package.
