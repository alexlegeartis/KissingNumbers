# Independent mixed-latitude 89D audit

Exact kissing lower bound: **13,958,643,714 = 13*8^10 + 2**. This improves the live 89D embedding bound 12,884,901,888 by 1,073,741,826, approximately 8.33333 percent. The earlier 12,690,155,704 construction is not the current comparison.

Use the same source-bound Z8 parent and three alphabets as the enclosing package. Keep prefixes (0,b), b=0,1,2,3 in colors 1 and 2. In color 0 keep b=0,1,2 with height zero, b=3 with height sqrt(8/3), and b=7 with height -sqrt(8/3). Delete the first two original blocks. The body has norm squared 88. The eleven zero-height fibers have norm squared 88, and the other two have norm squared 272/3. Normalize **each point by its own norm**, not by a common norm.

## All distinct pair types

The inherited homogeneous minimum 24 gives full inner product at most 48. Removing the prefix leaves:
- Within a fiber: body product at most 40. Adding the height squared gives 40 or 128/3, respectively, no more than half of the corresponding norm squared (44 or 136/3).
- Between the eleven zero-height fibers: at most 44 (55 fiber pairs), hence normalized product at most 1/2.
- Between a raised fiber and a zero-height fiber: at most 44 (22 pairs). The exact squared comparison is 4*44^2 <= 88*(272/3).
- Between the two raised fibers: their deleted prefix product is 4-4=0, so body product may reach 48. Their opposite heights subtract 8/3, leaving 136/3, exactly half of 272/3 (one pair).

Cross-color body bounds use the 22 blockwise products, each at most 2. Same-color bounds use the exact prefix products. This covers all 78 inter-fiber pairs and all 13 within-fiber types without sampling or numerical tolerance.

## Population and injectivity

The triangular odd-pivot generator is injective over Z8 and its first-two-coordinate projection is onto Z8^2; all 13 distinct labeled prefix fibers have exactly 8^10 members. Equality after pointwise normalization would yield normalized product 1; every pair of distinct indexed points has product at most 1/2. Thus neither puncturing nor unequal-radius normalization introduces duplicates, even across different colors.

The new independent checker reconstructs the generator-rank argument, prefix map, alphabets and homogeneous metric identity; it checks all cases with integers and rational numbers. It binds the generator by SHA256 to the fully replayed parent. The complete enclosing verify.py replays the inherited minimum proof and both geometry checks. Reusing the independently authored parent algorithm is explicitly distinguished from writing a second parent-minimum algorithm.

The external-source table comparison remains the derived Gamma72 plus 17D tail bound 6,218,181,330. Limited public searches did not identify a stronger external 89D construction, but no exhaustive worldwide-priority or global-optimality claim is made.

## Two pure poles

Add the two unit vectors along the new 89th coordinate, with both signs. They have mutual product -1, product zero with all zero-height fibers, and product magnitude 1/sqrt(34) < 1/2 with each nonzero-height fiber. The exact comparison is 4*(8/3) < 272/3. The nonzero body of every fiber point excludes equality with a pole. Thus the total is 13*8^10+2.
