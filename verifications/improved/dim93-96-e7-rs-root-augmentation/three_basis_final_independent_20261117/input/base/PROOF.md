# Three real crosspolytope alphabets from the audited Z8 parent

## Inherited premise and explicit alphabets

Let C be the frozen free rank12 Z8 code of length24 in `dim88-96-z8-golay`, whose Gray sign image has squared norm96 and minimum squared distance96. The homogeneous minimum24 is inherited from its previously audited proof; no new code or source theorem is assumed. The supplied generator has odd upper-triangular first12 columns, so |C|=8^12, and each first-two-coordinate fiber has size8^10.

For symbol a in0,...,7 define A(a) as the four signs of its usual affine Gray word on F2². Define B(a) by negating only the first physical coordinate of A(a). Define D(a)=2*(-1)^floor(a/4)*e_(a mod4). All three alphabets have the same8x8 Gram matrix: norm squared4, antipodes at symbols differing by4, and all other pairs orthogonal. For any two DIFFERENT alphabets, every cross inner product has absolute value2. All statements are verified by direct integer multiplication of the displayed matrices.

## Full96D construction

Map every ring word of C coordinatewise through A, through B, and through D, forming three copies. Within a copy, Gram preservation gives the source minimum squared distance96. Between copies, the inner product is at most24*2=48, while each vector has squared norm96. The copies are disjoint: A and B have different binary parity in every row, while D has one nonzero physical coordinate per row rather than four. Thus the union has

    3*8^12 = 206,158,430,208 points in R96.

The third copy is a new explicitly code-controlled24-support family. No earlier sparse family is retained without fresh compatibility checks.

## Compressed93–95D constructions

Retain22 rows, in R88. For head dimension d=5,6,7, label physical head axes0,...,d-1 by their three-bit expansions. Partition edges of K_d by the first bit in which their endpoint labels differ. For color c, let X be the axes with that bit0 and Y those with that bit1. Both lists have at most4 axes. Map the first ring symbol to its signed axis in X and the second to its signed axis in Y, selecting exactly the edge pairs assigned that color. The head is2*(+/-e_x+/-e_y), with norm squared8. Use alphabet A,B,D on the retained22 rows for colors0,1,2 respectively.

Within a color, the two row-axis lists are disjoint, so every prefix Gram is preserved exactly. The full new distance therefore equals the original source distance, at least96. Between colors, body inner product is at most22*2=44. The selected head supports are different edges; distinct D_d roots have inner product at most4. Thus total inner product is at most48. All points have norm squared88+8=96.

Every edge permits four signs, so every D_d root is used exactly once. Each fixed prefix contains8^10 source words. Deleting a fixed prefix is injective; different prefixes/colors have different heads. Therefore the exact dense counts are:

-93D:2*5*4*8^10 = **42,949,672,960**.
-94D:2*6*5*8^10 = **64,424,509,440**.
-95D:2*7*6*8^10 = **90,194,313,216**.

These counts include no sparse additions or poles. In particular, the previously constructed sparse22 layer is NOT presumed compatible with the third coordinate alphabet. No92D improvement is claimed here. No publication or worldwide-priority assertion is made.

## Exact replay

`verify.py` checks the three complete alphabet tables, source generator uniformity, explicit edge partitions and prefix/head maps, every within/cross-color head pair, uniqueness and all counts. `certificate.json` contains every selected prefix and head. The inherited source minimum24 is clearly separated from these fresh finite checks. Independent audit is requested separately.
