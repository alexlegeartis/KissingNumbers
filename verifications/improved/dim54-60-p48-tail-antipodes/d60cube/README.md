# 60D: a 12-cube refill with fourteen orthogonal P48p heads

This proves **58,635,154** (+3,583 over 58,631,571). From this directory,
run `python verify.py --siblings ../..` to check the pinned ancestor,
the integral orthogonal frame, the replacement tail code, and the sign
partition without optimization.

Keep the 420 disjoint signed head classes of the old construction.
Reassign their antipodal tail pairs bijectively using a literal 840-point
antipodal contact code. This costs one pure pole (840 in place of 841),
and costs no lifts: the old lifting rule only needs 420 antipodal pairs
from a 12D kissing code.

The 840-point code splits the 12 coordinates into two six-blocks. Take
the 24 signed axes, the sign cubes on the four-subsets complementary to
one of three fixed matching edges in each block, and the sign cubes on
four-subsets formed by joining equally colored edges in the two standard
one-factorizations of the complete graph on six vertices. The resulting
51 supports and 24 axes give 840 points. The verifier uses integer dot
products. Every unit point has ℓ¹ norm at most 2.

Partition the even 12-bit sign cube by the XOR of indices 0 through 11
at which the word is negative. The 16 classes have 128 words each.
Words within a class differ in at least four positions; words across
classes differ in at least two. Assign 14 classes to the 14 literal
orthogonal unit heads. For each head u and normalized sign word z in
its class, add the two points (±u/2, √3 z/2).

- The equator bound follows from the new head amplitude 1/2.
- Against the reassigned poles, the tail product is bounded by
  2/√12 before scaling by √3/2, giving 1/2.
- Both signs of every new head are outside all lifted signed classes,
  so the unit head products with lift heads are at most 1/2.
  Against an old lift the resulting bound is √3/8 + 1/4 < 1/2.
- With the same tail, opposite head signs give inner product 1/2.
  With different tails in one class, the bound is 1/4 + (3/4)(1/3) = 1/2.
  Across classes the heads are orthogonal and the bound is (3/4)(2/3) = 1/2.

All original layers have different head amplitudes from the appended layer.
The append has 14 × 128 × 2 = 3,584 distinct points.
