# Complementary tail directions above the signed P48 configurations

| Dimension | Preceding certified bound | Added | Removed | New lower bound |
|---:|---:|---:|---:|---:|
| 54 | 52,976,623 | 162 | 6 | **52,976,779** |
| 55 | 53,393,546 | 168 | 42 | **53,393,672** |
| 57 | 54,760,693 | 288 | 24 | **54,760,957** |
| 59 | 56,952,414 | 192 | 24 | **56,952,582** |
| 60 | 58,631,571 | 3,584 | 1 | **58,635,154** |

Keep the [signed P48 class construction](../dim49-60-p48-signed-refill/)
and add directions using unused signed minimal heads. The head and tail
spaces are orthogonal. For an unowned unit P48 minimum direction u and
unit tail direction z, the basic inserted pair is
(±u/2, √3 z/2). It is compatible with the full equator. If a tail direction
has absolute inner product at most 1/√3 with the old tail code and pairwise
inner products at most 1/3 with other new tails, one such pair can be
inserted for every new tail using the same unowned head line.

- In 55D, the 56 tails are the signed even-parity six-coordinate patterns
  and 24 directions obtained from one of the first six axes and the
  difference of the final axes in the actual E7 presentation. Fourteen
  colors carry signed A2 head trines, each with its three equator copies retired.
- In 57D, the 12 unused triples of coordinates in the 9D tail form an
  affine plane of order three. All eight signs on each triple give
  96 tails. Different triples meet in at most one coordinate. The four
  parallel classes, split by sign parity, give eight colors. On each
  color use an unowned signed A2 head trine at the same latitude as
  in 54D, retiring its three equator copies. Heads are disjoint across colors.
- In 59D, eight explicit triples with intersection at most one give
  64 tails. An eight-color partition uses eight disjoint signed A2 trines.
  Literal products against the old code are checked in Q(√2).
- In 60D, replace the old 841 poles with a literal antipodal 840-code,
  preserving all lifted classes by a bijection. Partition the even 12-bit
  cube by XOR syndrome and assign fourteen classes to orthogonal heads.
  See the [60D construction](d60cube/).
- In 54D, use a different latitude: (u/√3, √(2/3) z).
  Project the 27 roots on one side of an E7 cut to the E6 subspace, giving
  the minuscule family. Same-side normalized tail products are at
  most 1/4; opposite-side pairs are at most 1/2. A single A2 triple
  of unowned minimum heads and its negatives label the two sides.
  Retire the six equator copies of those heads. The head trine
  handles the three rays at each tail; the 162 new rays are all
  distinct. The old E6 tail products hit at most √(3/8) in absolute
  value, so the mixed cap is exactly 1/2.

The replay pins the input Gram, exact tail data and signed ownership
files. New mixed and pair inequalities are checked by bounded integers,
rational comparisons and exact quadratic signs. No heuristic candidate
generator is used in replay. The already checked geometry of the
unchanged predecessor and the published P48 minimum/population
are inherited premises. The improvements are over the preceding
repository constructions. For external comparison, Kamenetsky,
[arXiv:2611.01906v1, Table 3](https://arxiv.org/html/2611.01906v1)
reports lower values in all five dimensions.

Run python3 verify.py in this package, with NumPy and assertions enabled.
