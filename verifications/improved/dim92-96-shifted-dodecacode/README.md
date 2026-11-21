# Shifted dodecacode phases in dimensions 92–96

## Exactly verified continuation

| Dimension | Latest published bound checked | Replacement bound | Gain |
|---:|---:|---:|---:|
| 92 | 30,141,844,072 | 30,151,281,256 | 9,437,184 |
| 93 | 43,988,421,272 | 44,038,752,920 | 50,331,648 |
| 94 | 74,053,323,468 | 74,103,655,116 | 50,331,648 |
| 95 | 110,560,545,860 | 110,610,877,508 | 50,331,648 |
| 96 | 208,869,002,048 | 209,198,254,912 | 329,252,864 |

Replace the entire prior Golay–dodecacode augmentation in each rank. The earlier dense, six-support and root layers are retained. The new construction and the previous Golay–dodecacode layer are not asserted to be mutually compatible. Geometry from the earlier witness still applies to the earlier dense, six-support and root layers. The new phase choices leave the block support, tag weights, norm and coordinate magnitudes unchanged. This is an exact lower-bound construction over inherited published proofs. [Earlier construction](https://github.com/alexlegeartis/KissingNumbers/tree/dd7b8d2124086954e557a90ea218a2a70577e580/verifications/improved/dim92-96-golay-dodecacode).

## Three independently assignable phases

The earlier layer used 12 occupied rows with two entries of magnitude 2 per row. A quaternary weight-12 pattern prescribed the occupied rows and one of three pair matchings per row. Ordinary compatibility required outer distance at least 12. A pattern built from two binary Golay words a,b has three nonzero tag sets T1=a minus b, T2=b minus a, T3=a intersect b. For patterns with a,b,a+b all octads, these three sets are tetrads. Each union of two is an octad. The only other pattern type used is a dodecad with all 12 symbols the same.

Selectively allow outer distance 8 between two tetrad patterns with exactly two tag sets identical and the remaining tetrads disjoint. The previous worst-case dot for these pairs was 64, from the eight common rows all with the same tag. It falls to at most 48 if their two axis words differ in at least two of those eight rows.

Order inner coordinates by ascending physical row within T1, T2 and T3, in that order. Then each common pair of tag sets always aligns two fixed groups of four inner coordinates. The additive dodecacode has projection rank 12 on each of these three eight-coordinate sets. Its quotient there has dimension four. Each quotient admits four translates at mutual Hamming distance at least 2. In quotient coordinates calculated by binary elimination with highest pivots, they are [0,5,10,15] when tag 1 or 3 is absent, and [0,5,9,12] when tag 2 is absent. The syndrome map from an arbitrary 12-coordinate shift to the three projected quotients has binary rank 12. Thus a single translation of the 12-coordinate dodecacode independently realizes any three requested projected cosets.

For each pair of retained tetrads, there are exactly four possible completions in the ambient Golay pattern pool. Indeed the two retained tetrads have an octad as their union. Either retained tetrad lies in five octads, one of which is that union, and the other four intersect it in exactly this tetrad. The four possible third tetrads are distinct and disjoint. Counting directly confirms the five-octad property for each of the 10,626 tetrads. Sort the four possible third tetrads numerically as 24-bit masks. Use that index as the phase color on the corresponding pair of retained tag sets. Two distance-8 patterns of the allowed type receive different colors in their shared pair. Two differing coordinates of their axis projection then contribute zero there because the common row frame is orthogonal, giving dot at most 48 regardless of signs. The translation leaves same-pattern axis distances unchanged. Sign code and count per pattern remain 4,096 times 256.

The verifier checks shifts, all binary ranks, the four-completion property, all pattern pairs, geometric row bounds and inner code sizes. The selected lists have 337, 1,037 and 2,897 patterns at masks 15, 3 and 0. There are respectively 32, 133 and 478 rescued distance-8 pairs. The same mask-3 list works in 93–95D. The full 96D construction adds 3,037,724,672 points to the dense/sparse predecessor. A separate larger outer distance-12 search reached 2,736 patterns, so the phase construction is still larger than the best unshifted witness in this pass. The pattern optimization is not a premise of the proof.

## Replay

Run `node verify.js`. The base-36 pattern lists record pairs of 12-bit Golay messages. The verifier independently computes projected quotients, finds the translations, checks every relaxed distance-eight pair and all other outer pairs, and repeats the finite geometry checks. The unchanged cross-layer proof is pinned in the earlier construction; the superseded augmentation is not added again.
