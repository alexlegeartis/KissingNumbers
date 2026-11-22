# 60D: a nested sign cube and sparse tetrahedral tails

This exact tail module gives **58,636,594** over the pinned 58,631,571
signed P48p ancestor. It also composes with a separate owner exchange if
the exchange does not use any of the marked head lines.

Start with the [classical 840-point antipodal palette](../dim54-60-p48-tail-antipodes/d60cube/),
used for both poles and old lift directions. This retains the 420 old
classes and costs one pure pole. Use a literal 16-frame, rather than 14
orthogonal heads, for the even 12-bit sign cube. Assign the 128 words in
each XOR-syndrome class to one head, and use both signs of that head as in
the preceding construction.

Two changes use tetrahedra of unit minimal P48p heads. A tetrahedron
has distinct-head products -1/3; all four heads sum to zero.

1. **Sparse tails.** There are exactly 16 triples not contained in the
   51 four-supports of the 840 palette. In each six-block choose one
   coordinate of each pair of its fixed matching. Pair complementary
   triples, and reuse the resulting four colors in both six-blocks.
   Within each color the supports are disjoint. For each of the eight
   signs on a support, use one tetrahedron or its image under negation,
   according to the sign parity. The squared head amplitude is 3/8.
   Four tetrahedra are chosen with absolute cross products at most 1/6.
   This inserts 512 points and retires 32 equator rays.
2. **Nested dense tails.** The 128-word even cube syndrome code contains
   a 16-word linear subcode of minimum Hamming distance 6. It partitions
   each syndrome into eight subcosets. In two chosen syndromes, replace
   the old pair at each word by a tetrahedron at squared head amplitude
   2/5. Use a fresh, disjoint tetrahedron in each subcoset. Between
   distinct upgraded syndromes, impose absolute head product at most
   1/6 only on subcosets which have distance-2 cross pairs. For the two
   displayed syndromes that condition forms two disjoint complete 4-by-4
   bipartite graphs. Also impose this product bound against the old
   orthogonal heads of the other syndromes. Each upgraded syndrome
   adds 256 net tail points before the retirement of 32 equator rays.

The verifier checks the literal norm-six head rows, orthogonality,
ownership against the pinned old signed classes, all necessary discrete
tail relations, and the angular bounds using rational comparisons.
Sparse and dense directions have absolute tail products at most 1/2.
Distinct head lines also have absolute product at most 1/2, so mixed
sparse/dense pairs are safe by Cauchy–Schwarz in the two block amplitudes.
There is no need to enumerate the entire 58-million-point code.

Run `python verify.py` with NumPy and assertions enabled. The parent geometry
and lattice minimum are inherited from the pinned certificates, not
assumed for any new external data.
