# Conditional extensions by local Gray-block isometries

This is a research direction, **not a new kissing-number bound**. Let C be the
2^36-word binary Gray image of the Z8 Golay code in the DGLS construction,
and let v have one set bit in each four-bit block. The known code C ∪ (v+C)
has minimum binary distance 24. The eight even-parity four-bit words form four
antipodal pairs. An arbitrary permutation of the four pairs and an independent
flip of either member in each pair gives 4!·2^4 = 384 isometries of this
8-point metric space. Apply one such map independently in each block to C,
obtaining T(C). Each T(C) has the same internal distances and cardinality as C,
and every T(C) word has distance at least 24 from every word of v+C, because
their parity differs in all 24 blocks.

The unproved condition is d(T(C), C) ≥ 24. If an explicit T satisfying this
condition is found and independently certified, then C ∪ (v+C) ∪ T(C) has
3·2^36 = 206,158,430,208 words and gives that many 96-dimensional kissing
points. The research search finds counterexamples for individual T; failure
to find one is not proof of the condition.

Moreover, the first m ≤ 12 ring coordinates of C are jointly uniform over
Z8^m, because the first 12 columns of the generator form an invertible
triangular matrix over Z8. Each block isometry permutes the even words: one
has block sum +4, one has sum −4, and six have sum 0 after sign embedding.
Thus the same uniform section counts apply to T(C). Conditional on the
unproved compatibility, the three dense layers *alone* give the following
counts, without assuming compatibility with any other accessory layer:

| Dimension | Conditional count | Existing repository count |
|---:|---:|---:|
| 86 | 7,739,670,528 | 6,225,165,512 |
| 87 | 10,319,560,704 | 6,226,194,610 |
| 88 | 13,759,414,272 | 6,921,251,584 |
| 89 | 18,345,885,696 | 9,229,354,136 |
| 90 | 24,461,180,928 | 12,306,670,328 |
| 91 | 32,614,907,904 | 16,409,599,116 |
| 92 | 43,486,543,872 | 21,880,490,288 |
| 93 | 52,143,587,328 | 31,040,958,640 |
| 94 | 75,698,798,592 | 47,697,508,076 |
| 95 | 115,964,116,992 | 77,565,848,332 |
| 96 | 206,158,430,208 | 142,162,362,560 |

For dimensions 86–92, impose zero block sums on the first m=96−d blocks.
For 93–94, impose equality of the first m=97−d block sums. For 95, impose
that the first two block sums add to zero. In each case the constraints are
independent, so the remaining ambient dimension is d. Exact count arithmetic
is reproducible with `python derive_projection_counts.py`. The cited 86–87 repository counts come from `RESULTS.md`; the 88–96 counts
come from `research/discovery_campaign_20261104/full96/replay.log`.
