# Conditional index-eight extension in 96D

The published DGLS construction has an even-parity Gray layer C of 2^36
words and an odd-parity companion v+C of the same size. Let H be the index-eight
subcode obtained by requiring the first three coefficients in a free Z8 basis
to be even. The triangular generator proves |H|=4^3*8^9=2^33, and C is the
union of the eight indicated H-cosets.

A separately chosen metric isometry of the eight even four-bit words in each
of 24 blocks maps H to a code T(H) with the same cardinality and internal
minimum distance as H. Its words are at distance at least 24 from every word
of v+C, because parity differs in every block. If and only if an explicit
choice also satisfies distance at least 24 from C, the union C ∪ (v+C) ∪ T(H)
would contain 2^37+2^33 = 146,028,888,064 words. Disjointness would follow
from that cross-distance. This is conditional, not a bound or a claimed find.
The current repository comparison used here is 142,162,362,560 in 96D.

`run.py` uses the general block-isometry discovery engine to search all 384
isometries independently in each block, sampling from eight reference cosets of C. The
large cross-distance oracle is heuristic and does not exhaust any coset. Every counterexample it finds is an
exact rejection witness, but absence of one does not prove compatibility.
Any candidate needs an independent exhaustive or mathematically sufficient
cross-distance certificate, plus current literature comparison, before a claim.

    OPENBLAS_NUM_THREADS=1 python research/dgls96_index8_20261105/run.py --seconds 60 --rounds 8
