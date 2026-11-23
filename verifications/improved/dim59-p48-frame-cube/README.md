# 59D: a weighted Pless-frame cube

This construction gives **56,969,091** kissing points, 1,520 more than
the preceding [59D construction](../dim50-60-p48-global-owner-refill/).
The parent certificate supplies the size and validity of its unchanged stock.

Use 24 mutually orthogonal minimum-vector axes of the 48D Pless lattice.
Use the first 16 as head colors for the 1,024 even sign words on 11 sites.
The color of a word is the bitwise XOR of its negative site indices
(numbered 0 through 10). A cube word has tail proportional to its signs
times the weights (50,50,50,50,50,50,50,50,49,49,49), normalized to
squared length 3/4. Include both signs of its head, with head length 1/2.
There are 2,048 such points. For pairs of different words, the literal
weighted sign test checks the needed inequality both within a color and
between colors. The two points associated to the very same word have
inner product exactly 1/2.

An exact computation in the real quadratic field with sqrt(2) shows
that the old tail-direction dot products are bounded in absolute value
by 1417/2500. Together with the maximum normalized head contribution
sqrt(3)/8, this allows the complete old mixed stock to be retained. The
uniform bound is not used to charge deletions: the verifier identifies
the 576 actual conflicting old pure-tail poles, with their signs.

On the other eight axes add both signs for each of the 32 rational
weighted tails in the witness (64 points), deleting their 16 old
equatorial points. The verifier checks the supplement angles, the
cross pairs with the weighted cube and the exact interactions with
the retained tetrahedral additions. It makes conservative exact overlap
checks against 4,586,146 entries of the old signed owner superset.
All other old points are retained.

The verified accounting is 56,967,571 + 2,048 + 64 - 576 - 16 =
**56,969,091**.

Replay from this directory (Python 3.10+ and NumPy):

    OPENBLAS_NUM_THREADS=1 python3 verify.py

The verifier pins its inherited inputs by SHA-256. Its computations
are exact integer or quadratic-field comparisons; bounded integer dots
in the owner check may use binary64 with all possible intermediate
integer sums below 2^51.
