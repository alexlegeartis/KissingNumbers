# Exact four-site alphabet cap in 65 dimensions

The [explicit grouped four-site construction](../dim65-66-grouped-four-site-designs/) in dimension 65 has 10,726 four-site supports and 419,000,034 kissing points. Replace just that signed support layer. The optimal quadruple-packing theorem of [Bao and Ji, Theorem 4.9](https://arxiv.org/html/1401.2022#S4), gives a family of 10,741 four-subsets on 65 sites with no common triple, since the Johnson bound here is
`floor(65/4 * floor(64/3 * floor(63/2))) = 10741`.
Take all sixteen sign patterns on every support. They are mutually compatible and remain compatible with every unchanged vector by the same magnitude-only bounds as the explicit construction. This gives **419,000,274**, an additional 240 points.

The signed support layer also attains its exact ceiling. Given any signed selection, choose a uniform global sign pattern and retain exactly the vectors that agree with it on their support. The surviving supports form a quadruple packing, and each vector survives with probability one sixteenth. Thus at most `16 * 10741` points may lie in this alphabet, not a global upper bound for the kissing number.

The upgrade is an exact **existence proof using the cited design theorem**, not a file of 10,741 literal supports. The predecessor remains the fully explicit fallback. Run `python3 verify.py` for the predecessor and exact composition arithmetic, or use `--composition-only` to check the composition alone.
