# Coded-layer search in dimensions 32 through 35

The 33D search increased the certified lower bound from 370,416 to **370,528**. Its independent-set calculation used both signs of individual octad points, signed stars, and the roots they conflict with. The certificate contains the exact point data, an expanded pair check, and a frozen search instance. [33D construction](../../verifications/improved/dim33-joint-signed-refill/)

| Dimension | Starting bound | Best certified count in this search |
| ---: | ---: | ---: |
| 32 | 352,992 | 352,992 |
| 33 | 370,416 | **370,528** |
| 34 | 391,004 | 391,004 |
| 35 | 415,360 | 415,360 |

## When reweighting a larger binary code cannot work

If coordinate $i$ of every sign vector has common magnitude $a_i$, the squared norm is $L=\sum_i a_i^2$. A binary difference with support $B$ gives inner product $L-2\sum_{i\in B}a_i^2$. The kissing condition requires $\sum_{i\in B}a_i^2\ge L/4$ for every nonzero difference.

For each of the specified $[34,18,8]$ and $[35,19,8]$ shortenings in CodeTables, there is a weighted collection of weight-eight words that covers every coordinate equally. Under the probability distribution proportional to the integer multiplicities in the witness, the mean of $\sum_{i\in B}a_i^2$ is exactly $8L/n$. It is below $L/4$ in both dimensions, so at least one pair violates the kissing angle for every common magnitude assignment with $L>0$. The minimum weighted distance is maximal at equal magnitudes within each of these two codes. [CodeTables generator for the parent](https://www.codetables.de/BKLC/BKLC.php?q=2&n=34&k=18) · [Integer witnesses](weighted_sign_witness.json)

Run `python verify_weighted_signs.py` to enumerate the two shortenings from the parent and check the balanced collections using integer arithmetic. The obstruction concerns the specified binary word sets under common coordinate weights.

## Other finite searches

- **Projective orbits.** The `PSL(2,11)` action on three projective lines recovers the 2,739-octad structured part of the 36-coordinate construction. Replacing a projective-line action by the degree-11 action gives 1,584 octads in the tested 35-coordinate orbit packing. The tested group actions on 32 through 34 coordinates also fell below the incumbents. [36-coordinate source](../../verifications/improved/dim37-38-qr-merge/supports37.py)
- **Split tail signs.** On a fixed set of marked coordinates, let signed octad tails be $t$ and $u$. The maximum body overlap allowed for complete body sign families is $\lfloor 4-|t\cdot u|/4\rfloor$. Assigning separate body families to orthogonal tail lines allows a larger overlap than requiring the same body support on both lines. The tested marked-coordinate sets and rational tail rotations produced no larger witnesses.
- **Whole-octad exchanges.** Low-deletion pools around the published support lists contained 19,033, 61,541, 114,948 and 185,052 nonincumbent octads with up to six blockers in 32D through 35D, respectively. Searches that opened one coordinate, several coordinates, or parts of a local support family retained the published support counts of 1,713, 1,852, 2,013 and 2,202. These are restricted search results. [Published support counts](https://aeb.win.tue.nl/codes/Andw.html)
