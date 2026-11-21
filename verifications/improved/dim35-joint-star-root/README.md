# 35D: exact joint signed-star/root repacking

This construction has **415,360** kissing directions, two more than the
[415,358-point repository construction](../dim35-octad-switches/) and 52 more than [Qihong Lin's
415,308](https://arxiv.org/html/2611.03682v1#S2).

| Layer | Directions |
| --- | ---: |
| Dense affine binary code | 131,072 |
| 2,202 full signed octads | 281,856 |
| Retained signed roots | 2,158 |
| Signed stars | 262 |
| Full-support heads | 12 |

The dense affine code, 2,202 octad supports and 12 full-support heads are
unchanged. The replacement NPZ stores all selected stars, roots and unchanged
heads as integer rows. The first four layers have raw squared norm 32;
the full-support heads have raw squared norm 51. A finite conflict graph
lets stars and signed roots enter or leave jointly, rather than requiring
entire fixed-parity groups of stars to move together.

Run `python verify.py` from this folder in the repository checkout
(Python 3 with NumPy, without `-O`). It checks source hashes, all explicit mixed-layer products,
compatible sign parity on each octad, and the entire binary difference
code using exact integer inequalities. No optimizer, previous
verification result or floating-point calculation is needed. This is a
lower bound, not a claim of optimality.
