# 34D: coupled octad swaps and signed stars

The construction gives **391,004** kissing directions in 34 dimensions:
12 more than the [390,992-point certificate](../dim32-35-signed-layers/rotated34/)
and 24 more than [Qihong Lin's 390,980](https://arxiv.org/abs/2611.03682).

| Layer | Points |
| --- | ---: |
| Dense affine binary code | 131,072 |
| Full even-parity signed octads | 257,664 |
| Signed roots | 2,152 |
| Signed stars | 112 |
| Rational heads | 4 |
| **Total** | **391,004** |

The [recipe](witness.json) preserves the dense code and four rational
heads, and replaces four of the 2,013 octad supports. Each star template
lists its center (both signs, magnitude 4) and its four leaves
(magnitude 2) followed by the prescribed parity of negative leaves.
Relative to the preceding certificate, 64 stars are added and 16
removed; 48 roots become blocked and 12 are restored. The final root
ledger removes 92 signed points on 46 antipodal lines.

Run `python verify.py` for the integer support/angle proof and
`python audit.py` for the independent expanded-coordinate proof, which
also checks uniqueness and exact rank 34. Both require NumPy and replay
from the compact recipe with no solver. All points have squared norm
5,780,000 at the common integer scale 425. The first pair of octad
edits and first six star templates independently give 391,000; to
replay it run `python -c 'from replay import verify; print(verify(intermediate=True))'`.
