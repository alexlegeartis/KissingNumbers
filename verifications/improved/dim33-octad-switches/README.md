# 33D: two full-octad swaps and two signed star families

This gives **370,416** kissing directions, eight above the
[370,408-point predecessor](../dim32-35-signed-layers/nonantipodal33/)
and 176 above [Lin's 370,240](https://arxiv.org/abs/2611.03682).

The [recipe](witness.json) changes two zero-based octad supports;
each takes all 128 even-parity signings. The two star templates use
both center signs of magnitude 4 and all eight even-negative
signings of magnitude 2 on each leaf quartet. Of all 2,112 signed roots,
180 conflict with a star, 24 more than before. The dense code, the
other complete octads, all 2,064 individual exceptional octad
points and all 180 old stars are unchanged.

| Layer | Points |
| --- | ---: |
| Dense affine binary code | 131,072 |
| Full signed octads | 235,136 |
| Individual exceptional octads | 2,064 |
| Retained signed roots | 1,932 |
| Signed stars | 212 |
| **Total** | **370,416** |

Run `python verify.py` (Python 3.10+, standard library), or
`python audit.py` (also needs NumPy). Both rebuild the point set from
the pinned source and run independent exact proofs. The latter checks
86,520,912 explicit signed integer products and verifies full rank.
The code retains 500 unpaired points, without antipodal completion.
This is a lower bound, not an optimality claim.
