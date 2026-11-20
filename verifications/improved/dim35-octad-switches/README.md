# 35D: seven octad switches and a signed-star exchange

This witness gives **415,358** kissing directions, 22 more than the
[415,336-point construction](../dim35-shared-root-heads/) and 50 more
than [Qihong Lin's 415,308](https://arxiv.org/html/2611.03682v1#S2).

| Layer | Points |
| --- | ---: |
| Dense affine binary code | 131,072 |
| 2,202 signed octads | 281,856 |
| Retained signed roots | 2,172 |
| Signed stars | 246 |
| Full-support heads | 12 |
| **Total** | **415,358** |

The [witness](witness.json) replaces seven octad supports, indexed from
zero in the historical `code_data.json`. The new supports take all
128 even-parity signings, and retained supports keep their certified
parities. It removes the eight positive-center stars at indices 62–69
and supplies 104 replacement or additional stars, independently signing
the centers as prescribed and using even parity on the four leaves.
The dense layer and all 12 heads are unchanged. Of the full set of 2,380
signed roots, the final stars and heads block 208. Relative to the
preceding construction, 80 roots become blocked and six are restored.
The net change is 96 stars for 74 roots, gaining 22.

Run `python verify.py` with Python 3.10+ and NumPy. The exact verifier
reconstructs the dense affine code, verifies all support intersections
and signed layer products or integer bounds, and reconstructs the
complete signed-root ledger. The construction is a lower bound, not an
optimality claim.

Run `python audit.py` for the independent coordinate replay. It checks
every octad's selected parity by expanding all 256 signings, uses signed
integer products throughout, and independently reconstructs the code
and root ledger.
