# Complementary six-row parity packing in dimensions 94–96

Append these points to the [paired-group BCH construction](../dim94-96-bch-paired-group-append/). No incumbent point is removed or moved.

| Dimension | Preceding bound | Gain | New lower bound |
|---:|---:|---:|---:|
| 94 | 74,137,471,692 | 1,245,184 | **74,138,716,876** |
| 95 | 117,087,145,028 | 1,245,184 | **117,088,390,212** |
| 96 | 209,285,286,720 | 1,794,048 | **209,287,080,768** |

A support selects six ordinary four-coordinate rows, with zero elsewhere, including the inherited head coordinates in 94D and 95D. It carries a row parity mask and the complement of that mask relative to the support. Each fixed row parity gives eight cube vertices, a signed orthogonal frame with four axes. Two rows of differing parity have absolute inner product two. Use the explicitly generated 128-word length-six quaternary axis code, obtained by puncturing the length-seven construction from the [91–93D parity pack](../dim91-93-flat-parity-extensions/), and the 32 even six-bit signs. There are 8192 points per support in both modes together, of squared norm 24. With the same mode, differing axis words have at most three coincident axes; equal axis words have at least two signs changed. With opposite modes, each row contributes at most two. In every case the raw dot is at most 12.

Across different selected supports, let their overlap have size k and let d of these rows have opposite parity in the stored modes. Allowing each supported parity mask or its complement gives an upper bound of 4k−2 min(d,k−d). The verifier checks it is at most 12 on every pair. It also reconstructs both inherited BCH support types in the correct ordinary coordinates and checks that each meets a new support in at most three rows. Therefore those cross products are at most 12 without any assumptions on the incumbent signs.

The other inherited layers have squared norm 96. The three dense body alphabets and the shifted Golay–dodecacode row patterns have absolute product at most four against a new occupied row. There are just six such rows, so the product is at most 24, as required for squared norms 24 and 96. The inherited Reed–Solomon stock has six individual entries of magnitude four and the root stock has two entries of magnitude sqrt(48), also giving the required cross bound. The new points have different support profiles from those layers even after normalization. The old flat supports are distinct from the new ones by the overlap bound.

Run `python3 verify.py` (Python 3.10 or newer). It independently checks the finite axis and parity geometry, the support patterns and every support pair, and the exact counts. The inherited construction and its previously checked geometry remain separate premises. These are ordinary kissing-number lower bounds, not optimality claims.
