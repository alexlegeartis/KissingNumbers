# 38D: 652,028 by optimizing the affine offsets after compression

The earlier [QR merge](../dim37-38-qr-merge/) gave 643,836. This
continuation adds **8,192**, reaching **652,028** without changing
the signed octads or roots.

The original 39-bit kernel of Yang et al. is restricted to words with
the first two bits zero, then shortened. Its 37-bit shortening has
dimension 13, hence 8,192 words per coset. On the 37 ordinary
coordinates each word has entries ±2. One extra coordinate is fixed
separately for each fibre at +3, 0 or −3. The 8,192 points of a fibre all
use the same scalar. The first two old coordinates do not appear in
the output.

The old merge's fibres came from the published full-length parent.
Here the offsets are chosen *after* shortening, independently in each
scalar tier, using a 24-bit syndrome quotient. Dot products specify
the exact Hamming constraints between fibres: at least 10 within a
tier, at least 9 between zero and either nonzero tier, and at least
8 between tiers of opposite sign. Enumerating possible new offsets
and solving the compatibility graph gives 11 positive, 10 zero, and
10 negative fibres, instead of 10 each. This is an existence proof,
not a claim of optimality over all such offsets.

| Layer | Count |
|---|---:|
| Shortened QR fibres | 253,952 |
| Signed octads on 3,088 supports | 395,264 |
| Signed roots | 2,812 |
| **Total** | **652,028** |

The unchanged signed octad supports are Qihong Lin's CC BY 4.0 list
already included in the [original package](../dim37-38-qr-merge/THIRD_PARTY.md),
as are the length-39 generator masks. The new pair geometry uses the
same integer scale as before. The explicit syndrome list is
`offset_witness.txt`. Both checks enumerate the relevant coset
differences, and prove the middle/root and cross-layer angle bounds.

Run from this directory:

```sh
python verify.py               # standard library
python audit.py                # vectorized difference checks (NumPy)
python audit.py --full-dense   # all explicit dense–dense pairs
```

The optional audit expands every dense point and checks every unordered
dense–dense pair. Other pair types are covered by the exact layer
bounds, not an explicit quadratic cross-layer enumeration.
