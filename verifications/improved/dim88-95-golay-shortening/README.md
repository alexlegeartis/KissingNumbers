# Shortened Golay–dodecacode augmentations in 88D, 89D and 92–95D

| Dimension | Preceding lower bound | Increase | New lower bound |
|---:|---:|---:|---:|
| 88 | 12,884,901,888 | 92,749,824 | **12,977,651,712** |
| 89 | 13,958,643,714 | 92,749,824 | **14,051,393,538** |
| 92 | 30,121,921,128 | 19,922,944 | **30,141,844,072** |
| 93 | 43,934,943,896 | 53,477,376 | **43,988,421,272** |
| 94 | 73,999,846,092 | 53,477,376 | **74,053,323,468** |
| 95 | 110,507,068,484 | 53,477,376 | **110,560,545,860** |

These are lower bounds from explicit spherical codes, not upper bounds or optimality claims. The comparisons are with the preceding repository certificates.

## 88D and 89D

Start with 629 length-24, weight-12 quaternary words from pairs of binary Golay words. They have Hamming distance at least 12, vanish in original blocks 0 and 1 (zero-based), and are all active in block 3. Set the block-3 symbol to zero. Each projected pattern has weight 11, and distinct patterns retain distance at least 11. Thus they remain distinct and have at most 11 total double-occupancy-plus-equal-tag incidences. Using the [three perfect matchings, scaled axes and pattern bound](https://github.com/alexlegeartis/KissingNumbers/tree/27247d2c57e82d589f688b2996ec2dfabf8ac4fa/verifications/improved/dim92-96-golay-dodecacode), each pair of points from distinct patterns has inner product at most 44.

Shorten the additive dodecacode by taking the 1,024 words with last symbol zero, then deleting it. This preserves distance six. Use 144 binary signs from the punctured Julin–Golay code of length 11 and distance three, reproduced as cyclic shifts and negations of quaternary seeds, Gray mapped and then punctured. See [the quaternary construction](https://neilsloane.com/doc/Me185.pdf). The verifier checks all sign pairs, not an assumed minimum-distance theorem. For one pattern, different axis words can match in at most five rows, contributing at most 40. Equal axis words with distinct signs have at least three sign disagreements, also contributing at most 40. The result is 629 × 1,024 × 144 = 92,749,824 new norm-88 points.

The [88D and 89D parents](https://github.com/alexlegeartis/KissingNumbers/tree/27247d2c57e82d589f688b2996ec2dfabf8ac4fa/verifications/improved/dim88-89-z8-prefix-anticode) keep all 22 body blocks after deleting exactly blocks 0 and 1, each populated by one of the three incumbent dense alphabets. Each of eleven new nonzero blocks has maximum dot 4 against any old row, giving raw dot at most 44. In 89D, the new last coordinate is zero. The eleven zero-height parent fibres have squared norm 88, and the two raised parent fibres have squared norm 272/3. In both cases, a raw dot of at most 44 gives normalized dot at most 1/2. The two pure poles are orthogonal. The new body has just eleven nonzero rows, while every old dense body has all 22, so none of the points duplicates an inherited point.

## Shortened outer-code improvements in 92–95D

Keep the geometry and inner codes of the [previous augmentation](https://github.com/alexlegeartis/KissingNumbers/tree/27247d2c57e82d589f688b2996ec2dfabf8ac4fa/verifications/improved/dim92-96-golay-dodecacode) but replace its outer pattern lists with larger ones. In 92D the zero mask 15 list has 328 patterns, up from 309. In 93–95D the zero mask 3 list has 989 patterns, up from 938. The inner layer size stays 2²⁰. The old pre-augmentation codes are retained; the new lists replace the previous appended layers, and are not asserted compatible with those previous layers. The pinned norm-96 geometry compares each new occupied row against the complete old dense row alphabet. The same sparse Reed–Solomon and root estimates apply unchanged.

Run `node verify.js`. It checks all Golay-pair lists, also after puncturing, all exact small-code distances, and the resulting counts. A large numerical optimizer is not needed to replay the lists or the proof.
