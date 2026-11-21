# Larger shifted dodecacode outer lists in dimensions 92–96

| Dimension | Pinned prior lower bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 92 | 30,151,281,256 | 30,156,524,136 | 5,242,880 |
| 93 | 44,038,752,920 | 44,067,064,472 | 28,311,552 |
| 94 | 74,103,655,116 | 74,131,966,668 | 28,311,552 |
| 95 | 110,610,877,508 | 110,639,189,060 | 28,311,552 |
| 96 | 209,198,254,912 | 209,274,800,960 | 76,546,048 |

These are verified lower bounds, not optimality or unrestricted priority claims. The
geometric construction and inherited dense, Reed–Solomon and root layers are
unchanged from the [preceding shifted dodecacode proof](https://github.com/alexlegeartis/KissingNumbers/tree/defb2d72324ecff3f578e066e2c8acab136c33ae/verifications/improved/dim92-96-shifted-dodecacode).
The new outer lists **replace** the former outer lists; the two augmentation
layers are not unioned. Root and dense cross-layer checks inherit the pinned
proof because support weight and coordinate magnitudes do not change.

The verified lists contain 342, 1,064 and 2,970 quaternary weight-12 patterns,
respectively with masks 15, 3 and 0. Each contributes 4,096 dodecacode axis words
times 256 sign words. Two distinct patterns of Hamming distance at least 12
clear the angular threshold directly. At distance eight the accepted pairs
share exactly two same-label tetrads, the third tetrads are disjoint, and the
axis cosets differ in at least two positions on the eight shared rows by the
pinned quotient construction. These special pairs number 31, 143 and 518,
respectively. Each layer is verified against all its actual unordered outer
pairs, not just an edge list or the result of a heuristic.

Run `node verify.js` from this directory (Node.js 18 or later; no external
modules). The verifier recomputes the binary Golay words, assigns the three
orthogonality-preserving phases, checks projected quotient distances and all
outer pairs, and prints the exact counts in `claims.json`. For a second
implementation run `python verify_shifted.py` (Python 3.10+); it builds the
projected quotient via exhaustive coset enumeration rather than Gaussian
elimination, checks the inherited body receipts, and repeats all pair checks.
The search optimizer is provenance only: its seed was the pinned prior bound,
and it used one-to-two and two-to-three exchanges in the outer conflict graph.

The independent replay uses the pinned inherited checks in `audit.py`. The `SHA256SUMS` file covers the core proof inputs.
