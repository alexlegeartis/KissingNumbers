# 40D: 1,139,120 kissing points from Smart's octads

Combine [B. Smart's \(A(40,8,8)\ge4267\) code (2026)](https://aeb.win.tue.nl/codes/cwc/d8/a40.8.8.4267H) with the length-40 dense sign layer used by [Edel, Rains and Sloane, *On Kissing Numbers in Dimensions 32 to 128*, 1998](https://arxiv.org/abs/math/0207291), and all \(D_{40}\) roots. The result improves the [previous repository certificate](../dim40-p48q-affine-refill/) by **12,755**. Smart is credited for the constant-weight code and Edel–Rains–Sloane for the layered construction and the published \(A(40,10)\ge589824\) input. This is a lower bound, not an optimum proof.

| Layer | Size | Unnormalized squared norm |
| --- | ---: | ---: |
| Dense binary signs | 589,824 | 40 |
| All even-parity signs on each of 4,267 weight-eight supports | 546,176 | 8 |
| Signed \(D_{40}\) roots | 3,120 | 2 |

Normalize each vector separately. The dense Hamming distance at least ten gives within-layer dot at most 20 before normalization. Distinct supports meet in at most four positions; on the same support, different even sign patterns differ in at least two signs, so the middle-layer unnormalized dot is at most four. Distinct roots have dot at most one. Across the three layer pairs (dense/middle, dense/root, middle/root), dots are at most 8, 2 and 2, respectively; their normalized values are at most one half. Different layers have different support sizes, so their cardinalities add without duplication.

The [Brouwer table](https://aeb.win.tue.nl/codes/Andw.html) indexes Smart's code. The replay verifies all **9,101,511** support pairs and the integer inequalities for all six geometric pair types. The 589,824 dense words are inherited from the cited published lower bound, not an independent dense-code enumeration. The source file is read from the original host, not mirrored here. Run either command from this directory:
```sh
python verify.py --download
python verify.py /path/to/a40.8.8.4267H
```
The replay pins the exact downloaded SHA-256. `python final.py` independently recomputes the registered count from its ledger, without network.
