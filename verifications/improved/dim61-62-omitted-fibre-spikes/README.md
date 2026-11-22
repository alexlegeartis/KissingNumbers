# Omitted-fibre weighted spikes in 61D and 62D

| Dimension | Previous bound | Inserted | Retired poles | New lower bound |
|---:|---:|---:|---:|---:|
| 61 | 89,082,322 | 33,792 | 66 | **89,116,048** |
| 62 | 122,668,916 | 1,081,344 | 66 | **123,750,194** |

Keep the original cyclic dense code, together with the exact old dense
latitudes and all incumbent sparse layers. Its 60 ordinary coordinates
remain after deleting coordinates 24, 42, 53 and 55 from a 64-coordinate
binary code of minimum Hamming distance 16. The former dense points
occupy deleted labels 0,1,2,3 in 61D and 0,1,3,7,6,4 in 62D.
Use the omitted label 6 in 61D, and omitted label 2 in 62D.

Pick 33 ordinary axes whose marks are not used in the existing 35
weighted capsule classes. On a marked axis put the magnitude `a b` and
on the 59 other ordinary axes put magnitude `b`, with the unaltered signs
of a binary word in the omitted fibre. In 61D take `a = 11/2` and
head `r = 2/5`; in 62D take `a = 9/2` and head zero. Scale the ordinary
part by `b² = (64 − r²)/(59 + a²)`.

In each dimension, an explicit code of distances at least 24 or at least
20, respectively, lies in the common shortened kernel. Select one of its
distinct affine cosets per marked coordinate to get disjoint sets of signs.
The result has pairwise gaps in one axial class by the verified Hamming
distance, and between classes because the marks are distinct and all
selected words have the same deleted label. There is no dense-owner
removal; only the two old poles per chosen ordinary axis are retired.

All old weight-sixteen sign classes have odd sign parity on supports
dual to the parent code. A weighted body copies an even sign restriction,
so at least one sign must disagree with every old sparse sign word on
its support. One such disagreement suffices even when the old support
contains the mark. The exact verification also checks the old
four- and nine-support layers, incumbent capsule classes, retained poles,
and interactions with the original 61D/62D dense latitudes.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy. The verifier
reconstructs the actual cyclic generator from a pinned sibling, checks
the literal shortened-code spans and disjoint cosets, and compares all
necessary dot-product envelopes with exact rational arithmetic. The
verification does not depend on the heuristic used to discover the
shortened subcodes.
