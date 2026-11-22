# Omitted-fibre weighted spikes in 61D and 62D

| Dimension | Construction base | Inserted | Retired poles | New lower bound |
|---:|---:|---:|---:|---:|
| 61 | 89,072,802 | 1,966,080 | 120 | **91,038,762** |
| 62 | 122,659,076 | 1,966,080 | 120 | **124,625,036** |

These improve the preceding marked-capsule bounds by 1,956,440 and
1,956,120, respectively. The stock for these constructions is the one
*before* the marked-capsule exchanges, including the original append-only
nine-support layer. We trade the older capsule gains for much larger
weighted dense-fibre populations on all 60 ordinary axes.

Keep the original cyclic dense code, together with the exact old dense
latitudes and all incumbent sparse layers. Its 60 ordinary coordinates
remain after deleting coordinates 24, 42, 53 and 55 from a 64-coordinate
binary code of minimum Hamming distance 16. The former dense points
occupy deleted labels 0,1,2,3 in 61D and 0,1,3,7,6,4 in 62D.
Use the four omitted labels 4, 5, 6 and 7 in 61D, and omitted label 2 in 62D.

Use all 60 ordinary axes as marks. On a marked axis put the magnitude `a b` and
on the 59 other ordinary axes put magnitude `b`, with the unaltered signs
of a binary word in the omitted fibre. In 61D the profiles are
`a = 23/4` with heads `+3/2`, `−3/2` for labels 4, 5 and
`a = 11/2` with heads `+2/5`, `−2/5` for labels 6, 7. In 62D
take `a = 9/2` and head zero. Scale the ordinary
part by `b² = (64 − r²)/(59 + a²)`.

For 61D, relabel the original code by the explicit six-bit affine
coordinate list in the verifier. It checks that every quadratic
Reed–Muller row belongs to the actual original binary code. The
Delsarte–Goethals construction gives 131,072 words in this subspace
with minimum distance 24. Each of labels 4–7 selects 8,192 words.
Assign each marked coordinate one of 60 distinct prefix-zero shifts
modulo this quadratic subspace. The four exact angular profiles
account for the coordinate differences lost at each pair of deleted
labels and the product of their tails. The two new labels add
983,040 points without removing any of the previous points.

For 62D, an explicit rank-15 binary subcode of distance 20 in the common
shortened kernel supplies 32,768 points on each selected marked
coordinate. Choose one distinct coset per mark. The verifier generates
all the differences and confirms their distance. Different marks in
both dimensions use different words, so the inherited binary minimum
bounds their ordinary comparisons. No dense owner is retired; only
the two old poles of each marked coordinate are removed.

All old weight-sixteen sign classes have odd sign parity on supports
dual to the parent code. A weighted body copies an even sign restriction,
so at least one sign must disagree with every old sparse sign word on
its support. One such disagreement suffices even when the old support
contains the mark. The exact verification also checks the old
four- and nine-support layers, retained poles,
and interactions with the original 61D/62D dense latitudes.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy. The verifier
reconstructs the actual cyclic generator from a pinned sibling, checks
the literal shortened-code spans and disjoint cosets, and compares all
necessary dot-product envelopes with exact rational arithmetic. The
verification does not depend on the heuristic used to discover the
shortened subcodes.
