# Punctured cubic Reed–Muller codes in dimensions 60 and 61
 
The verified final bounds after [simultaneous support rephasing](GUARDED.md) are **68,033,576 in 60D** and **118,639,402 in 61D**. Run `python3 verify_extra.py` for the final proof. The table below documents the base construction before additional guards.

| Dimension | Previous repository bound | Exact new lower bound | Gain |
|---:|---:|---:|---:|
| 60 | 58,665,998 | **66,657,320** | 7,991,322 (+13.62%) |
| 61 | 91,038,762 | **116,711,210** | 25,672,448 (+28.20%) |

These configurations replace their predecessors; they are not unions with the old P48 or cyclic-code configurations. Both reuse the exact binary parent of the [published self-sparse family](../dim62-67-rm2-selfsparse/): 80 distinct cubic cosets of RM(2,6) in 64 coordinates, minimum Hamming distance at least 16. All parts of the construction use ordinary Euclidean metric and then normalize each vector by its own length.

## Dense signs

Delete the four bit sites (0,20,41,61), in this order for the four-bit integer label. Their evaluation map on RM(2,6) has rank four, so each label has 80 × 2^18 words. In 60D take labels 0 and 1 and use ±1 on all 60 coordinates. Parent distance 16 becomes at least 15, giving the required inner product bound 30.

In 61D use the same 60-coordinate body and append the following height:

| label | 0 | 1 | 2 | 3 |
|---:|---:|---:|---:|---:|
| height | 1/2 | 2 | −1/2 | −2 |

For parent labels j,k, let d be their four-bit Hamming distance. The body dot is at most 28+2d, and the squared norms including the heights h,t are 60+h² and 60+t². All ten unordered label pairs satisfy
`28+2d+ht ≤ sqrt((60+h²)(60+t²))/2`.
The unequal magnitudes are essential: at adjacent same-sign or diagonal opposite-sign labels, the worst dot 31 compares with half the product of radii sqrt(964), which is strictly larger than 31.

## Sparse signs

Every weight-16 parent word is orthogonal modulo two to every dense parent word. Products involving quadratics have degree at most five and hence even evaluation parity on F2^6; the verifier directly checks the cubic representatives are pairwise orthogonal as well. Distinct weight-16 parent supports meet in at most eight sites by parent distance 16.

In both dimensions take the **12,005** weight-16 supports avoiding all four marks. Put ±1 on each support, with the **odd** coset of a binary [16,11,4] extended-Hamming code: its parity-check columns are (1,a), a∈F2^4; the syndrome is zero in the four a rows and one in the overall parity row. This gives 2,048 signs each. A dense word has even sign parity on the support and therefore disagrees somewhere, giving mixed dot at most 14, strictly below the mixed limit for squared norms 16 and ≥60.

In 61D also take the 3,956 supports hitting only mark 0, and the 3,956 hitting only mark 20. Put ±1 on the 15 ordinary sites and use height +1 for the first class and −1 for the second. On the 15 ordinary sites use the shortened extended-Hamming coset with overall parity opposite that of the **unique** threatening dense label on the removed mark. This retains 1,024 signs for each support. If the height faces the dense height ±2, its unguarded bound 17 exceeds the threshold 16, but the parity condition forces a mismatch and drops the dot to at most 15. For the small facing height, the bound is 31/2 ≤ sqrt(241). All other dense/head pairs are angularly safe without parity. On a common mark, two distinct supports have ≤7 shared ordinary sites and the head adds 1. On different marks, the ordinary overlap is ≤8 and the opposite heads subtract 1. The mixed dot against an unmarked support is ≤8. Within a fixed signed bank the sign distance is at least four.

## Affine planes and coordinate poles

Add all 16 sign choices on every affine four-set contained in the physical coordinate names, and both pure poles in each coordinate. The actual fields are all 64 original coordinate names except the four marks in 60D (7,995 affine four-sets); in 61D only three names are missing, because name 0 is repurposed as the new height axis (8,555 affine four-sets). Affine planes meet in at most two positions. Against dense vectors the maximum raw product is four for a plane not using the tail, or 3+|h|≤5 for a plane using it; the mixed threshold is at least sqrt(60). Against sparse vectors all nonzero coordinates have magnitude one, giving at most four, the mixed threshold. Against a pole every family has normalized absolute coordinate ≤1/2. These checks also exclude duplicate normalized points.

The exact decompositions are
`41,943,040 + 12,005×2,048 + 7,995×16 + 120 = 66,657,320`
and
`83,886,080 + 12,005×2,048 + 7,912×1,024 + 8,555×16 + 122 = 116,711,210`.

## Replay

Run `python3 verify.py PATH_TO_EXISTING_DIM62_67_RM2_SELFSPARSE_PACKAGE` with Python and a C++17 compiler. The parent source files and premise are hash-pinned; the parent independently regenerates all 37,996 weight-16 supports twice. Then the new verifier checks prefix rank, sign parity, all dense height and mixed gates, shortening, both exact counts and the explicit labels. Neither optimizer nor floating point is needed. The repository layout uses the pinned sibling package by default.
