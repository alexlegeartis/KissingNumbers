# Signed support augmentation in dimensions 64 and 65

Checkpoint: 24 October 2026.

| Dimension | Constructed lower bound | Gain over ERS64, plus two poles in 65D |
|---|---:|---:|
| 64 | **332,170,112** | 432,128 |
| 65 | **332,238,338** | 500,352 |

These explicit lower bounds improve the 331,737,984-point construction in
[Edel, Rains and Sloane](https://arxiv.org/abs/math/0207291).
The historical 65D comparison is its direct sum with two poles.
This package does not establish unrestricted optimality or certify global priority.

## Replay

```bash
python final.py
python verify.py
```

Python 3.9+, NumPy and g++ are required. The verifier compiles the supplied
enumerator in a temporary directory and freshly enumerates both binary
[64,28,16] codes: 536,870,912 words altogether. It reconstructs all owner
codes and designs, checks every added support against all 30,828 original
supports and every other added support, and checks the sign certificates.
No files outside the repository or network access are needed.

## Dimension 64

Every point has squared norm 64 and distinct points have inner product at
most 32. The original layers contain 2^28 dense sign vectors, 30,828
weight-16 supports carrying 2,048 extended Hamming sign patterns apiece,
10,416 affine four-point supports carrying all 16 signs, and 128 axial points.
The amplitudes in those four layers are respectively 1, 2, 4 and 8.

Add 211 further weight-16 supports, each carrying 2,048 signs. Any two new
supports intersect in at most eight coordinates. Each new support meets
exactly eight old supports in ten coordinates, and all other old supports
in at most eight. The 1,688 affected old supports are distinct.

For each ten-coordinate intersection, the certificate supplies a shared
eight-set that is a parity check of the new support's Hamming sign code.
Give the old support a coordinate order making this eight-set an affine
hyperplane of F2^4, and shift its sign code so this check has odd parity.
The new signs have even parity there. At least one sign disagrees, so the
cross inner product is at most 4(10-2)=32. No old points are lost. Every
other support pair has inner product at most 4*8=32 regardless of signs.
Dense, four-support and axial cross bounds remain 32, 32 and 16.

Thus the gain is 211*2,048=432,128.

## Dimension 65

Choose the dense code independently of the support code. It is the doubly
even radical of the extended binary BCH [64,30,14] code. Its dimension is
28, and its minimum distance is 16, checked by exhaustive enumeration.
It contains the extended BCH [64,16,24] code.

An antipodal owner set U consists of the 65,536 dense vectors from that
subcode and eight even-parity sign vectors on each of the 336 affine
F4-lines of F4^3. The lines belong to the affine four-support design.
Distinct non-antipodal owners have absolute inner product at most 16:
the dense subcode has weights 24 through 40 away from its antipodal pair;
distinct F4-lines meet in at most one coordinate; and a dense vector meets
a four-support vector with absolute inner product at most 16.

For each of the 68,224 owners u, replace (u,0) by
(sqrt(3)/2*u,+4) and (sqrt(3)/2*u,-4), then add the poles +/-8e65.
Cap versus retained equator products are at most 16sqrt(3)<32.
Distinct non-antipodal owner caps have product at most 28; the two caps
of one owner have product 32. Antipodal owners and the poles also meet
the bound. All points retain squared norm 64. This adds 68,224+2 points
to the new 64D construction.

## Files

- `final.py` derives the two totals from the number of supplied blocks.
- `verify.py` replays both enumerations, every support incidence and sign check.
- `construction.py` rebuilds the BCH codes, affine designs, owners and cap inequalities.
- `enumerate_binary.cpp` performs exhaustive Gray-code enumeration.
- `data/construction.json` specifies all generators, added supports and shared checks.
- `data/SHA256SUMS.txt` covers the complete package except itself.

The support family was found by a deterministic greedy search with seed
642410 among the weight-16 words of the original support code's dual.
Search optimality is not needed by the certificate and is not claimed.
The research scripts are in `research/dim64_65_72_20261024/` in the surrounding
workspace; replay of this package does not require them.
