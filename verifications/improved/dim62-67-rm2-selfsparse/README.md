# One binary code supplies both the dense and sparse layers in 62–67D

The explicit 80 cubic cosets of RM(2,6) from the
[previous package](../dim62-67-rm2-cubic-cosets/) give a binary code C of
size 335,544,320 and minimum Hamming distance 16. Here the *same* code
also supplies the sparse layer. No old 16-support list is needed.

| Dimension | Previous certified | New lower bound | Increase |
|---:|---:|---:|---:|
| 62 | 147,704,844 | **150,561,804** | 2,856,960 |
| 63 | 207,510,788 | **214,085,892** | 6,575,104 |
| 64 | 402,819,968 | **413,526,912** | 10,706,944 |
| 65 | 402,819,970 | **413,526,914** | 10,706,944 |
| 66 | 409,380,228 | **421,522,308** | 12,142,080 |
| 67 | 422,204,310 | **437,513,094** | 15,308,784 |

## 64D: minimum words become supports

Use squared norm 64 and inner product limit 32. All words of C supply
dense ±1 signs. Extract the **37,996** words of weight 16 from C.
For each support put ±2 on its sixteen sites, with signs from one copy
of the extended Hamming [16,11,4] code. An explicit parity check is
given by the sixteen columns (1,a) with a in F2^4. It has rank 5, and
no nonzero kernel word has fewer than four ones, giving 2,048 signs
and inner products at most 32 within each support.

Distinct supports meet in at most eight sites because their binary
distance in C is at least 16. They therefore coexist regardless of
signs. Every such sparse vector has L1 norm 32 and is compatible with
any dense sign. Finally take all signings of the amplitude-4 affine
4-point planes in F2^6 and the two amplitude-8 axial poles per coordinate.
Distinct affine 4-sets meet in at most two sites; their signs, the
weight-16 stock, the dense signs, and the poles satisfy the same
limit by the indicated amplitudes. The total is
335,544,320 + 37,996×2,048 + 10,416×16 + 128.

## 62D and 63D

In 62D mark binary positions 0,20,41,61 and fix the last to zero.
The other three bits take the six labels 0,1,3,7,6,4 in cyclic order.
Replace their signs with the regular hexagon
(2,0),(1,√3),(−1,√3),(−2,0),(−1,−√3),(1,−√3).
Every cubic coset projects at rank four on the marked sites, so each
tag carries 80×2^18 body words. The separation for each tag pair follows
from minimum distance 16 in C and the two head coordinates, with all
36 angular inequalities checked exactly. Keep the 12,005 sparse supports
avoiding all four marks. The four-support stock consists of the affine
planes avoiding physical coordinate names 41 and 61; the two remaining
marked names are repurposed as the two head axes. Their mutual and
sparse interactions use amplitude bounds, and the verifier checks
their dense-head interactions and the poles. This gives 9,145 four-sets
and 124 poles.

In 63D intersect the entire 64D code with the linear hyperplane
x_0=x_20. Exactly half the dense words survive. Among the weight-16
supports, 21,361 avoid both sites, and 2,363 contain both; discard the
rest. On a support containing both, the equality selects half its
signs (the dual has no weight-two word). The corresponding 4-set
counts are 9,145 and 31, plus 124 poles. The section inherits every
angle exactly.

## 65–67D: tag the fourteen shells

For 65D simply embed and add two poles. For 66D and 67D start with the
64D code and the extra poles. Assign the three cubic shifts
822310, 213640, 779738 to new two-axis patterns 011, 101, 110
(one pattern in 66D, all three in 67D). For each pattern, take the
3,904 words of weight 14 in C translated by that cubic. Append the two
indicated sites to each body support; on the resulting sixteen sites
use the same 2,048 extended-Hamming signs.

The minimum within each translated binary code remains 16. The minimum
between any two of these translated codes, or between one of them and
C, is at least 14, checked by independent exact bit enumeration.
For two distinct tail patterns the bodies meet in at most 7 and the
tails at most 1; against the zero-tail layer the body overlap is at
most 8. In a common tail pattern the bodies meet in at most 6. All
full 16-support intersections are therefore ≤8; an additional direct
integer support-pair loop independently verifies this. The upper sparse
cross product with the original dense signs is still 32.

## Replay

Run `python3 verify.py` with Python 3 and a C++17 compiler on the path.
The two dense distance enumerators use different encodings of RM(2,6);
the sparse stock is reconstructed by both Gray traversal and a Plotkin
split and checked for exact equality. This package needs no external
fixture or optimization software. `claims.json` records the six lower
bounds; no upper bound is claimed.
