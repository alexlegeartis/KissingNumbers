# Exterior extended BCH caps in dimensions 65–67

| Dimension | Previous | New lower bound | Gain |
|---:|---:|---:|---:|
| 65 | 335,929,116 | **336,256,796** | 327,680 |
| 66 | 342,381,366 | **342,389,558** | 8,192 |
| 67 | 355,183,103 | **355,195,391** | 12,288 |

Use squared norm 64 throughout. Let P be the parity-extended narrow-sense
BCH code from designed distance 13 at length 63. Over the binary field P
has dimension 30 and distance at least 14 (the BCH bound plus even parity).
Its doubly-even radical D has dimension 28 and is the dense equatorial
code in the preceding constructions. Distinct P bodies, written as sign
vectors, have raw dot at most 36.

## 65D: a distance-22 subcode meets all four radical cosets

The parity-extended BCH subcode E from designed distance 21 has dimension
18. Exact enumeration finds minimum distance 22 and exactly 65,536 words
in D. For every word u of E, insert (√3·u/2, ±4). Remove only the 65,536
equatorial copies in E∩D. A distinct equatorial dense body gives dot at
most 36√3/2 < 32. Same-sign cap pairs have dot at most
(3/4)(64−2·22)+16 = 31; opposite signs have dot at most 32.

Keep the prior [sparse line lifts](../dim65-dg-dense-owners/) and two poles.
Every old equatorial sparse point has L1 at most 32, so a dense cap
contracts that product. A dense cap meets a sparse line lift in body dot
at most 16 before the common scaling and tail, making the total at most 28.
The poles give at most 32. Relative to the verified 64D stock the count is
335,795,354 + 2·262,144 − 65,536 + 2,688 + 2 = **336,256,796**.

## 66D and 67D: exterior translates restore all equatorial owners

For all the prior 4,096-word [Kerdock-axis groups](../dim66-67-kerdock-axis-owners/),
XOR every body with v = 0x800000037cd0eb67, a word of P outside D.
This preserves body distances and the disjointness of the different
axis groups, and moves every cap body outside D. Retain the entire dense
equator. With the inherited scale s = √37/8, cross-dots with distinct dense
equator words are at most 36s < 32. Within-axis products are unchanged,
across-axis tails are orthogonal, and all sparse bounds are uniform over
ordinary sign choices. Each axis still retires two poles but restores its
4,096 previously removed equator words.

## Replay and scope

Run `python verify.py` from this directory with Python 3 and NumPy. The
predecessor certificates provide the verified sparse stock; the new
replay checks code nesting, distances and counts. Optional
`node coset_audit.cjs` reconstructs both BCH generators by a distinct
polynomial method and checks a weight-at-most-20 witness for all 4,095
nonzero P/E cosets. Therefore at the fixed 65D cap scale and sign,
adding any other P body to the full E cap conflicts with an E word.
This does not cover replacements, other shapes or bodies outside P.
