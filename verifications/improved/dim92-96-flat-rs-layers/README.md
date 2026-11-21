# Flat Reed–Solomon sign layers in dimensions 92–96

| Dim. | Inherited lower bound | Appended points | Combined lower bound |
|---:|---:|---:|---:|
| 92 | 30,156,524,136 | 1,327,104 | 30,157,851,240 |
| 93 | 44,067,064,472 | 2,097,152 | 44,069,161,624 |
| 94 | 74,131,966,668 | 2,097,152 | 74,134,063,820 |
| 95 | 110,639,189,060 | 2,097,152 | 110,641,286,212 |
| 96 | 209,274,800,960 | 4,194,304 | 209,278,995,264 |

This is a direct append to the
[preceding shifted dodecacode outer extension](../dim92-96-shifted-dodecacode-outer-extension/).
Its dense, sparse Reed–Solomon, root and shifted-pattern layers are all
retained. It does not require the apparently related 96D-only draft.

In 96D, divide the 24 original four-coordinate blocks into six consecutive
groups of four blocks. In ranks 93–95 the first two blocks are replaced by the
inherited separate head; leave it zero and use the same groups, with group 0
reduced to two blocks. In 92D the first four blocks are replaced by its
own head. Leave it zero and partition the 20 available blocks into six
consecutive groups of sizes 4,4,3,3,3,3. These group conventions affect only
which symbols are permitted.

Let GF(16) use the polynomial x^4+x+1. Fix either the four projective even-parity
flat axes in a four-coordinate block (the sign vectors with first coordinate
+1 and even parity), or the odd-parity axes obtained by negating first
coordinates. A GF(16) letter selects one block in its group and one axis out
of four there. Use the length-six RS axis words, evaluating polynomials of
degree at most 3 at elements 0,1,2,3,4,5. In a shorter group of b blocks restrict
the field letter to 0,...,4b−1. Give each retained axis word all 32 sign words of
even parity.

Every new point has squared norm 24. Distinct axis letters in a given parity
family have dot zero. Distinct axis words have at most three matching letters,
hence dot at most 12. An identical axis word with differing sign words has at
least two sign flips and dot at most 8. The two parity families have cross dot
at most 2 in each group, hence at most 12. The block support and flat signs are
recoverable, and their parity distinguishes the two families. The one,
four, or zero shorter symbol restrictions leave respectively the counts
shown by the verifier: in ranks 93–95 the constant term has eight choices,
in rank 92 the last four evaluations form a complete information set with
twelve choices each, and rank 96 has all 16^4 cubic words.

Only six source blocks are active, and cross dot with any incumbent dense
alphabet or shifted-dodecacode row is at most 4 in each, so total at most 24.
The incumbent sparse RS points have six entries with magnitude 4; root
points have two entries with magnitude sqrt(48). Since every new coordinate
has magnitude at most 1, both cross dots are bounded by 24, which is the
correct threshold sqrt(24*96)/2. The new points are disjoint from the dense
layers (whose ordinary blocks are all occupied), the shifted layers (with
two entries in each of twelve blocks rather than four in each of six),
and the other sparse layers (support size 2 or 6).

Run `python verify_flat.py` (standard-library Python 3.10+). It checks field
arithmetic, all allowed axis words, the sign distances, and the small
row-alphabet dot tables directly. The old code and its proof are inherited.
