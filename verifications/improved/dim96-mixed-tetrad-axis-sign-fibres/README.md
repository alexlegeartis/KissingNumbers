# Mixed three-tetrad axis and sign fibres in dimension 96

The exact lower bound is **209,299,303,328**, up **4,054,112** from
209,295,249,216. All previously certified points are retained. The new
layer comprises 27 additional modes, using 16,156 distinct
(mode, axis word) pairs with nonempty sign fibres.

A mode is specified by two binary Golay octads whose intersection has
size four. Their three disjoint tetrads, (a\setminus b), (b\setminus a),
and (a\cap b), use three distinct D4 matching frames, in that order.
Each of the twelve rows has two entries of magnitude 2, hence squared
norm 8, and each row's unsigned projective axis is a symbol in a
four-letter alphabet. Start from the 4,096-word dodecacode of length 12
and minimum quaternary distance six. The witness gives a coordinate
permutation and an XOR shift separately for each mode. Thus distinct
axes within one mode leave at most six aligned rows. On any one axis, an independently chosen sign code of minimum
binary Hamming distance three supplies the signs.

## Exact old interfaces

Across two modes sharing q occupied rows, of which r use the same
matching frame, the product is at most (4(q+r)) before using axis or
sign information. Every active pair with the literal 2,970 old
shifted patterns has q+r=16. The old quaternary shifts are regenerated
from the old Golay tetrads and the published projected-quotient phase
rule. The 1,615 newer monochromatic modes and every new-new pair are
already safe by the raw geometric cap.

For an active interface, express row products in units of four.
An equal-axis, equal-frame row contributes ±2; a differing axis in the
same frame contributes zero; a different frame always contributes ±1.
If the old projected axis differs in two or more matched slots,
the mixed product is at most 12 units. Otherwise its raw absolute
sum is 16 or 14 units. To exceed 12, the signs must then align everywhere
with a nonzero row product, except that the 16-unit case also allows
one mismatch in a different-frame row. The verifier enumerates exactly
these possibilities. Projecting the old shortened Hamming code is
equivalent to asking whether the prescribed syndrome lies in the
span of the erased one-based coordinate labels. For each new axis
word, remove every sign word that could exceed the cap. Two exact sign-code recipes compete. The first takes the largest remaining
subset of the 16 shortened Hamming cosets. The second checks which one-
and two-bit changes can actually occur within the remaining sign set.
Make a graph on the twelve coordinate labels and an extra zero sentinel.
Join two coordinates when a sign pair can differ in those two bits,
or join a coordinate to the sentinel when it can differ by itself.
Give the sentinel color zero and color all remaining vertices
deterministically by saturation degree, then degree, then slot order.
A single XOR-syndrome coset of the resulting labels cannot contain any
of the one- or two-bit differences. Retain the largest such coset.
The replay uses whichever recipe gives the larger set, with first
cosets winning ties within either recipe and the Hamming recipe winning
an inter-recipe tie. This coding step is per axis word.

The remaining interfaces with the dense, grouped-root, amplitude-four
Reed–Solomon and six-row triality/flat layers do not use the axis
phases and follow from the same per-row bounds as the shifted stock.

## Replay

Run `python3 verify.py` with assertions enabled and g++ available.
The wrapper builds the C++17 integer replay in a temporary directory.
There are no external Python packages or architecture-specific compiler
flags. Every line of `witness.txt` specifies the two decimal octad masks,
the permutation of slot labels 0–9ab, a decimal 24-bit axis XOR shift, and
the exact point count. Within each tetrad, rows are in increasing physical
order, and the twelve packed axis symbols use two bits per slot. The
search for good permutations is not used in replay.
