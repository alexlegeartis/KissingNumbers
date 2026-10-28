# Dimensions 88–96: nonlinear Z8 Golay codes and exact sections

This package corrects a missed published input and constructs stronger repository
bounds in nine dimensions. The nonlinear binary **(96, 2^37, 24)** parent is due
to **I. M. Duursma, M. Greferath, S. N. Litsyn and S. E. Schmidt (2001)**,
[IEEE Transactions on Information Theory 47(4), 1596–1598](https://doi.org/10.1109/18.923742).
It is not a new code discovered here. The former use of the linear [96,33,24]
code missed this stronger nonlinear input. The corrected published-code baselines
are recorded separately from the gains over the preceding repository checkpoint.

| dimension | verified lower bound | gain over preceding repository certificate |
|---|---:|---:|
| 88 | **6,921,251,584** | +689,012,056 |
| 89 | **9,229,354,136** | +2,993,829,100 |
| 90 | **12,306,670,328** | +6,064,643,188 |
| 91 | **16,409,599,116** | +10,157,396,566 |
| 92 | **21,880,490,288** | +15,608,058,713 |
| 93 | **31,040,958,640** | +24,739,322,556 |
| 94 | **47,697,508,076** | +41,339,850,680 |
| 95 | **77,565,848,332** | +70,464,307,200 |
| 96 | **142,162,362,560** | +128,849,018,880 |

The comparison checkpoint is commit `2b187e0`. No new bound in 73–87D is
claimed by this package. No unrestricted optimum or comprehensive worldwide
priority claim follows from this source review.

## Replay and certificate files

```sh
python verify.py
```

Run without `-O`. The wrapper checks the complete manifest, runs both proofs
read-only, compares every field with the pinned receipts, and runs malformed
certificate controls. It does not enumerate 2^37 binary words.

- `final.py` checks the arithmetic in `bounds.json` and supplies the registry rows.
- `source/replay.py` uses canonical Z4 cosets, scalar weight lookup, exact
  coordinate support products and signed-tail profile counts.
- `source/verify_independent.py` imports no proof code from the first replay.
  It uses shortened dual checks, an independently generated binary error ball,
  six-plus-six Z4 enumeration, packed support products and expanded integer signs.
- `source/negative_controls.py` checks positive examples and rejects incorrect
  ring lifts, duplicate supports, missing checks, overlapping checks and wrong
  trace parity in both implementations.
- `data/z8_generator.json` gives the complete ring generator and polynomial factors.
- `data/supports.npz` pins the 26,022-support family from the preceding verified
  dual-obstruction construction, including both generator matrices and all witnesses.
- `data/replay.json` and `data/independent.json` are the two complete receipts.
- `data/comparison.json` distinguishes published inputs, derived comparisons and
  preceding repository values. `data/SHA256SUMS.txt` covers every package file.

Each replay independently checks all **16,777,216 Z4 words**, all **759** binary
Golay octads against **5,489** small errors each, and all **338,559,231** distinct
extra-support pairs. The coordinate products are exact sums of at most 96 zeros
and ones in binary32; no approximate separation margin is accepted. The second
support replay uses unsigned integer masks. Section counts and field operations
use integer arithmetic.

## Why the compressed Z8 distance proof is complete

Let C be the free rank-12 Z8 code given by the cyclic degree-11 polynomial

```text
7 + x + 4x^2 + 4x^3 + 2x^4 + 7x^5 + 3x^6 + 3x^7
  + 4x^8 + 7x^9 + 2x^10 + x^11,
```

extended by minus the coordinate sum. The supplied complementary factor
certifies the factorization of 1+x+...+x^22 modulo eight. Its first twelve
columns are triangular with odd diagonal, proving freeness and uniform
projection onto any initial m≤12 ring coordinates.

The Gray map sends the three binary digits of a ring symbol to its affine
Boolean function on two variables. Its eight images are exactly the even
four-bit words. Directly checking all 64 symbol pairs proves isometry for the
homogeneous weight: zero has weight zero, four has weight four, and every
other nonzero symbol has weight two.

For a nonzero ring word:

1. If its binary residue is zero, it is twice a word of the Z4 reduction.
   A complete enumeration proves that reduction has minimum Lee weight 12,
   giving homogeneous weight at least 24.
2. If the binary residue has weight at least 12, its odd coordinates alone
   give homogeneous weight at least 24.
3. The remaining residues are the 759 weight-eight Golay words. For each,
   choose its binary-coefficient lift u. Outside the octad, every possible
   ring word has the form u+2v with v in the Z4 reduction. Both replays exclude
   Lee weight at most three in that sixteen-coordinate coset. These exclusions
   give at least eight additional homogeneous weight units.

The first proof computes canonical remainders for all low-weight errors. The
second constructs five shortened dual checks and shows that every such error
has a different syndrome from the target; it only needs the verified
orthogonality of those checks. Thus neither proof relies on a search solver or
an unenumerated list of ring words.

The Gray image has 8^12=2^36 words at distance at least 24. Flipping the first
bit in every four-bit row gives a disjoint second code. Every cross pair differs
in at least one position of each of the 24 rows, so their union has 2^37 words
and minimum distance 24. Mapping bits to signs gives squared norm 96 and
maximum distinct inner product 48.

## Sections and additional layers

Write S_i for the sum of the four physical coordinates in row i. The first
m≤12 ring-coordinate projection is uniform. In the unflipped coset, one row
has sums −4,0,4 with multiplicities 1,6,1; in the flipped coset it has sums
−2,2 with multiplicities 4,4.

- In dimension 96−r, r=4,...,8, impose S_0=...=S_(r−1)=0. The dense count is
  6^r * 8^(12−r). This gives dimensions 88–92.
- For dimensions 93–94 impose equality of the first m=97−d row sums. The
  dense count is (6^m+2+2*4^m) * 8^(12−m).
- In 95D impose S_0+S_1=0, retaining 70*8^10=75,161,927,680 dense points.
- In 96D use the entire 2^37-word union.

The displayed constraint matrices have rank 96−d, certified by nonzero minors
modulo two and independently modulo three. Orthogonal identification of their
kernels with R^d preserves norms, distances and distinctness.

The weight-24 points have entries ±2. Transversal supports come from the
explicit GF(4) [24,9,12] code; the additional 26,022 supports carry the
explicit binary [24,14,6] sign code. Every additional support pair has
intersection at most twelve. Its one or two disjoint odd trace obstructions
also bound its intersections with all transversals by twelve. All those facts
are checked again here; earlier search receipts are not used as proofs.

For 88–94D, signs are assigned to each support in increasing physical-coordinate
order and retained exactly when their row sums satisfy the section equations.
The first replay counts profiles of row occupancies. The second expands the
16,384 sign words into integer coordinates and checks the row-sum conditions.
For 95D use the documented support-dependent sign permutation, aligning an
intersection of size eight with an exhibited dual octad. The verified retained
counts for intersection sizes 0,...,8 are
`[16384,0,8192,0,6144,0,5120,0,8960]`.

A final weight-six layer has entries ±4 and 32 signs per support. Reed–Solomon
[6,4,3] codes over GF(9), GF(11), GF(13) or GF(16) give q^4 supports. In 88–94D
these use only coordinates outside the selected rows. Add the twelve A3 roots
in each selected four-cell row, scaled to squared norm 96, and coordinate poles
outside those rows. Their number is 192+4m. In 95–96D keep the earlier exact
GF(16) section and coordinate poles. Both implementations verify all field
and cross-layer conditions. Layers with different numbers of nonzero physical
coordinates are disjoint, and the transversal and additional weight-24 supports
have different occupied-row counts.

The coordinate recipes are finite and explicit. A dense point is the sign
image of a ring word, optionally with the first bit of every row flipped.
For each additional weight-24 support, use its occupied physical coordinates
in increasing order as the sign-code positions. In the 95D case with eight
selected cells, map those cells to the lexicographically least dual octad;
map the remaining positions in increasing order. For a six-support word
(a_0,...,a_5), use cell q*r+a_r of six disjoint q-cell rows. In 88–94D place
those rows after the first 4m coordinates. In 96D use q=16 and the natural
ordering. In 95D order the eight cells selected by [2,2,1,1,1,1] first, then
all remaining cells, to define the required physical permutation. Six-support
signs have odd parity. In 88–94D pure roots are +/-sqrt(48)*(e_i-e_j) within each selected
four-cell row; the outside poles are +/-sqrt(96)*e_j. In 95–96D only
the stated coordinate poles are included.

## Attribution and comparison scope

The large dense-code increase is a recovery of **published prior work**. Simple
shortenings already give 2^33, 2^34, 2^35, 2^36 and 2^37 sign-code points in
92–96D, respectively, and these now replace the obsolete comparison floors.
The older repository claims in those dimensions did not beat these inputs.

For an additional comparison, substituting that published 96-bit parent into
the standard Edel–Rains–Sloane chain (96,24,6,1) gives **141,736,018,112**.
The present 96D bound exceeds this derived count by **426,344,448**, exactly
the already verified additional support layer. The parent code, standard chain
and additional support family have separate provenance; none is presented as
an independently discovered new binary code.

The row-zero template stops improving the present bounds below 88D. In 87D
its dense layer has only 5,159,780,352 points. Even retaining all 426,344,448
additional-support points, all 2,097,152 old six-support points and all 228
prescribed pure points would give at most 5,588,222,180, below the existing
87D bound. This is a limit of that fixed template, not a kissing-number upper
bound or a claim about other section choices.

The section choices, complete counts and compatible small layers are explicit
here. The manuscript and PDF remain at their pinned earlier snapshot.
