# Dimensions 88–96: nonlinear Z8 Golay codes and exact sections

This package corrects a missed published input and constructs stronger repository
bounds in nine dimensions. The nonlinear binary **(96, 2^37, 24)** parent is due
to **I. M. Duursma, M. Greferath, S. N. Litsyn and S. E. Schmidt (2001)**,
[IEEE Transactions on Information Theory 47(4), 1596–1598](https://doi.org/10.1109/18.923742).
It is not a new code discovered here. The former use of the linear [96,33,24]
code missed this stronger nonlinear input. The corrected published-code baselines
are recorded separately from the gains over the preceding repository checkpoint.

| dimension | verified lower bound | gain over preceding certificate (88–94: November row sections; 95–96: original update) |
|---|---:|---:|
| 88 | **6,967,552,256** | +25,315,328 |
| 89 | **9,270,596,604** | +27,257,856 |
| 90 | **12,341,351,672** | +24,709,396 |
| 91 | **16,435,507,572** | +19,134,548 |
| 92 | **21,896,439,920** | +12,449,792 |
| 93 | **31,056,799,856** | +12,272,924 |
| 94 | **47,703,347,340** | +4,408,320 |
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
- `source/adapted_common.py` supplies structural pattern guards and physical sign-position maps.
- `source/adapted_replay.py` and `source/adapted_independent.py` independently
  count the support-adapted zero-row sections and check their support maps.
- `source/adapted_equal_replay.py` and `source/adapted_equal_independent.py`
  independently count the support-adapted equal-row sections and check their maps.
- `source/adapted_negative_controls.py` rejects malformed adapted patterns,
  nonbijective maps and corrupted section claims.
- `source/section_negative_controls.py` rejects malformed row sections and
  corrupted counts, including an invalid uniform-projection count for the octad section.
- `nonlinear/source/guards.py` checks nonlinear source families, patterns and support maps.
- `nonlinear/source/independent.py` independently reconstructs the nonlinear
  source code and verifies its distances, patterns, support maps and section counts.
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
these use only coordinates outside the selected rows. Add the twelve A 3 roots
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


## November 2026 alternative sections

The current 88–94D results replace the original first-row sections below with the explicitly listed row sets in `data/sections.json`. The original parent proofs and their receipts remain unchanged. `data/section_receipt.json` supplies the new seven results; `final.py` combines these with the unchanged 95D and 96D results. These are improvements over the repository baseline, with no claim of worldwide priority or optimality.

### Construction and proof

## Verified improvements against the fetched repository baseline

The construction is the existing 96D parent from `alexlegeartis/KissingNumbers/verifications/improved/dim88-96-z8-golay`, with different linear sections. The published nonlinear binary parent remains Duursma–Greferath–Litsyn–Schmidt (2001), DOI 10.1109/18.923742. No new parent code is claimed.

All coordinates are grouped into 24 rows of four. For the rows T listed in `data/sections.json`, impose that each selected row sums to zero. This imposes exactly |T| independent real linear equations. The construction consists of:

1. Unflipped Z8-code Gray words whose selected ring coordinates are nonzero modulo 4, mapped to signs in R^96. Flipped-parent words cannot satisfy a zero row sum.
2. Every archived weight-24 support/sign vector satisfying these equations. Signs are assigned to occupied positions in their original increasing-coordinate order. No alteration of the source sign code is used.
3. The source Reed–Solomon weight-six layer, moved to any 6q distinct coordinates outside all selected rows, using increasing outside-coordinate order. The size is unchanged. This move preserves all compatibility: its supports intersect one another in at most three, and its inner product with any weight-24 vector is at most 6*4*2=48.
4. All source pure roots: 12 A 3 roots in each selected row, plus both coordinate poles outside selected rows. Total 192+4|T|.

All vectors have squared norm 96 and pairwise inner product at most 48 by the source parent proof and unchanged layer compatibility bounds. Taking a linear section never weakens these bounds. An orthonormal basis of the section gives the claimed Euclidean dimension. The transversal weight-24 layer contributes zero in these sections.

## Count checks

The selected witnesses were found by exact support/sign enumeration and a superset zeta transform over row masks. Search optimality is not needed for validity and the large search histograms are not shipped.

`source/section_replay.py` processes all 26,022 supports without profile grouping using Python big-integer intersections of sign-condition masks. `source/section_independent.py` independently uses matrix-expanded signs and coordinate prefix sums. Both enumerate complete projected ring images. Since the generator is free over Z8 of rank 12, every projected image has equal fiber size.

In 88D the selected eight rows form a binary Golay octad. Their Z4 image has size 32,768, with 3,281 fully nonzero vectors. This adds 1,048,576 dense points relative to an independent eight-row projection. For the four through seven selected rows the Z4 projections are all surjective.

Results:
- 88D: 6,922,487,106, increase 1,235,522.
- 89D: 9,229,785,598, increase 431,462.
- 90D: 12,307,197,352, increase 527,024.
- 91D: 16,410,365,876, increase 766,760.
- 92D: 21,881,382,480, increase 892,192.

The comparison baseline is the source `bounds.json`, also checked against the repository RESULTS.md. These are repository improvements, not an independent claim of worldwide novelty. Both inherited parent proofs and both new section count proofs are replayed by the production verifier.

### Equal-row sections and production verification

For 93D impose S0=S1=S4=S8. For 94D impose S0=S1=S8. Their dimensions are 96-(4-1)=93 and 96-(3-1)=94. Selected Z8 projections are surjective; explicit enumeration checks both Gray-parent cosets. The transversal layer contributes 536,870,912 and 1,073,741,824, respectively. These counts are also checked directly on every sign word. The RS layer uses the first 6q outside cells in increasing order; A 3 roots lie within selected rows, and poles lie outside.

Run `python verify.py` from this package. It checks the manifest, both complete parent proofs, both new exact section-count implementations, all layer counts and the explicit unit minors proving section dimensions. It compares computed receipts with pinned data. The candidate-specific corruption controls reject 64 malformed counts/sections, including a duplicate selected row, wrong rank, wrong sparse/dense/transversal/RS/pure counts and the erroneous uniform-projection formula for the 88D octad.

Both new counting implementations are independent: one uses integer sign masks and set-based ring image closure; the other uses matrix-expanded sign words, coordinate prefix sums and unique-array ring images. They share only final geometry/claim consistency checks. The inherited parent implementations remain fully independent. See `SECTION_AUDIT.md` for the independent mathematical audit; its packaging-gap paragraph describes the pre-integration candidate and is resolved by these production files and controls. The construction is an injective isometric identification of a linear section, not coordinate deletion.

## Support-adapted signs: latest November continuation

The top table and `bounds.json` now report the further adapted-sign improvements, compared with the preceding published row sections. The earlier chapters and receipts are retained as historical proofs. Current witnesses are `data/adapted_zero.json` and `data/adapted_equal.json`; their sign-position patterns are `data/adapted_patterns.json` and `data/adapted_equal_patterns.json`.

For each weight-24 support, inspect occupancy in the selected rows. For zero sections reject odd occupancy. Order the nonempty selected rows by occupancy (two before four), then by row index. Use the corresponding pattern's ordered groups of sign positions on those occupied cells, sorted by physical coordinate within each row. For equal sections order all selected rows by occupancy then index, including empty rows; mixed odd/even occupancies cannot have equal sums. Map remaining occupied cells in increasing order to remaining sign positions in increasing order. This is a full bijection of 24 positions, explicitly reconstructed for every included support by the independent proof.

Retain every codeword satisfying zero/equal sums after this permutation. The same [24,14,6] sign code is used on each support. Within-support distances are preserved by the bijection. Across distinct supports, intersection at most12 guarantees product at most4*12=48 for any signs, so different supports need not share a permutation. Dense/sparse, sparse/RS and sparse/pure bounds likewise hold for arbitrary signs. All vectors still have squared norm96. Support sizes separate layers, and within-support bijectivity ensures distinctness. The operation remains an isometric linear section, not coordinate deletion.

Dual octads explain why adaptation helps: four pair constraints can have a redundant parity equation, doubling their retained sign count. More general patterns are counted exactly, not inferred from parity heuristics. The final88D selected rows are [0,5,8,10,13,18,19,21]; their projected Z4 image is surjective, with 65,536 images and6,561 all-nonzero images. Thus the dense count is6,879,707,136, while62,319,616 adapted sparse vectors give the new total6,942,236,928. The previous octad bonus is not mistakenly retained.

The latest93D equations are S0=S1=S8=S16;94D uses S0=S1=S8. The unchanged transversal layer is counted directly. RS points use the first6q outside coordinates. A 3 roots and outside poles certify the section span; the equal sections additionally contain the all-positive dense word for the common-sum direction.

`verify.py` runs both full parent proofs, both preceding section proofs, both primary adapted proofs and both independent physical-permutation audits. It checks immutable claims, every pattern count, every support mapping, dense images, all remaining layer counts, constraints and inherited pair-category arguments. Parent support data and parent proof code are unchanged. New rejection controls cover repeated or out-of-range sign positions, wrong group sizes, nonbijective physical assignments, parity-incompatible claims, duplicate rows and corrupted layer counts. All artifacts are manifest-pinned. `AUDIT.md`, `FINAL_88_AUDIT.md` and `EQUAL_AUDIT.md` preserve the independent research audits; their packaging-gap notes refer to the pre-integration state and are resolved by these production files.

The pattern search was heuristic. Even exhaustive row scans establish optimality only for the fixed stored patterns, not among all sign permutations or kissing configurations. No worldwide-priority claim is made.

## Latest continuation: nonlinear sign fibers

The current top-table 88–94D values are the nonlinear continuation in [nonlinear/README.md](nonlinear/README.md). It preserves every earlier parent/section/adapted-sign certificate and uses a classical Nordstrom–Robinson input, not a new parent-code discovery. The latest immutable claims are `nonlinear/data/verified.json`. `final.py` overlays them, while95–96D remain unchanged. All earlier count tables and gains elsewhere in this README are historical stages.

The full production verifier now additionally runs two exact nonlinear-code proofs, both checking all 8,386,560 pairs of the4,096-word source family, every pattern, every physical support map, all dense projections and ranks. Candidate-specific controls reject corrupted octacode generators, translates, partitions, span bases, source words, families, translations, pattern counts/groups and final layer claims. The original support file and both parent proof programs remain unchanged. The root SHA256 manifest covers the complete nonlinear subpackage.
