# Larger compatible sign families in dimensions 61–67

**13 November 2026.** This package gives exact finite improvements to the
code constructions. The preceding values are the repository certificates at
commit `d595c660a8672e0ee436935b74790dadb5d9a838`.

| Dimension | Preceding bound | New bound | Gain |
|---:|---:|---:|---:|
| 61 | 88,260,394 | **88,975,146** | 714,752 |
| 62 | 121,824,268 | **122,539,020** | 714,752 |
| 63 | 184,286,468 | **184,393,988** | 107,520 |
| 64 | 335,496,064 | **335,711,104** | 215,040 |
| 65 | 335,564,290 | **335,779,330** | 215,040 |
| 66 | 341,994,884 | **342,271,364** | 276,480 |
| 67 | 355,090,822 | **355,095,446** | 4,624 |

The finite construction, its size and every relevant interaction are checked
exactly. The search is not an optimality proof. Dated literature comparisons
remain separate from the geometric proof.

## Reproduce

Python 3, NumPy, threadpoolctl and a C++17 compiler named `g++` on the path are
required. No optimizer, candidate pool or discovery checkpoint is needed.
From this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python rejection_controls.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python final.py
```

`verify.py` checks the manifest and pinned source files, independently rebuilds
the binary parent generator, enumerates all of its 268,435,456 words in the
scalar checker, checks the sign families in two representations, and replays
the dimension transfers. It also freshly reconstructs the 65D owner class and
the Netto four-support design. Its output is recorded in `verification.json`.
`final.py` derives the registered counts; it is a count driver, not a substitute
for the full proof.

## 1. Points, sign codes and the repair rule

Work at squared norm 64, with allowed inner product at most 32. A weight-sixteen
support carries entries ±2. Its 2,048 signs form an affine translate of an
extended Hamming `[16,11,4]` code. The five parity checks consist of the all-one
row and four further rows whose sixteen columns are the sixteen distinct odd
five-bit vectors. This proves rank five and excludes every word of weight one,
two or three from the kernel. Distinct signs on the same support therefore
give inner product at most `64−8·4=32`.

For two supports with intersection at most eight, the amplitude bound gives
inner product at most `4·8=32`, independently of signs. For a ten-coordinate
intersection, we require a common weight-eight dual check contained in the
intersection, with opposite affine parities on that check. At least one sign
must then differ in the intersection, and

\[
  x\cdot y\le4(10-2)=32.
\]

This same rule is applied to **new/new pairs as well as old/new pairs**.
Every ten-overlap is enumerated and checked; larger overlaps are rejected.
The code does not assume that a support's score or a searcher's compatibility
flag establishes this condition.

The original 30,828 weight-sixteen supports come from a reconstructed binary
`[64,28,16]` code. The scalar replay enumerates all `2^28` words, proves the
minimum, and confirms the exact weight-sixteen list. Thus any two original
supports intersect in at most eight coordinates, even when their Hamming sign
codes have been reordered or shifted.

## 2. Three explicit full-class families

Each JSON row gives the physical 64-bit support, five local parity checks and
a local affine shift. An unlisted original support retains a standard Hamming
sign class, with the total parity specified by the construction below.

| Certificate | Original supports retained | Added supports | Use |
|---|---:|---:|---|
| `data/short61_62.json` | 9,673 | 937 | 61D and 62D |
| `data/global64.json` | 30,828 | 1,940 | 63D, 64D and 65D |
| `data/expanded66.json` | 30,828 | 1,938 | 66D |

The shortened family was optimized after removing coordinates
`24,42,53,55`; old supports meeting those coordinates impose no constraints
on that search. All retained sparse sign classes have odd total parity.
This replaces the earlier 588 added supports by 937, a gain of
`(937−588)·2,048=714,752` points in both dimensions.

The global family has exactly `30,828+1,940=32,768` supports. Its sparse layer
therefore contains exactly `2^15·2^11=2^26` points. The full 64D count is

\[
2^{28}+2^{26}+10\,416\cdot16+128=335\,711\,104.
\]

The eight initially imposed old/new checks on each discovered support have
rank four, including the all-one check. Completing that system to rank five
leaves a useful choice of the last Hamming check. The discovery search uses
that freedom to satisfy additional new/new repairs. The certificates contain
the completed codes, so replay does not depend on this search description.

## 3. Exact transfers to 61D and 62D

Write the four removed bits in the order `24,42,53,55`. Their projection from
the parent code has rank four. The remaining 60 coordinates form the body.

For 61D, shorten the last two removed bits to zero. A word whose first two
removed bits are `(i,j)` gives a unit vector with signed body amplitude
`1/8` when `j=0`, and `sqrt(17)/32` when `j=1`. Its last coordinate is
`(−1)^i/4` or `(−1)^i/16`, respectively. Norms are exactly one. Parent distance
16 bounds every body distance, and all resulting radical inequalities are
checked by exact rational squaring. There are `2^26` dense points. Since each
retained sparse support is dual to the parent code and has odd sign parity,
its dense/sparse products have the stronger fourteen-term bound required by
the second amplitude.

For 62D, keep six fibres with three-bit labels `0,1,3,7,6,4` and put their heads
at the regular hexagon `(2,0),(1,sqrt(3)),(−1,sqrt(3)),(−2,0),
(−1,−sqrt(3)),(1,−sqrt(3))`. The last removed bit is zero. The norm-64 body
has entries ±1. Each fibre has `2^24` points and all 36 head-pair types satisfy
the exact distance bound.

The independent four-support layers and poles give the totals

\[
\begin{aligned}
\tau(61)&\ge2^{26}+10\,610\cdot2\,048+8\,555\cdot16+122,\\
\tau(62)&\ge6\cdot2^{24}+10\,610\cdot2\,048+9\,145\cdot16+124.
\end{aligned}
\]

Every four-support block is constructed from affine planes in `F_2^6`.
Unique triples prove that distinct blocks intersect in at most two positions.
All dense, sparse, four-support and axial cross bounds are checked explicitly.

## 4. The 63D balanced section

For each of the 72 supplied dual octads `H`, take the hyperplane perpendicular
to its indicator. The dense projection has rank seven and is the even-weight
code on eight positions, retaining `70·2^21=146,800,640` points.

On a sparse support `S`, let `m=|S∩H|`. An odd `m` gives zero. An even `m<8`
gives `binomial(m,m/2)·2^(11−m)`. If `m=8`, the contribution is 560 unless
the octad is a dual check; in that case it is 1,120 for even phase and zero
for odd phase.

**Untouched original Hamming classes are reordered separately for the chosen
octad.** When `m=8`, place that intersection on a dual hyperplane to retain
1,120 signs. This changes no repaired pair: every original support with a
ten-overlap is explicitly listed among the modified classes. The other
original/original and original/new intersections are at most eight, so their
validity is independent of sign order. The four-support design may also be
labelled independently, sending this octad to coordinates `0,...,7`.

The chosen section has terms

\[
146\,800\,640+35\,503\,776+1\,986\,560+102\,900+112
=184\,393\,988.
\]

`section63.py` independently enumerates the actual affine sign words on every
modified and added class, and checks the standard sign words for every
untouched class type. It does not rely only on the projection formulas.

## 5. The 65D owner lift and the 66D extension

For 65D, use the existing BCH radical dense code in place of the 64D dense
code. This substitution preserves all dense/sparse products because their
bound depends only on the amplitudes. The owner construction is freshly
replayed: 65,536 dense owners and 336 sets of eight sparse owners give 68,224
owners. None is in the weight-sixteen layer. Replacing each owner by two caps
and adding two poles gives a net gain of `68,224+2` above the new 64D count.
All new cap/equator bounds are checked exactly.

For 66D, the selected 73-bit weight-sixteen support list leaves 34,028 supports
after shortening. The first 64 coordinates reproduce the original 30,828
supports; 3,200 further supports have nonzero tail. Each of the 1,938 added
supports meets each of those 3,200 supports in at most eight coordinates.
The vector and scalar representations check all **6,201,600** such pairs.

The scalar checker also checks all **578,935,378** unordered pairs of retained
parent supports directly. Consequently the 66D proof does **not** assume the
minimum distance of the full 73-bit code or its full weight enumerator.
Membership is checked, but completeness of that parent list is unnecessary
for a lower bound. The Netto construction is regenerated from its explicit
projective orbit; its 11,088 retained four-support blocks and all poles give

\[
2^{28}+(34\,028+1\,938)\cdot2\,048+11\,088\cdot16+132
=342\,271\,364.
\]

## 6. Smaller affine classes give a further 67D improvement

The full 2,048-sign search did not improve the existing 1,792 added support
classes in 67D. Allowing smaller affine subcodes gives **49 additional classes**
with a total of **4,624 points**. Each is an affine subcode of a Hamming code;
all its linear differences have weight at least four. Their sizes range from
16 to 512. The certificate explicitly gives every check and affine shift, and
the verifier enumerates the entire kernel and checks its minimum distance.

The added classes impose further weight-eight parity constraints on 357
previously untouched parent classes. Their checks are jointly completed to
full Hamming duals with consistent phases, so each affected parent class
retains all 2,048 signs. Every pair of the 49 new supports has intersection at
most eight. Every larger intersection with a parent support or a previously
added support has a common eight-coordinate check of opposite parity.
The verifier checks every affected support pair and all changed sign codes.

`baseline67.py` also freshly verifies the unchanged 1,792-class construction,
all 17,203,200 interactions with its 9,600 nonzero-tail supports, and all
**817,191,378** retained parent-support pairs. Together with the reconstructed
Netto four-support layer, this proves the baseline 355,090,822 before the
4,624 new points are added. The full package therefore does not rely on an
unreplayed 67D population claim or the minimum of the full 73-bit code.

This extension changes the optimization variable from the number of full
Hamming classes to the sum of the sizes of compatible affine sign classes.
The retained 49-class witness is a finite lower bound, not a claimed optimum
within that wider family.

## 7. Verification files and scope

- `verify_family.py` reconstructs the parent generator and checks sign kernels,
  distinctness, exact packed/coordinate intersections and every repair.
- `check_family.cpp` is a separate scalar implementation; it freshly enumerates
  the parent words and checks all old/new and new/new pairs.
- `verify_independent.py` independently constructs sign kernels, checks both
  selected families, explicitly recounts the 63D section and checks the
  compression geometry. `independent_verification.json` records its result;
  `INDEPENDENT_REVIEW.md` states the reviewed proof scope.
- `geometry_counts.py` proves the 61D/62D radical geometry and common cross-layer
  bounds, constructs the affine four-support layers and derives counts.
- `section63.py` checks the octad sections and independently enumerates signs.
- `expanded_geometry.py` and `check_expanded.cpp` check the 66D support transfer
  and every retained parent-support pair.
- `replay_parents.py` reruns the unchanged BCH owner construction and constructs
  the Netto four-support design, without modifying source files.
- `baseline67.py` freshly verifies the old full-class 67D construction.
  `subcode67.py` verifies the additional affine classes in `data/subcode67.json`,
  every changed Hamming class and every affected support interaction.
- `verify_subcode67_independent.py` brute-filters all 65,536 sign words for each
  changed code and independently checks 115,998,720 exact sign-pair products
  on the 758 critical overlaps. It supports `--controls` for eight adversarial
  checks. `independent_subcode67.json` and `subcode67_rejection_controls.json`
  are its deterministic positive and rejection receipts.
- `rejection_controls.py` rejects duplicated supports, broken check matrices,
  lost odd parity, wrong old/new and new/new phases, and optimized Python.
  Both phase mutations exhibit actual norm-64 points of inner product **40**,
  exceeding the allowed 32. The scalar and Python checkers reject them.
- `verify.py`, `final.py`, `verification.json`, `parent_verification.json` and
  `rejection_controls.json` are the entry points, count driver and replay records.
- `data/selection.json` chooses the explicit families. `data/source_inputs.json`
  pins the exact source files. `SHA256.json` covers the complete package.

Discovery programs and finite unsuccessful searches are archived under
[`research/dim32plus_continue_20261113/codes`](../../../research/dim32plus_continue_20261113/codes/).
The proof uses the explicit data and exact checks above, not an optimization
solver's status or a proposed search optimum.
