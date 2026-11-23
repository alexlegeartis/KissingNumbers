# Atlas-compatible Delsarte–Goethals replacements in dimensions 65–68

| Dimension | Before the DG split | New lower bound | Increase |
|---:|---:|---:|---:|
| 65 | 419,000,274 | **422,858,578** | 3,858,304 |
| 66 | 429,290,020 | **432,812,452** | 3,522,432 |
| 67 | 448,616,838 | **451,803,398** | 3,186,560 |
| 68 | 481,044,376 | **481,928,984** | 884,608 |

Relative to the simpler DG split that deletes all of the 65D atlas, the
65D bound also increases from 422,678,354 by 180,224. This is a lower
bound, not an optimality claim.

## Construction

Start with the [complete 126-hyperplane atlas](../dim65-68-complete-hyperplane-atlas/)
and its actual RM80 dense and sparse stock. Retain all sparse support
masks and their 2,048 words per mask, but jointly rephase their signs.
Retain the later best four-site layers, using
[the optimal 65D packing](../dim65-optimal-quadruple-layer/),
[the grouped 66D design](../dim65-66-grouped-four-site-designs/), and
[the nested 67–68D trades](../dim66-68-nested-affine-face-trades/).
Keep the poles on the tail axes only.

Partition the 64 body sites into the 32 disjoint pairs in
`assignment.txt`, each assigned a distinct cubic label from the old
RM80 palette. In each assigned cubic label select 1,024 Delsarte–Goethals
quadratic parts and all 128 affine shifts. This is the same length-64
minimum-distance-24 code as in the basic DG split. Replace each
selected owner by two vectors, with swapped magnitudes `(a,b)` and
`(b,a)` on its pair, where `a=28/5` and `b=29/20`. Flip the
owner's sign at the large magnitude only. The other 62 body coordinates
keep magnitude one and the tail is zero. Thus the replacement vectors
all have squared norm `N=38185/400`.

The exact donor–donor caps for the same owner, distinct owners at the
same pair with matching or swapped profiles, and owners at different
pairs are respectively `62−2ab`, `16+a²+b²−2`, `16+2(ab+1)`,
and `32+2(a+b)`. Against retained dense, four-site, and pure tail
pole vectors the caps are `32+a+b`, `2+a+b`, and 0. Against a
weight-16 sparse vector with zero tail containing the large site,
the support is orthogonal modulo two to every old owner, but the
donor has one sign flipped. Its parity differs from the new even
sparse parity, so at least one ordinary sign disagrees. The caps
with both sites or only the large site present are `12+a+b` and
`13+a`. The small-site-only and neither-site caps are `15+b`
and 16. A sparse vector with nonzero tail has exactly 14 body sites,
so the maximum raw dot is `12+a+b` with no parity condition.

## Retain 44 complete hyperplane fibres

The last line of `assignment.txt` identifies the retained 44
hyperplanes in the order of the published atlas witness. All head
types and all sign words on those hyperplanes are kept. The new
angular test compares every allowed body-sign difference between a
donor owner and an atlas vector on its 32-site hyperplane. The
difference word is in the coset of RM(2,5) determined by the
published hyperplane shift and the restricted cubic label. (All
selected owner additions and atlas head phases restrict to
quadratic polynomials.) On every retained hyperplane this coset
has weight at least four. For any weight-four word, an exhaustive
check shows neither site of the designated pair for that cubic
label is among the four disagreements.

For weight four, the unweighted body dot in units of sqrt(2) is
24. At a large site inside the hyperplane the original sign
necessarily agreed, so flipping it decreases the dot. A small
site inside the hyperplane increases the bound by at most
`b−1`. For weight at least five, the unweighted bound is 22
and the two modifications together add at most `a+b`.
Both `sqrt(2)(24+b−1)` and `sqrt(2)(22+a+b)` meet
the angular gate for squared norms `N` and 72.
The donor has no tail term in either comparison.

The checker also verifies the inherited hyperplane–hyperplane
parity differences and hyperplane internal distances, and
constructs the signed sparse-to-atlas guard rows. Each
`edgebits` bit in `phase68.json` simultaneously toggles the
two opposite parities of one shared sparse overlap row, without
changing their opposition. The checker starts the new phase with
even total parity for every zero-tail support, retains the
inherited total parity when a tail is present, and checks that
every sparse overlap and atlas check belongs to a consistent
balanced local binary dual. Where the hyperplane contains all
16 sparse sites, it preserves the full-mask parity as well, so
the required one disagreement is forced to be at least two.
As in the original proof, balancing the local coordinate
classes extends each imposed span to a rank-five dual, giving
exactly 2,048 sign words of minimum distance four. All 67,923
masks and 31,512 shared sparse guards survive jointly.

## 352 four-tail sparse banks in 68D

The complete atlas had dropped 128 earlier four-tail banks. With just
44 hyperplanes, 112 of those old banks are compatible again. Find
further banks by classifying cubic labels under invertible affine
coordinate changes modulo quadratic terms. Select labels whose
binary cosets have weight-12 words and are sufficiently far from
all seven original sparse cubic stocks. An exact scan against the
additional signed witnesses and retained hyperplanes leaves another
240 compatible masks. Append all four tail sites to each of the
352 distinct 12-site body masks in `restored68.txt`.

The replay needs no orbit classification or search. It verifies
every sparse overlap including these 352 added supports. Each new
bank has no dangerous sparse overlap and admits all required head
guards in its even-parity balanced local dual. Each thus supplies 2,048
sign words. Their body has only twelve sites, so their mixed donor
dot is at most `10+a+b` (the conservative earlier gate also
suffices). The old non-donor dense, four-site and pole magnitudes
continue to satisfy their gates. This adds 720,896 points beyond the
44-hyperplane variant alone and leaves all lower-dimensional sections
unchanged.

## Replay

Run `python3 -B verify.py` with a C++20 compiler on a BMI2-capable x86-64
machine. The input hashes pin the existing cubic palette, sparse and
atlas witnesses, and relevant predecessor claims. The C++ checker
reconstructs the old and new signed checks and every sparse
overlap, and tests the exact atlas/cubic cosets. The Python driver
reconstructs the DG family and checks every mixed bound with integer
or rational cross multiplication. Each replaced owner adds one net
vector. The deletion costs are the 128 body poles and 82 full
hyperplanes, at 2, 4, 6 and 24 heads per hyperplane in 65–68D.
No numerical optimization, network access or solver is needed for replay.
