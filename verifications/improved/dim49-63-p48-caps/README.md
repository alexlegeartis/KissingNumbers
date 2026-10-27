# Dimensions 49–63: the cap construction over P₄₈

**Retained construction and history.** The October registry takes the strongest
construction in each dimension. See the [current registry](../../../RESULTS.md),
[code-layer updates](../dim61-68-95-96-code-layers/) and
[nearest-coset updates](../dim73-88-gamma72-cosets/). Older counts below remain
results about the original construction.

The **28 October refilling continuation** improves every registered bound in
49–60D. Its base contains **7,678 lines** and its 577 classes retain
**4,047,491 distinct mother lines**. `october/family_refilled.npz` exhibits
the new base class, 577 exact isometries and every retained base index.
For compact storage, isometries use generator words or references to the retained
historical transforms, and index sets use packed bitmaps. `october_family.py`
expands them with bounded integer products before verification. Its exact counts and gains over commit
`a6ba3387f2b0b54942812790a1ddf19bb4e99431` are in
`october/family_refilled.verified.json`. `verify_refilled.py` checks every base
pair and all isometries twice, and checks disjointness using both canonical
integer hashes and independent full 48-coordinate comparisons. The binary64
replay first bounds every product and partial sum within its exact integer range.
`verify_october.py` runs the historical replay and both new representations.

The new replay and family generator require NumPy and threadpoolctl.

```bash
python verify_october.py
python verify_refilled.py
python scripts/refill_family.py --output /tmp/family_refilled.npz
```

The last command deterministically regenerates every array of the installed
family from its exhibited base. The family search considers the original 577
images plus 2,000 distinct seeded random-word images, chooses by marginal
coverage, and refills after every selection. The base was discovered by filling,
one-for-two replacements and neutral swaps in a 12,854,287-line candidate pool.
The timed discovery search is not required to trust or replay the certificate.
The P48 minimum six and 52,416,000-point parent shell remain literature inputs.
The unchanged direction and pole codes retain their earlier verification scope;
60D uses the published 756-point Coxeter–Todd shell and the 841-point 12D pole code.
The new result is a verified improvement over the preceding repository
construction, with no unrestricted optimality or new worldwide priority claim.

The strategy follows the fill/correct principle in
[PackingStar's published method](https://arxiv.org/abs/2511.13391v4).
Its newer announcement discusses joint removals, midpoint insertions and
auxiliary layers, but does not specify a complete implementation of that new
scheme. Our fixed lattice coordinates enforce Gram feasibility automatically;
the search operates on the exact conflict graph. This is a CPU combinatorial
adaptation, with no reinforcement-learning reproduction claim.
[The pinned source review](../../../common/source_review_20261028_p48.json)
records that distinction and their larger reported 26–31D values.

Choosing whole image families jointly matters even before enlarging the base:
reordering the original 577 images and restoring overlaps recovers up to 2,223
lines in the prefixes used by 50–60D. A line removed because an earlier image
owned it may be retained when that earlier image is not selected. We rebuild
disjointness from the raw images rather than append those lines to an already
selected family. `october_family.py` compares complete families separately.

The 24 October family has a **7,387-line base**, **577** explicit isometric
images and **3,901,694** distinct mother lines. `verify_october.py` checks
all base pairs, isometries and disjointness; `october_family.py` selects
the best prefix of either family, without mixing them. The three
`october/family_55_61.*` files preserve the certificate and source reports.
`UPDATE_PROVENANCE.json` records source hashes; `data/SHA256SUMS.txt`
covers this package. The generated table below now uses the stronger family.

**Current scope.** This package supplies the registered bounds in 49–60D.
Its cap counts in 61–63D remain valid constructions and are now superseded by
the [October code constructions](../dim61-68-95-96-code-layers/). The older
[ERS chain](../dim62-63-ers-chain/) also beats its caps in 62–63D.
`scripts/ceiling.py` bounds families of images of the shipped bases, whose
largest class size is fixed by the chosen shipped certificate. Larger base classes remain open; no absolute
ceiling for this cap construction is proved.

```
python verify.py                 # ~1 min: the lattice, the class, the arithmetic
python final.py                  # instant: the result table
python scripts/regenerate.py     # ~15 min: rebuilds and re-verifies all 1400 classes
```

| dim | k | λ(k) | Σ|C_i| | gain | **this work** | previously |
|---|---|---|---|---|---|---|
| 49 | 1 | 1 | 7 678 | +15 358 | **52 431 358** | 52 416 002 |
| 50 | 2 | 3 | 23 034 | +46 074 | **52 462 074** | 52 416 006 |
| 51 | 3 | 6 | 46 062 | +92 136 | **52 508 136** | 52 416 012 |
| 52 | 4 | 12 | 92 078 | +184 180 | **52 600 180** | 52 416 024 |
| 53 | 5 | 20 | 153 320 | +306 680 | **52 722 680** | 52 416 040 |
| 54 | 6 | 36 | 275 391 | +550 854 | **52 966 854** | 52 416 072 |
| 55 | 7 | 63 | 480 069 | +960 264 | **53 376 264** | 52 416 126 |
| 56 | 8 | 120 | 906 389 | +1 813 018 | **54 229 018** | 52 416 240 |
| 57 | 9 | 153 | 1 149 624 | +2 299 554 | **54 715 554** | 52 416 306 |
| 58 | 10 | 255 | 1 884 834 | +3 770 178 | **56 186 178** | 52 416 510 |
| 59 | 11 | 302 | 2 215 274 | +4 431 152 | **56 847 152** | 52 416 604 |
| **60** | **12** | **378** | **2 738 706** | **+5 478 253** | **57 894 253** | 52 416 841 |
| 61 | 13 | 577 | 4 047 491 | +8 096 136 | 60 512 136 &nbsp;*(superseded by 71 310 488)* | 52 417 154 |
| 62 | 14 | 966 | 5 899 945 | +11 801 822 | 64 217 822 &nbsp;*(superseded by 100 371 468)* | 52 417 932 |
| 63 | 15 | 1 282 | 7 473 594 | +14 949 752 | 67 365 752 &nbsp;*(superseded by 182 407 428)* | 52 418 564 |

The retained caps in 49–61D clear the evaluated ERS chain, whereas the
62–63D caps do not. The comparison is computed by
[`ers_exposure.py`](../../closed/dim96-ers-takeover/ers_exposure.py).
The current code certificates are stronger still in 61–63D. These comparisons
do not exclude a larger cap base or a different construction.

Every value above stays below the dimension-64 record 331 737 984, as monotonicity of τ
requires — `final.py` checks this.

## The idea

In ℝ⁴⁸ ⊕ ℝᵏ, with every point a unit vector and every pairwise inner product ≤ 1/2, at cap
level *t*:

| family | points | count |
|---|---|---|
| equator | (v/\|v\|, 0), v a minimal vector of P₄₈ not used as a cap head | 52 416 000 − 2Σ\|C_i\| |
| caps | (s√t·û, √(1−t)·z), s ∈ {±1}, u a line of class *i*, z ∈ Z_i | 2Σ\|C_i\|·\|Z_i\| |
| poles | (0, w) | τ(k) |

A **class** is a set of minimal lines pairwise at \|cos\| ≤ γ. The four binding constraints
are

| pair | condition |
|---|---|
| same direction z, two lines of a class | t·γ + (1−t) ≤ 1/2 → t ≥ 1/(2(1−γ)) |
| same line, two directions | ⟨z,z′⟩ ≤ (1/2 − t)/(1 − t) |
| two classes, two directions | ⟨z,z′⟩ ≤ 1/2 |
| cap · pole | √(1−t)·⟨z,w⟩ ≤ 1/2 |

so the total is `52416000 + 2Σ|C_i|(|Z_i| − 1) + τ(k)`, and what matters is how many
directions one class line can serve.

### P₄₈ gets antipodal pairs, not zero-sum triples

Its class threshold is |⟨u,u′⟩| ≤ 2 of 6, i.e. **γ = 1/3**, which forces t ≥ 3/4; and at
t = 3/4 two directions on the same line need ⟨z,z′⟩ ≤ (1/2 − 3/4)/(1 − 3/4) = −1, i.e.
z′ = −z. So a line serves an antipodal **pair** and never a triple, giving 2 points per
line, and λ(k) is the number of pairs available — that is, half the size of the largest
ANTIPODALLY SYMMETRIC 60° code of ℝᵏ. That is not τ_lattice(k)/2: Cohn's record
configurations are antipodally symmetric for k = 9, 10, 11, 13, 14, 15 (checked by
[`scripts/verify_capdirs.py`](scripts/verify_capdirs.py)), so λ(k) = τ(k)/2 there, and only
k = 12 is excluded, τ(12) = 841 being odd.

**k = 12, checked rather than assumed.** λ(k) is half the size of the largest antipodally
symmetric 60° code, and the record being odd only rules out taking *all* of it — a symmetric
SUBSET could still have beaten the lattice's 756. It does not: of the 841 points of the
dimension-12 record, exactly **36 have their antipode in the set**, so the largest symmetric
subset is 36 points and λ would be 18 against the lattice's 378. Each pair is worth about
2 × 6 650 = 13 300 in dimension 60, so this was worth 480 000 had it gone the other way. The
841 points are Cohn's, at `../../closed/dim09-19-record-maximality/data/c841.npy`. Using τ_lattice(k)/2 cost seven of the fifteen
dimensions, dimension 63 most of all: λ(15) = 1282 rather than 1170. The Leech and Γ₇₂ both have γ = 1/4,
which is exactly the threshold that permits t = 2/3 and hence zero-sum triples worth 4
points per line — that is the one place they do better. See
[`../dim73-95-gamma72-caps/`](../dim73-95-gamma72-caps/).

### The class size is the whole game

The gain is 2Σ|C_i|, so everything turns on how large a class one can find. The conflict
graph on P₄₈'s 26 208 000 minimal lines is regular of degree n[3] = 36 848, and **Caro–Wei
guarantees only 712** lines — with no coordinates at all. An actual search on explicit
coordinates finds **7069**, a factor **9.93**. That factor is remarkably stable across
lattices: 11.3 for the Leech, 9.58 for Γ₇₂, 9.93 here.

The Caro–Wei version of this package gave dimension 63 a gain of +1 614 600; the explicit
class gives +13 847 064, a factor **8.58**.

### What a better family could still be worth, and why it is capped

The original 7 069-line family was compared with a ceiling restricted to
images of that base. The historical calculation found that 94–105% of the
remaining deficit was the overlap two random images have by chance. This is a diagnostic for the original fixed base, and it is
the reason for the `classes` status. It does **not** close a family built by greedy on the
whole line set, which is not limited by image overlap at all: P₄₈ has only 26 208 000 minimal
lines, all of which fit in a GPU, so the Γ₇₂ treatment of
[`../dim73-95-gamma72-caps/scripts/gpu/`](../dim73-95-gamma72-caps/scripts/gpu/) transfers
directly and would give a nearly flat family instead of one decaying from 7 069 to 6 483.

**For the old fixed base, the remaining gain is at most 1.1%.**
This conditional ceiling uses 7 069 lines per class. It is not an upper
bound on larger classes; the October 7 387-line base exceeds that input:

| dim | λ(k) | Σ\|C_i\| now | at 7 069 each | best possible gain | on the claim |
|---|---|---|---|---|---|
| 56 | 120 | 834 576 | 848 280 | +13 704 | +27 408 (+0.05%) |
| 59 | 302 | 2 044 532 | 2 134 838 | +90 306 | +180 612 (+0.32%) |
| 60 | 378 | 2 530 141 | 2 672 082 | +141 941 | +283 882 (+0.49%) |
| 61 | 577 | 3 740 770 | 4 078 813 | +338 043 | +676 086 (+1.13%) |

so no family restricted to the older 7 069-line base moves dimension 61 past 60 574 780, and the dimensions below 56 are
capped under a tenth of a per cent. The quantity that would really move this range is the
BASE CLASS: 7 069 lines is a factor 9.93 over Caro–Wei, against the Leech's 11.6, and a class
of 8 250 would carry every one of these thirteen dimensions up by 17%. That is the class
problem of [`../../closed/class-problem-upper-bound/`](../../closed/class-problem-upper-bound/),
whose published upper bound is 425 lines against 248 found in the Leech — i.e. wide open, and
not a matter of more compute on the family.

## History of the lower bounds in dimensions 49–63

| when | value | who |
|---|---|---|
| 1971 | 52 416 000 | Leech–Sloane [LS71] construct P₄₈, giving τ(48) ≥ 52 416 000; dimensions 49–63 inherit it |
| 1998 | (not competitive) | Edel–Rains–Sloane [EdRS98] — their construction is below 52 416 000 throughout 48–63 and only overtakes at 64 |
| — | 52 416 000 + τ(k) | the direct-sum floor, never written down anywhere as a kissing number because the table stops at 48 |
| this work | **52 430 140 … 67 365 752** | the cap construction with explicit P₄₈ classes |

**Historical literature note from the original release.** No earlier evaluated
claim in these dimensions was known at that time. Later 7 077-line constructions
superseded 49–55D before the October update regained them. The range has simply
not been looked at, because the table stops at 48 and resumes at 64. If someone has
evaluated a construction here, we have not found it.

## Files

```
verify.py                     the exact verification: lattice, class, arithmetic
final.py                      the result table
README.md                     this file
data/p48p_gram.npy            48x48 Gram matrix of P_48p, Catalogue of Lattices
data/p48p_gens.npy            3 published automorphism generators
data/P48p.html, P48q.html     the catalogue pages themselves, for provenance
data/class_p48_7069.npy       the base class: 7069 minimal lines as coefficient vectors
data/class_sizes_p48.npy      the sizes of the 1400 disjoint classes
rows49-63.json                final.py's output, one record per dimension;
                              ceiling.py and scripts/make_readme_table.py both
                              read it, which is how the table above stays
                              identical to what the driver computed
scripts/parse_p48.py          catalogue HTML -> Gram + generators, with exact checks
scripts/regenerate.py         rebuild and re-verify all 1400 disjoint classes: every one
                              a 60-degree code, pairwise disjoint, sizes as claimed
scripts/verify_classes.py     standalone exact checker for saved class files
scripts/verify_extension.py   [superseded] the 130-class extension of a 1170-class
                              base, checked without trusting the builder. The
                              family is 1400 now and regenerate.py verifies all of
                              them; nothing writes the extra_classes_p48.npz this
                              reads, so it exits with a message. Kept for the record
scripts/capdirs.py            the cap-direction line systems for dimensions 57-63
scripts/verify_capdirs.py     those line systems, verified exactly: antipodal, and every
                              pair of distinct lines at 60 degrees or more
scripts/ceiling.py            ceilings restricted to the shipped base classes --
                              what it would give if S reached lambda * alpha exactly
scripts/settle_4963.py        which of the three ingredients are provably at their
                              ceilings, and what that leaves open
scripts/settle_family.py      whether S(lambda) is improvable or the family has already
                              landed where the method must: 82.3% of the ceiling
scripts/make_readme_table.py  regenerates the result table above from the driver's output,
                              so it cannot drift from it again
scripts/gen_lines.py          [pipeline stage] enumerate the 26 208 000 minimal lines
scripts/fast_greedy.py        [pipeline stage] the class search, over that pool
scripts/improve_class.py      [pipeline stage] local search, ~6000 -> 7069 lines
scripts/orbit_classes.py      [pipeline stage] the original disjoint-class run, superseded
                              by regenerate.py
```

**The `scripts/` entries above marked as pipeline stages need the output of the stage before
them**, and those intermediates (the 26 208 000-line pool, the wide search pool) are hundreds
of megabytes and are not shipped. Run bare, they exit with `FileNotFoundError` on a missing
`.npy` — that is expected, not a defect. `verify.py`, `final.py` and `scripts/regenerate.py`
are the three that run from what is here, and between them they cover every claim.


### What is shipped, and what is regenerated

The 1400 classes are 8 023 414 lines of 48 coordinates - 165 MB even compressed - so they
are **not shipped**. They do not need to be: they are a deterministic function of the Gram
matrix, the published generators, the base class and two fixed seeds, all of which are here
and all of which are small. `scripts/regenerate.py` rebuilds them in about eight minutes,
asserts that every image is an isometric image with every row minimal and every pairwise
inner product <= 2, asserts **pairwise disjointness**, checks that no class is empty, and
checks that the size vector it derives is exactly the shipped `data/class_sizes_p48.npy`.

It runs in **two stages**, and the second one matters. |Aut(P_48p)| >= 145 728, so 1400
independent random words collide by the birthday bound about 1400^2/(2*145 728) = 6.7 times,
and a repeated group element reproduces a class already present - which the overlap deletion
then empties. Stage 1 draws 1400 words of length 25 from seed 7, taking for each of the
first 400 the least-overlapping of 16 candidates; stage 2 replaces every class
left below 3000 lines by a fresh image from an independent stream (length 30, seed 99).
Without stage 2 the family has one empty class and 6 918 581 lines.

`scripts/orbit_classes.py` is the original single-stage run, kept for reference; if you want
the classes written out as files, use it and expect 142 MB.

### The provenance chain

```
Catalogue of Lattices page P48p.html
  -> scripts/parse_p48.py     Gram + generators, checked: symmetric, even, det = 1 (Bareiss),
                              A G A^T = G for every generator
  -> scripts/gen_lines.py     the 26 208 000 minimal lines
  -> scripts/fast_greedy.py   a class
  -> scripts/improve_class.py the 7069-line class            <- shipped, verified by verify.py
  -> scripts/regenerate.py    1400 disjoint classes          <- sizes shipped, regenerable
  -> final.py                 the fifteen numbers
```

## A note on the inner-product histogram

`verify.py` reports the class's off-diagonal histogram as
`{0: 11337483, 1: 11979647, 2: 1664716}`, summing to C(7069,2) = 24 981 846. An earlier
version of this package printed 11 341 017 zeros; that count included the diagonal
(7069 entries, halved to 3534) and was wrong by exactly that. The shipped figure is
self-consistent with C(7069,2) and is what `verify.py` recomputes.

## References

See [`../../../common/published.py`](../../../common/published.py).
