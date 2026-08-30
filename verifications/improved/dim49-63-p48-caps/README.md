# Dimensions 49–63: the cap construction over P₄₈

> **Dimensions 62 and 63 are no longer claimed from this construction.** They are claimed
> from the Edel-Rains-Sloane chain, at **71 310 732** and **138 419 844**, in
> [`../dim62-63-ers-chain/`](../dim62-63-ers-chain/). This construction cannot reach those
> figures with *any* class family, and the point is sharper than the totals: even the chain's
> level-0 term alone, 2²⁶ = 67 108 864 and 2²⁷ = 134 217 728, is above its absolute ceiling.
> With every one of the λ(k) classes at the maximum size 7069 the ceiling is 66 075 240 and
> 70 543 480 — short of level 0 by 1 033 624 and 63 674 248.
> [`scripts/ceiling.py`](scripts/ceiling.py) computes the ceiling in all fifteen dimensions and
> checks it against the sign code. Dimensions **49 to 61** are what this package claims, and
> there it wins by a factor 1.7 or more over level 0. The values it reaches in 62 and 63 are
> still computed and shown below, as the weaker bound they are.

Fifteen dimensions, gains from **+14 138** to **+13 847 064**.

```
python verify.py                 # ~1 min: the lattice, the class, the arithmetic
python final.py                  # instant: the result table
python scripts/regenerate.py     # ~15 min: rebuilds and re-verifies all 1400 classes
```

| dim | k | λ(k) | Σ|C_i| | gain | **this work** | previously |
|---|---|---|---|---|---|---|
| 49 | 1 | 1 | 7 069 | +14 140 | **52 430 140** | 52 416 002 |
| 50 | 2 | 3 | 21 206 | +42 418 | **52 458 418** | 52 416 006 |
| 51 | 3 | 6 | 42 402 | +84 816 | **52 500 816** | 52 416 012 |
| 52 | 4 | 12 | 84 746 | +169 516 | **52 585 516** | 52 416 024 |
| 53 | 5 | 20 | 141 091 | +282 222 | **52 698 222** | 52 416 040 |
| 54 | 6 | 36 | 253 417 | +506 906 | **52 922 906** | 52 416 072 |
| 55 | 7 | 63 | 441 802 | +883 730 | **53 299 730** | 52 416 126 |
| 56 | 8 | 120 | 834 576 | +1 669 392 | **54 085 392** | 52 416 240 |
| 57 | 9 | 153 | 1 058 935 | +2 118 176 | **54 534 176** | 52 416 306 |
| 58 | 10 | 255 | 1 738 456 | +3 477 422 | **55 893 422** | 52 416 510 |
| 59 | 11 | 302 | 2 044 532 | +4 089 668 | **56 505 668** | 52 416 604 |
| 60 | 12 | 378 | 2 530 141 | +5 061 123 | **57 477 123** | 52 416 841 |
| **61** | **13** | **577** | **3 740 770** | **+7 482 694** | **59 898 694** | 52 417 154 |
| 62 | 14 | 966 | 5 899 945 | +11 801 822 | 64 217 822 &nbsp;*(superseded by 71 310 732)* | 52 417 932 |
| 63 | 15 | 1 282 | 7 473 594 | +14 949 752 | 67 365 752 &nbsp;*(superseded by 138 419 844)* | 52 418 564 |

so 49–61 clear it by a factor 1.5 or more and 62–63 do not clear it at all, with no class
family that could. The evaluation is
[`../../closed/dim96-ers-takeover/ers_exposure.py`](../../closed/dim96-ers-takeover/ers_exposure.py),
which runs the full chain in every claimed dimension; the two dimensions it displaces are
claimed from that chain instead, in
[`../dim62-63-ers-chain/`](../dim62-63-ers-chain/). The authors themselves say that at
dimension 48 their construction "gives less than half" of 52 416 000, which is what keeps
49–59 safe; what changes at 60 is that A(60,15) ≥ 2²⁵ becomes available, and at 62 that the
level-1 term stops being negligible.

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

`scripts/ceiling.py` measures the distance from the shipped family to the ceiling — every one
of the λ(k) classes at the maximum 7 069 — and finds that 94–105% of it is the overlap two
random automorphism images have BY CHANCE. That argument closes the *images* route, and it is
the reason for the `classes` status. It does **not** close a family built by greedy on the
whole line set, which is not limited by image overlap at all: P₄₈ has only 26 208 000 minimal
lines, all of which fit in a GPU, so the Γ₇₂ treatment of
[`../dim73-95-gamma72-caps/scripts/gpu/`](../dim73-95-gamma72-caps/scripts/gpu/) transfers
directly and would give a nearly flat family instead of one decaying from 7 069 to 6 483.

**It is worth at most 1.1%.** The ceiling is what caps it, and the ceiling is exact:

| dim | λ(k) | Σ\|C_i\| now | at 7 069 each | best possible gain | on the claim |
|---|---|---|---|---|---|
| 56 | 120 | 834 576 | 848 280 | +13 704 | +27 408 (+0.05%) |
| 59 | 302 | 2 044 532 | 2 134 838 | +90 306 | +180 612 (+0.32%) |
| 60 | 378 | 2 530 141 | 2 672 082 | +141 941 | +283 882 (+0.49%) |
| 61 | 577 | 3 740 770 | 4 078 813 | +338 043 | +676 086 (+1.13%) |

so no family whatsoever moves dimension 61 past 60 574 780, and the dimensions below 56 are
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

**No prior claim in any of these fifteen dimensions is known to us.** The range has simply
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
scripts/ceiling.py            the ceiling of this construction, dimension by dimension --
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
