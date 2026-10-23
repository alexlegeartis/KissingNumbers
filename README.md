# Kissing numbers: constructions and verification in dimensions 18 through 96

The kissing number τ(n) is the largest number of unit spheres that can touch a central unit
sphere in ℝⁿ; equivalently, the largest set of unit vectors with pairwise inner products at
most 1/2. This repository documents constructions in **58 dimensions**, each with a package
that states its verification scope, together with the negative results that guide further
searches. Some of our earlier records have since been superseded. Their constructions,
values and verification packages remain here as part of the research record.

> [!IMPORTANT]
> **October 2026 — new SOTAs in dimensions 41–46.** The October constructions
> join our earlier results in the table below. The reported 43D total of
> **2 911 990** awaits its certificate; the recovered **2 878 155**-point seed
> has passed fresh population and exact angle verification.

```bash
python audit.py --repository-only  # all 58 claims; check this repository without its surrounding workspace
```

`update.py` is the maintainer's entry point rather than the reviewer's: it runs every package
driver, then `audit.py`, then `audit.py --write-results`, then each package's table generator,
then the paper's tables and its checkers, then `run_all.py` -- in that order, so a number can
only travel outward from the driver that computed it. Run it after any change (`--check` to
report without rewriting anything). Nothing in it decides a value; steps 1-5 only move numbers
a driver already produced.

`audit.py` recomputes what is cheap to recompute — dimensions 47 and 71 from Venkov's
theorem, dimensions 49–95 from the shipped class-size tables — and checks them against
the pinned comparison data in `common/published.py`, for monotonicity and against records
in higher dimensions. These are checks against a recorded snapshot; they do not establish
that every packaged value remains the best known today. It takes a few seconds.

## Results and construction history

The comparison column records the value used when each result was obtained. Factors
refer to that comparison. **Superseded** marks results subsequently improved by others;
our original values and construction descriptions are retained. See
[record status](#record-status-and-retained-results) for the later values and sources.
The status notes use [Cohn's live leaderboard](https://cohn.mit.edu/kissing-numbers/),
retrieved **23 October 2026**. Its 41–46D entries agree with the October comparison
values below; its 18D entry is our 8 358-point result, credited to this repository.

| dim | recorded comparison | **this work**     | factor   | construction and record status |
| --- | -------------------- | ----------------- | -------- | ---------------------------------------------------------------- |
| 96  | 6 218 372 160        | **12 886 999 232**| **2.07** | no table reaches it, and Edel–Rains–Sloane was never evaluated there |
| 71  | 331 737 984          | **2 603 658 750** | **7.85** | the tables have no entry between 65 and 71                       |
| 70  | 331 737 984          | **1 249 778 250** | **3.77** | the same                                                         |
| 69  | 331 737 984          | **627 822 180**   | **1.89** | the same, with the Gram forced by Hermite                        |
| 68  | 331 737 984          | **361 275 480**   | 1.09     | the same, at k = 4                                               |
| 63  | 52 418 564           | **138 419 844**   | **2.64** | the tables stop at 48 and resume at 64                           |
| 62  | 52 417 932           | **71 310 732**    | 1.36     | the same                                                         |
| 46 | 12 309 600 | **13 387 797** | 1.09 | October: independently normalized P48n families with six-dimensional heads, plus 27 pure-head points |
| 45 | 7 380 090 | **7 650 601** | 1.04 | October: joint rearrangement of body families and five-dimensional heads, plus 16 pure-head points |
| 44 | 4 153 600 | **4 939 281** | 1.19 | October: 1 034 complete body families with rational four-dimensional heads |
| 43 | 2 553 792 | **2 911 990*** | 1.14 | October reported total: a switched triangular face, repacked auxiliary families and seven added points; certificate pending |
| 42 | 1 792 386 | **1 991 445** | 1.11 | October: jointly replaced auxiliary layers with two-dimensional heads, plus a pure-head hexagon |
| 41 | 1 324 472 | **1 460 036** | 1.10 | October: body families with one-dimensional heads, plus two pure-head poles |
| 38  | 566 652              | **591 612**       | 1.04     | **Superseded.** The Leech cap construction at codimension 14 |
| 31  | 238 350              | **238 662**       | 1.00     | **Superseded.** The norm-8 frame layer at height √2, 672 points deleting no equator point, carried by 42 pairwise disjoint type-B classes holding 10 328 of 10 416 owner lines, with an axis of 118 of 126 |
| 30  | 220 440              | **221 012**       | 1.00     | **Superseded.** The dimension-28 norm-8 frame layer at height √2: a whole Leech frame of heads carries the whole cross-polytope of ℝᵏ, 96k points, deleting nothing |
| 29  | 209 496              | **209 968**       | 1.00     | **Superseded.** The same frame layer one dimension up, on the whole cross-polytope of ℝ⁵ — 480 points, deleting nothing |
| 28  | 204 520              | **204 896**       | 1.00     | **Superseded.** The head of a deletion-free layer need not have norm 6: at height √2 a norm-8 head clears the cap threshold, and 48 of them are a whole Leech frame |
| 27  | 200 044              | **201 566**       | 1.01     | **Superseded.** The coset triangle of three norm-6 vectors on all four triangles of directions, with the class heads, the axial normals that remove nothing and the second cap layer all solved together rather than in sequence (first layer B. Lindow, second layer here) |
| 26  | 198 550              | **199 806**       | 1.01     | **Superseded.** The same triangle on both triangles of the hexagon, the second side the involution image of the first, plus a jointly-solved layer of free heads that remove nothing |
| 25  | 197 056              | **197 579**       | 1.00     | **Superseded.** 1016 heads in the lens of a minimal vector, each removing only its owner, plus one non-lattice equator point |
| 18  | 7 654                | **8 358**         | 1.09     | the base of Cohn–Li's hexagon had a free combinatorial parameter nobody had varied: which coset of RM(1,4) the three six-set families lie in.  Three bent cosets with pairwise-bent sums that do not sum to zero give the twelve tier-B slots three different 512-word spaces instead of one, and tier B rises from 256 to 960 |


The full table of all 58 registered constructions is [RESULTS.md](RESULTS.md).
Its comparisons use the repository's pinned data and include results that have since
been superseded. The earlier manuscript is in [paper/](paper/); it remains its earlier
snapshot. See [The write-up](#the-write-up).

The [P48n family/head package](verifications/improved/dim41-46-p48-heads/) supplies
the October constructions. The five registered totals in 41D, 42D and 44–46D pass
integer angle and exact rank checks using inherited body populations. A fresh replay
also passed for the older standalone 41D/42D/44D population certificates.
**\*43D:** 2 911 990 is a reported working total, excluded from the registry until
its certificate is supplied; the older 2 878 155-point seed is fully certified under
the stated P48n assumptions. The latest progress report supersedes older archive
totals. No October improvement is claimed in 47D or 48D.

## Layout

```
audit.py                  every claim, checked for mutual and external consistency
run_all.py                every verification script, with a single verdict
RESULTS.md                the full table (generated by audit.py --write-results)

common/                   shared modules and the two proofs
  kpoint_lp.py            the k-point moment LP with exact rational dual certificates
  theta.py                extremal theta series; the Venkov one-point distribution
  published.py            pinned published comparison data, with provenance
  leech.py golay.py m24.py
  PROOF-kpoint.md         the cross-section method, as numbered lemmas with proofs

verifications/
  improved/               57 of the 58 claims, fourteen packages, one per construction
                          idea; dimension 96 is the exception and sits in closed/
  recovered/              values already published, re-derived here and claimed nowhere
  superseded/             claims this project made and then lost, and why
  closed/                 mechanisms pushed to their exact ceiling: what NOT to retry

KNOWLEDGE.md              the full working record, 158 sections, including everything
                          that failed -- with a preamble on how to read it, which
                          sections supersede which, and where each script now lives
```



### `verifications/improved/` — fourteen ideas, 57 dimensions

Every claim but one is in this tier. The exception is dimension 96, whose package is
`closed/dim96-ers-takeover/`, because what it establishes is that the cap construction
cannot reach 96 — the claim there is Edel–Rains–Sloane’s.

These directory names preserve the classification at packaging time. A construction
in `improved/` may now be superseded by a later publication; its package is retained.


| package                                                                                      | dims   | the idea                                                                                                                                                                                     |
| -------------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[dim25-lens-heads](verifications/improved/dim25-lens-heads/)`                               | 25     | cap heads of squared length 3 in the lens of a minimal vector remove exactly their owner and never share a removal, so τ = 196560 + H + 2 + \|E\| with H = 1016 and one non-lattice equator point |
| `[dim26-27-iota-triangles](verifications/improved/dim26-27-iota-triangles/)`                 | 26, 27 | the axis forces the heads onto squared length 8/3 and onto fixed triangles of directions, so the layer is two or four head sets at cosine ≤ 1/4 cross-constrained at cosine 1/2; the coset triangle of three norm-6 vectors summing to zero is one such set, and the involution ι(u) = −v − u on each class puts S and ι(S) on two triangles at once |
| `[dim28-norm8-frame-layer](verifications/improved/dim28-norm8-frame-layer/)`                 | 28     | the head of a deletion-free layer need not have norm 6: at height √2 a NORM-8 head clears the cap threshold where a norm-6 head does not, and the 24 vectors 8eᵢ are a Leech frame, each carrying all eight directions ±e_k of ℝ⁴ — 384 points for 8 axis points |
| `[dim29-30-frame-layer](verifications/improved/dim29-30-frame-layer/)`                           | 29, 30 | the dimension-28 norm-8 frame layer is not about dimension 28: at height `√2` with `x = v/2`, `|v|² = 8`, the layer deletes nothing, one head carries the whole cross-polytope of `ℝᵏ` and all 48 vectors of a Leech frame carry every direction, so it is `96k` points in every dimension. The price is that all `⌊τ(k)/3⌋` owner classes must be type B for one frame — and in the coordinates where the frame is `8eᵢ` the type-B lines are the octad vectors of the Golay code, which the monomial group `2¹²:M₂₄` permutes, so the filter disappears and what is left is a packing solved by coordinate descent on the Golay sign words |
| `[dim31-frame-layer](verifications/improved/dim31-frame-layer/)`                             | 31     | the norm-8 frame layer at height √2 carried to k = 7: all 48 vectors of a Leech frame on the whole cross-polytope of ℝ⁷, 672 points, deleting no equator point.  It needs 42 pairwise disjoint type-B classes holding at least 10 252 of the 10 416 lines 42 full classes would give; the published attempt reached 10 183 and this one carries 10 328, found by screening candidates on the free slots their octads still have rather than by scanning uniformly — a sampled clique search needs pair density p ≥ 0.5097 to reach 42 at all, and flattening the octad coverage takes p from 0.4984 to 0.5467.  The axis is E7 rotated by 45°, 45°, 90° in the coordinate planes (0,3), (2,4), (5,6), keeping 118 of 126 against the 110 the published axis kept |
| `[dim38-leech-large-codimension](verifications/improved/dim38-leech-large-codimension/)`     | 38     | the same construction at codimension 14, where the binding constraint flips and the whole Leech shell partitions into 644 classes                                                            |
| `[dim39-ers-constant-weight](verifications/improved/dim39-ers-constant-weight/)`             | 39     | Edel–Rains–Sloane with the 2026 constant-weight codes, at n₀ = n rather than n₀ = 32                                                                                                         |
| [`dim41-46-p48-heads`](verifications/improved/dim41-46-p48-heads/) | 41, 42, 44–46 | independently normalized P48n body families with rational heads and pure-head extensions; inherited populations; latest 43D certificate missing |
| `[dim49-63-p48-caps](verifications/improved/dim49-63-p48-caps/)`                             | 49–61  | the cap construction over P₄₈, with explicit classes of 7069 lines where Caro–Wei guarantees 712                                                                                             |
| `[dim62-63-ers-chain](verifications/improved/dim62-63-ers-chain/)`                         | 62, 63 | the Edel–Rains–Sloane chain (n, 15, 2), levels 1–2 built here. Its level-0 term alone already beats the cap construction's **ceiling** in these two dimensions, so no class family could have reached them |
| `[dim68-69-gamma72-cross-sections](verifications/improved/dim68-69-gamma72-cross-sections/)` | 68, 69 | realizability of a k-point Gram settled by **exhibiting** the k-tuple in Γ₇₂ rather than through the LP's own witness lemma, which unlocks k = 3 and k = 4 |
| `[dim70-71-gamma72-cross-sections](verifications/improved/dim70-71-gamma72-cross-sections/)` | 70, 71 | the same k-point method on Nebe's Γ₇₂                                                                                                                                                        |
| `[dim73-95-gamma72-caps](verifications/improved/dim73-95-gamma72-caps/)`                     | 73–95  | the cap construction over Γ₇₂, whose class threshold γ = 1/4 is exactly what permits zero-sum triples                                                                                        |


Each package has a README with the idea, a description of every file, and **the history of
the lower bound in its dimensions** — who held it, when, and by what construction.

### `verifications/recovered/` — nothing claimed


| package                                                                               | dims   | what it is                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[dim46-47-p48-cross-sections](verifications/recovered/dim46-47-p48-cross-sections/)` | 46, 47 | Venkov's theorem applied to **products** of inner products with several fixed minimal vectors. Both values were already published; re-deriving them from a method that knows nothing of either is the sharpest check the cross-section machinery has, which is why the package is kept |




## Record status and retained results

**Later improvements.** Our earlier values remain valid lower bounds and stay in
this repository. The later values below are from
[Cohn's leaderboard](https://cohn.mit.edu/kissing-numbers/) as retrieved on
**23 October 2026**, and, for 49–55D, from Yang et al.,
[*Large-Scale Autonomous Discovery of Kissing Number Constructions*, Table 1](https://arxiv.org/html/2609.35051v1#S1.T1)
(28 September 2026). Cohn's table covers dimensions through 48; the higher-dimensional
comparison is given separately by its source.

| dim | our retained result | later result | source |
|---|---:|---:|---|
| 25 | 197 579 | 197 583 | Cohn, 23 October |
| 26 | 199 806 | 199 868 | Cohn, 23 October |
| 27 | 201 566 | 201 589 | Cohn, 23 October |
| 28 | 204 896 | 205 172 | Cohn, 23 October |
| 29 | 209 968 | 210 495 | Cohn, 23 October |
| 30 | 221 012 | 222 073 | Cohn, 23 October |
| 31 | 238 662 | 241 079 | Cohn, 23 October |
| 38 | 591 612 | 605 322 | Cohn, 23 October |
| 39 | 756 116 | 768 916 | Cohn, 23 October |
| 49 | 52 430 140 | 52 430 156 | Yang et al., Table 1 |
| 50 | 52 458 418 | 52 458 468 | Yang et al., Table 1 |
| 51 | 52 500 816 | 52 500 930 | Yang et al., Table 1 |
| 52 | 52 585 516 | 52 585 772 | Yang et al., Table 1 |
| 53 | 52 698 222 | 52 698 696 | Yang et al., Table 1 |
| 54 | 52 922 906 | 52 923 774 | Yang et al., Table 1 |
| 55 | 53 299 730 | 53 301 140 | Yang et al., Table 1 |

These dated sources determine the supersession notes; a new global record count
is not inferred from the repository's older comparison cache.
The October 41–46D constructions are presented alongside the earlier results above,
with the outstanding 43D certificate obligation stated explicitly.

**Historical recovered cross-sections in dimensions 46 and 47.** τ(47) ≥ 23 766 960 is due to
Boyvalenkov–Cherkashin, *Results in Mathematics* **80** (2025), Paper No. 3, and is
equation (3) of the preprint [arXiv:2312.05121](https://arxiv.org/abs/2312.05121)
(December 2023), which observes *in the same paragraph* that those A₀ vectors form a
47-dimensional kissing configuration.
τ(46) ≥ 12 309 600 is Table 3 of Sun–Wang,
[arXiv:2607.20359v3](https://arxiv.org/abs/2607.20359) (18 August 2026), and its ingredients
were in Ozeki's 2016 Siegel theta tables a decade earlier, though not stated as a kissing
number. Both were obtained here independently, and both are **recovered rather than
claimed**: they are not in `RESULTS.md`, not counted among the 58, and not presented as
results in the paper. They are kept, in
`[verifications/recovered/](verifications/recovered/)`, because reproducing a published
count from a method built without reference to it is the strongest validation available.
Cohn's live table now includes both recovered values. The pinned August comparison
cache still records the older 9 741 412 and 5 318 060 entries.

**Withdrawn earlier constructions.** The old lattice cross-section values in
dimensions 45 and 44 were beaten by Sun–Wang. The new family/head constructions
are separate, larger results. Both of this
project's values lie inside the LP brackets it computed, which is a sharp cross-check on the
machinery — see
`[verifications/superseded/dim44-45-p48-cross-sections/](verifications/superseded/dim44-45-p48-cross-sections/)`.

**Historical calibration.** The 19D construction studied here has exact ceiling
11 948, matching Ho's record at that time. Cohn's live table now lists 12 460;
the ceiling remains a result about that construction family. Dimension 12's record
likewise landed elsewhere while this work was in progress.

**Original comparison scope.** Dimensions 18, 25–31, 38, 39, 41–42, 44–46, 49–63, 68–71, 73–95.
These are the registered construction results compared with their recorded baselines;
dimension 96 is the evaluated Edel–Rains–Sloane chain. The list includes the superseded
results above. An earlier improvement remains a constructive result after a larger
configuration is found. The literature status of the 70D and 71D consequences of the
one-point distribution of Γ₇₂'s minimal shell still merits specialist review.

## The construction ideas

**Independently normalized body families and heads (41–46).** The October constructions
use a rank-40 complement of sqrt(3) E8 in P48n, with heads in one to six dimensions.
Each family is normalized separately, so families can have different original lengths.
Joint class replacements and head motions admit additional complete or partial families;
zero-body points in the existing head span give the latest 41D, 42D, 45D and 46D
extensions. The 43D working construction additionally uses seven nonzero-body points
and awaits its coordinate certificate. Exact checks, population assumptions and fresh
replay reports are documented in the
[verification package](verifications/improved/dim41-46-p48-heads/).

**Cross-sections of extremal lattices (46, 47, 68–71).** For an extremal even unimodular
lattice, Venkov's theorem makes the minimal shell a spherical 11-design, so every polynomial
of degree ≤ 11 integrates exactly. Applied to *products* of inner products with several fixed
minimal vectors this constrains the **k-point** distribution, and minimising the all-zero
cell by linear programming — with an exact rational dual certificate, plus translation and
lattice-range constraints that the design property cannot see — bounds the cross-section.
Dimension 46 comes out **exact** (the LP minimum and maximum coincide). The same code
reproduces every cross-section whose answer is known: E₇ = 126, Λ₂₃ = 93 150, D₆ = 60,
E₆ = 72, 43 164, 44 550, 49 896, 27 720. Mathematics in
`[common/PROOF-kpoint.md](common/PROOF-kpoint.md)`.

**The cap construction (25–27, 29–31, 38, 49–63, 73–95).** In ℝⁿ ⊕ ℝᵏ, put a lattice shell on the
equator, delete a few lines, and put "caps" above the deleted directions. The optimum uses
**zero-sum triples** of cap directions rather than antipodal pairs, which needs class
threshold γ ≤ 1/4 and cap level t = 2/3; three is the maximum possible, because *n* points
pairwise at cosine ≤ −1/2 force *n* ≤ 3. Run over the Leech with its 248-line class this
reproduces Cohn's table in dimensions 26, 28, 29, 30, 31 **exactly** and beats it in 25 and
27 — a strong check on the model, and it is run as such by
`[verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py](verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py)`.

**Explicit classes.** The gain is proportional to the size of a *class* — a set of lines with
pairwise |⟨u,u′⟩| ≤ μ/4. Caro–Wei bounds it with no coordinates at all, but real search does
about ten times better in every lattice tried: 221 → **2118** for Γ₇₂, 712 → **7069** for
P₄₈, and 32 000 pairwise-disjoint classes totalling 47 150 230 lines.

**The Edel–Rains–Sloane chain, evaluated where its authors did not (39, 62, 63).** Their count
is a sum over a chain of support sizes, `N(n) = Σ_ν A(n, n_ν, n_ν)·A(n_ν, ⌈n_ν/4⌉)`, and they
evaluated it only at 32, 36, 40, 44, 64, 80 and 128. At `n = 39` with the 2026 constant-weight
tables it gives 756 116; at `n = 62` and `n = 63` the chain `(n, 15, 2)` gives 71 310 732 and
138 419 844. Its levels above zero are **built** here rather than cited — a cyclic `[15,6,8]₄`
through the grid map, and a `[15,10,4]` — and the whole chain is checked in coordinates at
`n = 8` and `n = 16`, where it reproduces τ(8) = 240 (the E₈ root system) and τ(16) = 4320 (the
Barnes–Wall value, and the record) exactly.

## Standards

**Two suites exist specifically to attack the mathematics rather than the bookkeeping**, and
they are the ones to run if you doubt any of this:

```bash
python common/validate_lp.py    # 15 independently known truths, each inside its LP bracket
python common/capalgebra.py     # the cap construction's algebra and count, from scratch
```

`validate_lp.py` exists because the k-point LP produces *lower* bounds: one that is too small
is merely weak, but one that is too **large** is wrong, and it becomes too large the moment
any single constraint is not valid. It checks LP-min ≤ truth ≤ LP-max on fifteen
configurations whose size is known independently — from four lattices and three separate
sources of truth. If any constraint were too strong, some truth would fall below its own
minimum. None does. `capalgebra.py` re-derives the cap construction's seven inequalities and
its count in exact rational arithmetic, independently of every script that uses them, and
reproduces Cohn's table in dimensions 26, 28–31 from the count alone.

No floating-point arithmetic decides anything. Every linear-programming bound carries an
exact rational dual certificate, so an imprecise solver can only produce a *weaker* bound,
never an invalid one. Every Gram matrix used is proved realizable — a bound for a Gram that
does not occur is worth nothing, and dimension 69 was withdrawn once on exactly this ground
before the triple was exhibited in coordinates. Every construction is verified in explicit coordinates
where the configuration is small enough to write down, and by exact integer arithmetic on the
Gram where it is not, with outstanding certificate gaps identified explicitly.
Published comparisons are dated: **Cohn's table moves**, and two entries (dimensions
12 and 19) changed while this work was in progress. The supersession table above records
the later September publication separately from the pinned computational baselines.

Negative results are recorded as carefully as positive ones, in
`[verifications/closed/](verifications/closed/)`, because they are what stops the next
attempt from wasting time. Among them: dimensions 17–24 are a single construction, the layer
identity τ(Λ₁₆₊ₖ) = 4320 + 513·N₂(Lₖ) + 32·N₄(Lₖ), verified as a structure inside the Leech
lattice and with every component at a stated ceiling but one; the Cohn–Li mechanism is
provably exhausted in dimensions 20 and 21 and closed at exactly Ho's value in 19; Λ₂₁, Λ₂₂ and Λ₂₃ are maximal
spherical codes with exact min-max cosines √(8/29), √(3/11), √(4/15); every published record
in dimensions 9–19 is maximal; the antipode construction is a max-weight clique problem whose
arms lie in a 30° cap, which caps it structurally; and the maximum-class problem survives
every two-point method, the three-point SDP and the subconstituent split at 4680/11 = 425.45.

## Errata found during this work

Recorded because they are the kind of thing that survives a careful reading:

- **The dimension-40, -41 and -42 verdicts went stale**, and the failure mode is the one this
project keeps meeting: `verifications/closed/dim32-44-ers-audit/` concluded in August 2026 that
dimensions 32–44 sit *exactly* at the Edel–Rains–Sloane value and improve none. On 2026-09-14 a
re-read of the literature found I. Dorofeev, X. Sun and C. Wang, *Optimal Extensions of
Cross-Sections: Sphere Packings in Dimensions 38 to 43*, arXiv:2607.20359v4 (29 Aug 2026),
which beats it in three of them with **lattice** cross-sections of an extremal 48-dimensional
lattice — 1 092 000 / 1 324 472 / 1 792 386 in dimensions 40, 41, 42 — and raises dimension 45
to 7 379 838. Nothing here was claimed in those dimensions, so no claim is withdrawn, but the
package's stated verdict was wrong and now carries a correction banner. **Dimensions 32–37 are
unaffected, and now for a measured reason**: calibrating their own mechanism on their own three
results (it reproduces dimension 40 to 3 %) and applying their determinant floor caps
dimension 32 at about 205 000 — a factor 1.69 below its record. The lesson is the standing one:
a table is a function of the literature, and the literature moves. See `KNOWLEDGE.md` §131.

- **The k-point LP was missing its mixed moments** for a day — it imposed only exponent
vectors with every entry even, 56 equations at k = 3 instead of 161. The symptom was
visible in the output and unnoticed: two *equivalent* Gram matrices gave different answers.
Dimensions 46, 47 and 71 were unaffected (k = 1 has no mixed moments, and 46 was already
exact); dimensions 44, 45 and **70** all moved. Cheap decisive test: run the LP on a pair
of equivalent Grams.
- **Dimension 69's row was stale** in an earlier results table, quoting an LP value for a
Gram whose realizability is unproved. It is now in `closed/`.
- **The calibration script never calibrated anything**: it crashed at k = 6 (E₆ and E₇ are
returned as subsets of the E₈ roots, so they live in ℝ⁸) and used small greedily-found
classes rather than the real 248-line one. Rewritten; it now reproduces five table entries
exactly.
- **The August comparison table lagged four dimensions** — 32, 33, 34 and 37 — after
Echols' 2026 constant-weight improvements appeared on Brouwer's page. In particular,
496 232 improved the table's 494 312 in dimension 37. Yang et al.'s September
Table 1 subsequently reports 498 024, and Cohn's live table now lists 500 008.
These are comparison updates, not constructions claimed here.

- **The Edel–Rains–Sloane floor was taken at level 0 only**, twice over. First, dimensions
62 and 63 were claimed from the cap construction over P₄₈ at 64 217 822 and 67 365 752 while a
*single sign code* gives 2²⁶ and 2²⁷; `audit.py` did not catch it because it read its floors
from the same stale table, and it now checks the level-0 floor as its own step 4b. Then the
same mistake one level up: level 0 is not the whole construction, and the full chain
`(n, 15, 2)` gives **71 310 732** and **138 419 844**, so those two dimensions moved again.
The full chain is now evaluated in *every* claimed dimension by
`closed/dim96-ers-takeover/ers_exposure.py`, which also ranks each claim by how far the best
published constant-weight code would have to be beaten to overturn it — dimension 68 is the
most exposed, at a factor 1.47.
- **The GPU class certificate is not reproducible across software stacks.** `run.py` claimed
to be “bit-identical on any machine”; measured, the 31 000-class certificate made under one
Colab image did not regenerate under a later one, while two builds inside one process agreed
bit for bit. Every random *decision* is numpy's, but the pool is *assembled* by torch
primitives whose tie-breaking is unspecified, and the pool's row order is what the greedy
window is drawn from. This was not merely untidy: everything downstream is indexed by pool
POSITION, so a `--resume` against a different pool would silently reuse lines, and the
invariant `|alive| + Σsizes = |pool|` counts rows and would not notice. A resume done that
way had to be discarded and the family rebuilt in one session. `gamma72.pool_digest` now
fingerprints the pool, the state carries it, and a resume without a match is refused.

- **A plausible mechanism ended an investigation, and cost 24 spheres in dimension 38.**
The axis layer needs a 60° code of ℝ¹⁴ lying 30° from every cap direction; a rotated copy of
the cap directions kept 1908 of 1932, and the README explained the missing 24 away: the caps
cover 1.24% of S¹³, so “a random rotation loses about 24 points and best-of-N would recover
only a handful”. The first half is a measurement and is right — 200 random rotations lose 56
on average, sd 11. The second half was not measured and was false: SO(14) has 91 dimensions
and the softened objective is differentiable, so **descent** walks 1908 → 1932 in two hundred
steps. Dimension 38 is now 591 612, exactly its family ceiling 3K(24) + K(14). Nothing in a
repository can distinguish a measured statement from a plausible one when both are prose;
`scripts/axis_rotation.py` now carries the search so the claim can be re-derived rather than
argued.

- **A sufficient condition became the search space, and cost 221 025 spheres.** A pole in the
Γ₇₂ cap construction must be a 60° code and lie 30° from every cap direction. The rule used
looked for one **inside the same lattice**, where a shell argument makes the 30° condition
automatic — for a of norm M and w of norm m in one lattice, |a−w|² and |a+w|² are lattice
norms, so |⟨a,w⟩| ≤ M/2, and M ≤ 3m settles it. Correct, exact and beautiful, and it answered
“is this pole legal?” so cleanly that “which poles are there?” stopped being asked: at k = 23
it found 248 poles against τ(23) = 93 150. Nothing forces a pole into the lattice. A rotated
copy of the cap directions is a 60° code for free, a single 30° cap covers 2.3×10⁻⁸ of S²²
so all 93 150 together cover at most 2.1×10⁻³, and the rotation shipped loses 76 of the
93 150 — so the layer is 93 074. Twenty
dimensions carry a layer above the shell rule — 75 and 77–95 —
and four of them have since gone further still (KNOWLEDGE.md §105).  §106 then went past the layers entirely: the cap DIRECTIONS at k = 12 and 13 were not the largest available either, and replacing them with the Kappa sections K₁₂ and K₁₃ — plus two partitions made perfect by an order-3 isometry — moved six dimensions by 549 638, against §105's 1040.  `verifications/improved/dim73-95-gamma72-caps/scripts/axis_rotate.py`
derives the rotations and `scripts/verify_poles.py` checks every pair in exact integers.

- **A guard pinned the wrong number in place.** `factcheck.py` tested that the paper contained
the literal `2.3\times10^{-3}`. It did, and the check passed for a day while that number named
the wrong quantity — a single 30° cap on S²² covers 2.284×10⁻⁸, and 2.1×10⁻³ is the union of
93 150 of them. Worse, the guard made the error load-bearing: correcting the paper would have
failed the build. A check that tests for a literal cannot distinguish a right number from a
wrong one, and it converts the first into a requirement. The replacement derives the cap
measure from ∫sin²¹, derives the union bound and the expected loss from that, reads the
shipped layer from its own `.npz`, and checks the paper against all three — and was
perturbation-tested one figure at a time.

- **Seven axis layers had the poles in the wrong space, and every check passed.** The
construction puts the cap directions and the poles in one ℝ^k; the rotation must therefore be
an isometry of ℝ^k. The one built was an isometry of the coordinates the direction set is
*written* in, and Λ₂₁'s cross-section is 27 720 vectors written in Λ₂₄'s 24 coordinates
spanning 21. So the poles left the axis space, and dimensions 85, 87, 89, 90, 93, 94 and 95
needed up to ℝ⁹⁶ instead of ℝ^(72+k). Every inner product was correct, which is why nothing
saw it: the exact cap test, the isometry NᵀCN = D²C, the pole–pole inheritance and the
independent `Fraction` check are all statements about the Gram matrix, and a rotation of the
ambient coordinates preserves every entry of it. **A verification can be exhaustive over the
objects it compares and silent about what those objects are.** Rank is cheap, and no
inner-product check of any depth could have substituted for it.
`scripts/axis_block.py` derives the replacements from a generator built out of the direction
set, which fixes the orthogonal complement pointwise; `scripts/verify_poles.py` and audit
check 5o both check the space now.


- **Two files stated τ(20) two ways.** `axis_rotate.py` and `verify_poles.py` each carried a
private `TAU = {...}` with τ(20) = 17400 and τ(21) = 27720, the Λ₂₀ and Λ₂₁ values, both 2048
below the best known — the best configurations in 20 and 21 are not lattices. Nothing computed
used them; they set a printed column, so nothing failed and nothing could. Both now import
`common/published.py`, and audit check **5n** evaluates every integer-valued τ table anywhere
in the tree against it. Retyping the two numbers would have fixed the instance; the import
fixes the class.




## Reproducing everything

```bash
python run_all.py          # 68 scripts, about 40 minutes measured, one verdict
python run_all.py --changed   # skip jobs whose inputs have not moved since they passed
python run_all.py --only dim31 # just the jobs whose label or command matches
python run_all.py --full   # 86 scripts, budget about 9 hours: adds the all-pairs sweep, the
                           # negative controls, the class regeneration and the LP brackets
python run_all.py --list   # what would run, and roughly how long each takes
python run_all.py --only p48-heads  # October constructions: exact consistency replay
python update.py --skip-paper --repository-only  # maintain constructions without changing the manuscript
```

Or individually, each from its own directory:


| script                                                                        | what it establishes                                            | time   |
| ----------------------------------------------------------------------------- | -------------------------------------------------------------- | ------ |
| `audit.py`                                                                    | all 58 claims, mutual and external consistency                 | 30 s   |
| `common/theta.py`                                                             | E₇ = 126, Λ₂₃ = 93 150, dim 47, dim 71                         | 1 s    |
| `common/kpoint_lp.py`                                                         | the LP, validated on the Leech and P₄₈                         | 10 s   |
| `…/dim25-lens-heads/verify.py`                                                | τ(25) ≥ 197 579, every pair in exact arithmetic               | 4 min  |
| `…/dim26-27-iota-triangles/verify26.py`                                      | τ(26) ≥ 199 806, every pair in exact arithmetic               | 1 min  |
| `…/dim26-27-iota-triangles/verify27.py`                                      | τ(27) ≥ 201 509, every pair in exact arithmetic               | 2 min  |
| `…/superseded/dim27-triple-partition/scripts/verify_configuration.py`         | τ(27) ≥ 200 540, the superseded claim, from its coordinate file | 18 s   |
| `…/dim28-norm8-frame-layer/verify28.py`                                      | τ(28) ≥ 204 896, the norm-8 frame layer                       | 2 s    |
| `…/dim29-30-frame-layer/verify.py 29`                                        | τ(29) ≥ 209 968, the same layer at k = 5                      | 2 min  |
| `…/dim29-30-frame-layer/verify.py 30`                                        | τ(30) ≥ 221 012, and at k = 6                                 | 2 min  |
| `…/dim31-frame-layer/verify.py 31`                                           | τ(31) ≥ 238 662, the norm-8 frame layer at k = 7              | 7 min  |
| `…/dim28-norm8-frame-layer/ceilings.py`                                      | dim 28: the direction weight and the axis are at their ceilings | 3 s  |
| `…/dim29-30-frame-layer/ceilings.py 29` and ` 30`                            | the same report at k = 5 and k = 6                            | 5 s    |
| `…/dim31-frame-layer/ceilings.py 31`                                         | the same report at k = 7, 1120 vertices enumerated            | 20 s   |
| `…/dim38-leech-large-codimension/scripts/verify.py`                           | τ(38) ≥ 591 612                                                | 10 s   |
| `…/dim39-ers-constant-weight/scripts/verify.py`                               | τ(39) ≥ 756 116 (and τ(38) ≥ 570 236)                          | 10 s   |
| `.../dim41-46-p48-heads/verify.py` | current 41, 42, 44-46 totals and the older 43D seed; inherited populations | 5 s |
| `…/dim46-47-p48-cross-sections/derive.py`                                     | recovers 12 309 600 and 23 766 960, claiming neither           | 6 s    |
| `…/dim49-63-p48-caps/verify.py` + `final.py`                                  |  the P₄₈ lattice, class, and fifteen numbers                    | 3 s    |
| `…/dim49-63-p48-caps/scripts/regenerate.py`                                   | rebuilds and re-verifies all 1400 disjoint classes             | 15 min |
| `…/dim49-63-p48-caps/scripts/ceiling.py`                                      | the construction's **absolute** ceiling, against level 0        | 1 s    |
| `…/dim62-63-ers-chain/scripts/verify.py`                                     | τ(62) ≥ 67 108 864, τ(63) ≥ 134 217 728, mechanism verified      | 1 s    |
| `…/dim62-63-ers-chain/scripts/verify_chain.py`                               | τ(62) ≥ 71 310 732, τ(63) ≥ 138 419 844, levels 1–2 built here | 20 s   |
| `…/dim68-69-gamma72-cross-sections/scripts/leech_truth.py`                     | the k-point LP is **exact** on the Leech through k = 4          | 2 s    |
| `…/dim68-69-gamma72-cross-sections/scripts/recheck_cert.py`                    | both certificates re-derived by a second arithmetic path        | 40 s   |
| `…/dim68-69-gamma72-cross-sections/scripts/presentation_invariance.py`         | the bound does not depend on which basis presents the Gram      | 80 s   |
| `…/dim68-69-gamma72-cross-sections/scripts/k4_gram_min_det.py`                 | no four-tuple of Γ₇₂ minimal lines beats D₄                        | 2 min  |
| `…/dim70-71-gamma72-cross-sections/derive.py`                                 | τ(70) ≥ 1 249 778 250, τ(71) ≥ 2 603 658 750                   | 11 s   |
| `…/dim73-95-gamma72-caps/scripts/calibrate.py`                                | **the calibration**: reproduces Cohn's table in dims 26, 28–31 | 2 min  |
| `…/dim73-95-gamma72-caps/scripts/verify_classes.py data/disjoint_classes.npz` | the 120 Γ₇₂ classes, exactly                                   | 46 s   |
| `…/superseded/dim44-45-p48-cross-sections/brackets.py`                        | the withdrawn claims, and Sun–Wang inside the brackets         | 40 min |
| `…/closed/dim22-23-maximal-cross-sections/verify.py`                          | Λ₂₁, Λ₂₂, Λ₂₃ are maximal                                      | 6 s    |
| `…/closed/dim32-44-ers-audit/PIPELINE.py`                                     | dimensions 32–44 are exactly at the ERS value                  | 1 s    |
| `…/closed/dim17-23-cohn-li-mechanism/scripts/verify19.py`                     | the dimension-19 Cayley graph structure                        | 3 s    |
| `…/closed/dim17-24-layer-identity/layers.py`                                  | the layer identity, as a structure inside the Leech lattice     | 30 s   |
| `…/closed/dim17-24-layer-identity/flats.py`                                   | the flat layer exactly; dims 20 and 21 rebuilt from the Golay code | 80 s |
| `…/closed/dim17-24-layer-identity/equator.py`                                 | the deep-hole radius rule, and A(16,8,6) = 16                   | 20 s   |
| `…/closed/dim96-ers-takeover/ers96.py`                                        | dimension 96 belongs to Edel–Rains–Sloane, at 12 886 999 232    | 20 s   |
| `…/closed/dim96-ers-takeover/ers_sweep.py`                                    | **the level-0 floor in every claimed dimension** — read this    | 1 s    |
| `…/closed/dim96-ers-takeover/above96.py`                                     | why the table stops at 96, and not somewhere arbitrary          | 4 s    |


**Requirements.** Python 3.8+, numpy, scipy, and sympy for `paper/formulas.py`, which
re-derives the manuscript's algebra symbolically. The dimension-27 package needs neither
numpy nor scipy for its main checker. OR-Tools is needed only by `alpha512.py` in
`closed/dim17-23-cohn-li-mechanism/`, which re-proves α(G₀) = 160 optimal; its output is
shipped as a log so the claim is readable without it.

Two scripts want several gigabytes of memory — the all-pairs sweep and the class
regeneration — and are better run on their own than inside `--full`.

## What is not in this repository

Three things the working directory had that this repository deliberately does not:

- **The large intermediate files.** The P₄₈ and Γ₇₂ minimal-line pools, the 32 000 explicit
classes — several gigabytes. What is shipped instead is the *small seeds* they are built
from (Gram matrix, published automorphism generators, one base class, a size table) and
deterministic regeneration scripts that rebuild and re-verify them. `regenerate.py` for
dimensions 49–63 reproduces the shipped size vector exactly, and asserts pairwise
disjointness while doing it.
- **H. Cohn's coordinate data sets.** `dimensions1-24.txt` and `dimensions25-31.txt`, 120 MB,
from [https://hdl.handle.net/1721.1/153312](https://hdl.handle.net/1721.1/153312). Only two scripts use them and both take the
path as an argument: `superseded/dim27-triple-partition/scripts/rebuild_from_published.py`, which
regenerates the dimension-27 coordinate file bit-for-bit from the published block, and
`closed/dim09-19-record-maximality/sweep.py`. Everything else verifies from what is here.
- **The dead ends that produced nothing citable** — for instance the 249 MB of dimension-37
cover searches. Their *conclusions* are in `verifications/closed/`, with the measurements
that support them; the intermediate data is not.



## If you are picking this up

The repository is meant to be continued, not just read. A suggested route:

1. `python audit.py`, then `python run_all.py` — establish for yourself that it does what it
  says before trusting any of it.
2. `[verifications/closed/](verifications/closed/)` — **read this before starting anything.**
  Eight packages, each proving some mechanism to be at its exact ceiling. Several looked
   promising right up to the point where they were shown to be dead, and the READMEs say
   what the measurement was.
3. `[verifications/closed/dim17-layered-family/](verifications/closed/dim17-layered-family/)`
  — **dimension 17, reduced to a single integer and then settled**: τ(17) = 5346 + 2α with
   α = 192 exactly, so τ(17) = 5730 and Cohn–Li's configuration is optimal in this family.
   Worth reading for the method rather than the number: no amount of CP-SAT could separate
   192 from 193, and the graph turned out to be a Cayley graph on `F_2^10`, where Delsarte's
   LP settles it in fifteen seconds with a dual of six rational numbers.
4. `[verifications/closed/class-problem-upper-bound/](verifications/closed/class-problem-upper-bound/)`
  — the highest-leverage open question here. A class of **249** lines (the record is 248,
   the best upper bound 425.45) would immediately give +2, +8, +16, +32, +52, +96 and +168
   in dimensions 25 through 31 — one extra *line* is two extra vectors, so it is worth
   twice the per-vector figures quoted in KNOWLEDGE.md §1. Every two-point method, the three-point SDP and the
   subconstituent split have been tried and all give exactly 4680/11; closing the gap needs a
   four-point relaxation or a lattice-structural argument.
5. `[common/PROOF-kpoint.md](common/PROOF-kpoint.md)` §8 — what would take the cross-section
  method past dimension 44: the Siegel theta series of degree ≥ 3, which would replace these
   bounds by exact values in dimensions 44–47 and 68–71.

Two habits this project learned the hard way, both worth keeping:

- **Re-fetch the published table before claiming anything.** It moved twice during this work
(dimensions 12 and 19), and the range above dimension 48 — where no table has entries at
all — is the easiest place to be wrong about what "previously known" means.
- **Test new LP code on a pair of equivalent Gram matrices.** (3,−3,0) and 3·A₃ must give the
same answer. A missing-moments bug hid behind that for a day.



## The write-up

*New lower bounds for kissing numbers in dimensions 18 through 96* is a 47-page account of
the two mechanisms, with the results as numbered theorems: the cross-section method in §4
(dimensions 68, 69, 70 and 71, and dimensions 46 and 47 recovered) and the cap
construction in §5 (dimensions 25, 26, 27, 28, 29, 30, 31, 38, 49–61, 73–95), with the
code-theoretic dimensions 39, 62, 63 and 96 in §6; dimension 18, the bent-coset hexagon,
postdates the write-up and is documented in its package only.

The October 2026 family/head constructions postdate this manuscript snapshot.
They are documented in their verification package; the manuscript files have
not been updated. Use `--skip-paper` to keep them unchanged during maintenance.

**It is in this repository, at [`paper/`](paper/)** — source, bibliography, `.bbl`, PDF and
three checkers that read this repository and exit non-zero on a mismatch:

```bash
cd paper
python factcheck.py      # 686 checks: every figure in the write-up, recomputed from here
python unsupported.py    # no figure in the write-up lacks a home here
python formulas.py       # 257 checks: every formula re-derived with sympy
```

`run_all.py` runs all three unless `--skip-paper` is supplied. Until 2026-08-30 the write-up was kept in a sibling directory of
the working tree, on the grounds that it cites the repository — which is true, is the ordinary
arrangement for a paper with a verification package, and meant that those nine hundred checks
**skipped in every clone**, so no reader could run the part that ties the manuscript to the
computations. They still skip rather than fail in a checkout with no `paper/`.

The write-up is deliberately narrower than this repository — the obstruction results get one
page rather than eight packages — and it is not a substitute for the verification code. Every
number in it is cross-checked against [`RESULTS.md`](RESULTS.md).

It compiles with `pdflatex` twice; `kissing46.bbl` is tracked so that a build without `bibtex`
still resolves all forty-five citations. Its session-by-session work-log is **not** here: it
is a first-person record of what went wrong, and it stays with the author's private notes.

## Licence, citation, and how this was made

MIT (see `[LICENSE](LICENSE)`). `[CITATION.cff](CITATION.cff)` carries the metadata GitHub
uses for its "Cite this repository" button.

The constructions, the verification code and the write-ups in this repository were developed
by the author working with **Claude Code** (Anthropic). The repository is deliberately
structured so that someone else can pick it up the same way: every package states its idea in
prose before its code, every negative result says what was tried and why it failed, and
`KNOWLEDGE.md` is the full working record — 158 sections, most of them about things that
did not work.

Alexey Kravatskiy, MIRIAI (Moscow Independent Research Institute of Artificial Intelligence),
Moscow, Russia. [kravatskii.a@miriai.org](mailto:kravatskii.a@miriai.org)
