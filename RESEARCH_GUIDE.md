# Construction and verification guide

Detailed methods, verification scope and reproduction instructions. [Current results](README.md#results-and-construction-history) · [Weekly history](HISTORY.md#recent-updates).

## Reading the results table

The **recorded comparison** is the historical baseline used when each result was
obtained. **Comparison method / source** names its construction and links to the
paper or certificate supporting that baseline. Entries marked *derived* are
comparison floors computed by direct sum or embedding. In a direct sum, the
parent configuration and the indicated tail occupy orthogonal subspaces, so their
kissing counts add. The paper links document the parent constructions; the tail
links point to the [pinned configuration counts and references](common/published.py)
used in the calculation. The 80D direct sum is also a historical Nebe–Sloane table
entry, as recorded there.

**Gain over recorded comparison** measures the increase over the historical
baseline in the same row: `(this work / recorded comparison - 1) * 100%`,
or the ratio `this work / recorded comparison` when it is at least 2.
That baseline can precede an intervening record by another author. For example,
the 25D row compares 197,585 with Ma et al.'s historical 197,056, a gain of
529 points (+0.268%). The advance over Takhanov-Yun's subsequent 197,583-point
construction is **two points**, and the gain over this repository's preceding
197,584-point certificate is **one point**. In 31D the current 251,395-point construction gains
**13,045 points** over the historical 238,350-point comparison. Its advance
over PackingStar's November 11 release of 251,254 is **141 points**, and its
increase over the preceding published repository value 251,200 is **195 points**. The
[25D](verifications/improved/dim25-critical-cluster/),
[historical 31D](verifications/improved/dim31-dual-centres/) and
[retained 26–31D residual completions](verifications/improved/dim26-31-residual-completions/) packages
name these comparisons and the associated constructions.

**Later gain** compares a superseding SOTA with our retained result, using the
same percentage/ratio convention. Small percentages retain enough digits to
show improvements of only a few points. **Our method** links to our verification
package. **Superseding method** names the method behind the displayed
larger SOTA and links to its repository or article. Full construction notes are
expandable below.

The **SOTA value** column gives the larger current result
where ours has been surpassed, and **SOTA** where our value remains best in the
reviewed sources. A larger SOTA value marks a superseded result;
our original values and construction descriptions are retained. All registered
construction results are listed below in increasing dimension order, including small
improvements. See
[record status](HISTORY.md#record-status-and-retained-results) for the later values and sources.

The status notes use [Cohn's leaderboard](https://cohn.mit.edu/kissing-numbers/),
retrieved live on **17 November 2026**, with all 48 rows and explicit source credits captured in
[the dated comparison](common/cohn_20261117.json). It credits this repository in 18D and 25–46D;
all 23 listed values matched the repository at that checkpoint. The [25 October snapshot](common/cohn_20261025.json)
and [26 October source review](common/source_review_20261026.json) remain historical records.
All 70 constructions exceed their recorded historical comparisons. [Takhanov’s 17 November release](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1) supersedes the retained 26D and 27D values. The November 14 review records complete certificates
exceeding PackingStar's November 11 release in every dimension from 26D through 31D. See the [current comparison](common/source_review_20261114_26_31.json)
and the retained [October 29 review](common/source_review_20261029_27_31.json).
Higher-dimensional status uses the separate sources below.
The phase continuation separately records its [comparison scope](common/source_review_20261028_coxeter_todd.json).

## Construction notes by dimension

The 28 October nonlinear-code review corrects the earlier 92–96D comparisons.
Those pre-correction values were below the published DGLS code and its
shortenings, so their earlier record labels do not establish literature priority.
The new rows use the corrected baselines and the explicit section/layer proofs.

The historical SOTA period covers this project's run of published records in each
dimension, including our successive improvements before the displayed value. A
“since” date marks a current record in the reviewed sources; superseded records also
have an end date, taken from the first dated larger construction identified below.
Dates use Moscow time. The 24 October entries mark this integration; the earlier
dates describe repository publication history rather than private discovery dates.
The start dates follow the
[first release, 30 August](https://github.com/alexlegeartis/KissingNumbers/commit/75a4047),
[26D update, 8 September](https://github.com/alexlegeartis/KissingNumbers/commit/3e2f0fc),
[28-31D update, 16 September](https://github.com/alexlegeartis/KissingNumbers/commit/457c4d3),
[18D release, 17 September](https://github.com/alexlegeartis/KissingNumbers/commit/32c4e2a),
and [41-46D integration, 23 October](https://github.com/alexlegeartis/KissingNumbers/commit/5bed851).
The end dates are [Yang et al., v1, 28 September](https://arxiv.org/abs/2609.35051v1)
for 25D, 27D, 38D, 39D and 49–55D;
[Takhanov–Yun, 5 October](https://github.com/k-nic/kissing_d25-31/commit/b3768310ba566da5c8947546266c103fbcdce181)
for 26D, 29D and 31D; and
[K. Li, 10 October](https://github.com/Felpix-Studios/kissing-number/commit/365a76854e70eecc07dfefce7fb17df1617f2b9a)
for 28D and 30D. The later-result table below gives the latest retrieved values,
which can be larger than the first result that superseded ours.

The table above includes all 70 registered constructions. The October update adds
36D, 37D, 64D, 65D, 66D and 67D, raises 28D, 43D and 45D, installs the stronger P48 class family,
and registers code and coset improvements through 96D. [RESULTS.md](RESULTS.md) gives the generated registry
and the detailed construction descriptions.
Its comparisons use the repository's pinned data and include results that have since
been superseded. The manuscript in [paper/](paper/) retains its documented earlier snapshot,
including its historical appendix. See [The write-up](#the-write-up).

The [Coxeter–Todd package](verifications/improved/dim36-37-coxeter-todd/)
proves 501 813 points in 36D by exact short-shell enumeration, individual
base-triangle switches and joint equatorial refill, then reaches 505 629 in 37D
while retaining all latitude orbits. Two independent coordinate proofs check every
affected pair category. The
[P48n family/head package](verifications/improved/dim41-46-p48-heads/) supplies
the six October 41–46D constructions. In October the 43D certificate checked 3 450 individual small-body points, and the 45D construction used 2 187 complete families and 16 pure heads. The October angle checks use the documented population identities.

The November 41D and 42D constructions use P48q, and the 43D construction uses P48p. Their exact replays count every complete affine family from the Gram matrix and check all angular constraints. The 43D witness contains 3 165 individual small-body points. [Configurations and verification](verifications/improved/dim41-43-affine-sections/).

The 44D construction takes seventeen complete projected shells in the complement of a scaled D4 in P48m. The choice of normal vectors gives three empty projected shells of squared radius three. The resulting count is **4 940 961**. [Coordinates and moment proof](verifications/improved/dim44-d4-section/).

In 40D the distance-ten binary search covers all 7 497 quadratic-residue shortening types. Any three intact standard translates already limit an otherwise unrestricted binary sign layer to 589 824 words. In 45D the determinant-128 P48m section reaches 7 392 584 under the complete-family distance test. [Research calculations](research/40_48_affine_exchanges_20261122/).

<details>
<summary>Detailed construction notes</summary>

- **26–31D (13 November checkpoint):** **201 188, 209 884, 222 120, 231 814, 248 212 and 251 420**. The [26D/27D radical-orbit and head-frame exchanges](verifications/improved/dim26-27-radical-orbits/), [28D sextet unions](verifications/improved/dim28-29-november-sextets/), [29D native design](verifications/improved/dim29-30-native-partial-design/), [30D exact certificate unions](verifications/improved/dim30-certificate-unions/) and [31D annealed design](verifications/improved/dim31-annealed-design/) supply complete exact proofs. Earlier construction packages remain available.

- **29D** gains one point over Takhanov’s bound of 232,153 by replacing two E8 roots with four midpoints and repacking 19 auxiliary points into 18. [Complete exact proof](verifications/improved/dim29-takhanov-auxiliary-exchange/).

- **30D:** Two contact-sum rounds give **248,351**, 62 above the preceding 248,289-point witness. [Complete exact proof](verifications/improved/dim30-contact-sum-exchanges/).

- **18D:** the base of Cohn–Li's hexagon had a free combinatorial parameter nobody had varied: which coset of RM(1,4) the three six-set families lie in.  Three bent cosets with pairwise-bent sums that do not sum to zero give the twelve tier-B slots three different 512-word spaces instead of one, and tier B rises from 256 to 960

- **25D:** 1016 heads in the lens of a minimal vector, each removing only its owner, plus one non-lattice equator point

- **26D and retained 28D–30D constructions:** Independent signed-point choices and recentered candidate pools give 201 136, 221 694, 231 638 and 247 960. [Complete exact coordinate proofs](verifications/improved/dim26-30-asymmetric-exchanges/) handle the unpaired points explicitly.
- **27D (retained construction):** A 384-point quadratic normal module replaces PackingStar's 370-point sector, giving **209 600**, a gain of 14. [Full coordinates and exact all-pairs verification](verifications/improved/dim27-quadratic-normal-module/).
- **31D (retained construction):** A normal-coordinate signed simplex module and exact compatible-union selection give **251 110**, 108 above PackingStar's 28 October 251 002-point file. [Complete coordinates and all-pairs verification](verifications/improved/dim31-normal-module-exchange/) credit and pin the public baseline.
- **26D–30D and the preceding 31D construction:** Shared midpoint replacements with axiswise residual-head completion give 201 108, 209 532, 216 680, 226 520, 238 808 and 246 638. The final 30D/31D pure tails are represented exactly over Q(√3)/Q(√2), replacing—not supplementing—the ordinary D6/D7 tails. [Coordinate assembly and independent full replays](verifications/improved/dim26-31-residual-completions/) describe all components. Earlier midpoint, joint-Leech, frame-layer and dual-centre packages remain as historical constructions.

- **36D:** Three Coxeter–Todd lattices with 32 706 explicit three-point phase orbits and a 403 704-point base modified by 90 triangle switches give **501 822** points. The [latest orbit exchange](verifications/improved/dim36-37-orbit-exchange/) retains every base point and latitude body; fresh source reconstruction and exact integer products verify all changed interactions.

- **37D:** The 501 822-point 36D construction as an equator, plus two compatible latitude layers of 1 908 points each, gives **505 638**.

- **38D:** The Leech cap construction at codimension 14

- **39D:** [Edel–Rains–Sloane with the updated constant-weight codes](verifications/improved/dim39-ers-constant-weight/)

- **41D:** October: body families with one-dimensional heads, plus two pure-head poles

- **42D:** October: jointly replaced auxiliary layers with two-dimensional heads, plus a pure-head hexagon

- **43D:** October: 118 complete body families, 3 447 explicit tiny points selected by continuous head placement and shared-deletion trades, and five pure heads; 2 919 782 points with a supplied exact certificate

- **44D:** October: 1 034 complete body families with rational four-dimensional heads

- **45D:** October: coordinated movement of large, medium and small heads retains 2 187 complete families plus 16 pure heads; 7 655 056 points

- **46D:** October: independently normalized P48n families with six-dimensional heads, plus 27 pure-head points

- **49–60D:** [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) extend the cap construction by allowing individual signed owners without requiring each class to be antipodal. Exact replay checks classes, ownership, isometries and tails.

- **61–67D:** [Compatible support and sign families](verifications/improved/dim61-67-compatible-sign-families/) use shortened binary layers, exact head compression and balanced sections. The 67D construction also uses 49 smaller affine Hamming subcodes; full source replays and independent sign-pair checks verify the extensions.

- **68D:** Six affine Gamma72 slices with a D4 normal geometry give **471 886 208** points. An actual D7 normal extension strengthens the population count through exact moments and nested orientation averaging; [two independent proofs](verifications/improved/dim68-69-gamma72-extended-moments/) check the complete finite profile domain.

- **69D:** Ten Gamma72 slices are compressed from four normal coordinates to one rational-height coordinate. The [D7 moment extension](verifications/improved/dim68-69-gamma72-extended-moments/) raises the certified population mean to `6716925440/9`; its integer ceiling and two pure head poles give **746 325 051** points.

- **70D:** Three compatible affine Gamma72 slices give **1 360 986 112** points, with an [exact moment dual and projected geometry proof](verifications/improved/dim70-71-gamma72-cross-sections/).

- **71D:** Gamma72 one-point cross-section count; the recorded tables have no entry between 65 and 71

- **73–95D:** [Gamma72 cap construction with explicit disjoint classes and rotated axis points](verifications/improved/dim73-95-gamma72-caps/)

- **96D:** no table reaches it, and Edel–Rains–Sloane was never evaluated there

</details>

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
  improved/               all 70 claims, 178 packages; earlier constructions remain
                          available with their original proof data
  recovered/              values already published, re-derived here and claimed nowhere
  superseded/             claims this project made and then lost, and why
  closed/                 mechanisms pushed to their exact ceiling: what NOT to retry

KNOWLEDGE.md              the full working record, 160 sections, including everything
                          that failed -- with a preamble on how to read it, which
                          sections supersede which, and where each script now lives
```

[`research/`](research/) preserves the November campaigns, including the new
32+ search, restricted-family certificates and unsuccessful searches.
Its index lists the results and compact evidence replays.

Dimension-five upper-bound research is maintained separately in
[KissingUpperBounds](https://github.com/alexlegeartis/KissingUpperBounds).
The former `dim5upper/` snapshot was transferred to that repository's
`legacy/dim5upper-public-snapshot-20261023/` package, with its original files
and SHA-256 inventory, before being removed here. The unrestricted bound remains
`40 <= tau(5) <= 44`.


### `verifications/improved/` — 178 packages, 70 dimensions

Every current claim is in this tier. The historical 96D ERS evaluation remains
in `closed/dim96-ers-takeover/`; the new doubled-support construction is larger.

These directory names preserve the classification at packaging time. A construction
in `improved/` may now be superseded by a later publication; its package is retained.


| package                                                                                      | dims   | the idea                                                                                                                                                                                     |
| -------------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`dim62-partial-nonary-refill`](verifications/improved/dim62-partial-nonary-refill/) | 62 | Checked partial sign fibres and exact enumeration of the individual old blockers |
| [`dim71-gamma72-single-block-refill`](verifications/improved/dim71-gamma72-single-block-refill/) | 71 | An exact Leech-block refill of the four projected Γ72 slices |
| [`dim77-84-gamma72-half-root-shells`](verifications/improved/dim77-84-gamma72-half-root-shells/) | 77–84 | Tilted roots and orthogonal D4 extensions of Γ72; exact parity-resolved shortest-shell enumeration |
| [`dim36-37-orbit-exchange`](verifications/improved/dim36-37-orbit-exchange/) | 36, 37 | 501,822 and 505,638; nine-point gains from a compatible orbit exchange; fresh exact source and replacement replay |
| [`dim61-67-compatible-sign-families`](verifications/improved/dim61-67-compatible-sign-families/) | 61–67 | larger compatible support and sign families; fresh exact code and layer verification |
| [`dim68-69-gamma72-extended-moments`](verifications/improved/dim68-69-gamma72-extended-moments/) | 68, 69 | 471,886,208 and 746,325,051; actual D7 extension, exact population duals, independent direct-profile proof and two 69D poles |
| [`dim69-gamma72-head-compression`](verifications/improved/dim69-gamma72-head-compression/) | 69 | retained 743,522,304-point construction; ten slices with one reconstructed coordinate; two independent exact moment and geometry proofs |
| [`dim49-60-p48-signed-refill`](verifications/improved/dim49-60-p48-signed-refill/) | 49–60 | signed P48 class exchanges; twelve improvements with complete finite class, ownership and tail verification |
| [`dim29-30-native-partial-design`](verifications/improved/dim29-30-native-partial-design/) | 29,30 | retained 231,814 in 29D and 248,180 in 30D; native partial designs, exact exchanges and full geometry |
| [`dim30-leech-sloane-section`](verifications/improved/dim30-leech-sloane-section/) | 30 | retained 248,084; partial-fiber 1,356-point section, complete exact geometry and independent integer reconstruction |
| [`dim61-62-66-affine-reconstruction`](verifications/improved/dim61-62-66-affine-reconstruction/) | 61, 62, 66 | exact affine sign subcodes and independent replay |
| [`dim63-65-affine-reconstruction`](verifications/improved/dim63-65-affine-reconstruction/) | 63–65 | exact affine sign subcodes and independent replay |
| [`dim31-annealed-design`](verifications/improved/dim31-annealed-design/) | 31 | current 251,420; E8 orientation search, exact exchange and full geometry |
| [`dim31-joint-design`](verifications/improved/dim31-joint-design/) | 31 | historical 251,396; joint exact-cover design, exact exchange and full all-pairs verification |
| [`dim31-leech-sloane-module`](verifications/improved/dim31-leech-sloane-module/) | 31 | retained 251,395; 2,559-point module subset, complete exact geometry and independent integer reconstruction |
| [`dim30-bw14-module`](verifications/improved/dim30-bw14-module/) | 30 | retained 248,020 points; exact 932-for-960 Barnes–Wall module replacement, full geometry and independent integer checks |
| [`dim28-29-november-sextets`](verifications/improved/dim28-29-november-sextets/) | 28, 29 | current 28D 222,120 and retained 29D 231,776; exact finite geometry and independent direct-integer changes over the November 11 PackingStar files |
| [`dim31-latitude-contact-exchange`](verifications/improved/dim31-latitude-contact-exchange/) | 31 | retained 251,200-point four-for-eight exchange; two complete exact proofs, including independent integer-only enumeration |
| [`dim27-31-joint-refills`](verifications/improved/dim27-31-joint-refills/) | 27–31 | joint refills and reflected unions above the October 28 PackingStar release; complete exact proofs |
| [`dim26-30-asymmetric-exchanges`](verifications/improved/dim26-30-asymmetric-exchanges/) | 26, 28–30 | current 26D and retained 28–30D asymmetric integer configurations; complete independent all-pairs verification |
| [`dim27-quadratic-normal-module`](verifications/improved/dim27-quadratic-normal-module/) | 27 | retained 209 600-point quadratic code; full independent all-pairs verification |
| [`dim31-normal-module-exchange`](verifications/improved/dim31-normal-module-exchange/) | 31 | retained 251 110-point integer code; full independent all-pairs verification |
| [`dim26-31-residual-completions`](verifications/improved/dim26-31-residual-completions/) | 26–31 | retained residual-head completions and exact exceptional pure tails; independent full replays |
| [`dim25-critical-cluster`](verifications/improved/dim25-critical-cluster/) | 25 | critical-cluster replacement; 926 circles and 247 rational points; exact conic contacts and complete independent proof |
| [`dim25-point-split`](verifications/improved/dim25-point-split/) | 25 | one-for-two point splitting; 971 circles and 103 rational points; complete independent exact proof |
| [`dim25-collective-deformation`](verifications/improved/dim25-collective-deformation/) | 25 | retained 197,584-point predecessor: 971 circles and 102 moved points; exact biquadratic comparisons |
| [`dim31-dual-centres`](verifications/improved/dim31-dual-centres/) | 31 | shared midpoint replacements and dual E7 centres; exact coordinate replay |
| [`dim38-joint-latitudes`](verifications/improved/dim38-joint-latitudes/) | 38 | Orthogonal Hamming latitude replacement; two exact proofs with inherited Golay/Leech minimum-norm premise |
| [`dim39-deformed-orbits`](verifications/improved/dim39-deformed-orbits/) | 39 | reduced-symmetry dense sign orbits and retained sparse layers; two independent complete exact proofs |
| [`dim40-affine-shells`](verifications/improved/dim40-affine-shells/) | 40 | three affine lattice shells with independent exact minimum enumeration and complete bounded integer cross-products |
| [`dim40-affine-centre-refill`](verifications/improved/dim40-affine-centre-refill/) | 40 | alternate affine centre and thirteen individual short-shell points |
| [`dim18-bent-hexagon`](verifications/improved/dim18-bent-hexagon/) | 18 | bent-coset hexagon over odd Barnes–Wall; matches Cohn’s current table |
| `[dim25-lens-heads](verifications/improved/dim25-lens-heads/)`                               | 25     | cap heads of squared length 3 in the lens of a minimal vector remove exactly their owner and never share a removal, so τ = 196560 + H + 2 + \|E\| with H = 1016 and one non-lattice equator point |
| `[dim26-27-iota-triangles](verifications/improved/dim26-27-iota-triangles/)`                 | 26, 27 | retained construction: the axis forces the heads onto squared length 8/3 and onto fixed triangles of directions, so the layer is two or four head sets at cosine ≤ 1/4 cross-constrained at cosine 1/2; the coset triangle of three norm-6 vectors summing to zero is one such set, and the involution ι(u) = −v − u on each class puts S and ι(S) on two triangles at once |
| [`dim26-27-midpoint-replacements`](verifications/improved/dim26-27-midpoint-replacements/) | 26–27 | shared deletions for midpoint matchings on orthogonal tail directions |
| `[dim28-norm8-frame-layer](verifications/improved/dim28-norm8-frame-layer/)`                 | 28     | retained construction: the head of a deletion-free layer need not have norm 6: at height √2 a NORM-8 head clears the cap threshold where a norm-6 head does not, and the 24 vectors 8eᵢ are a Leech frame, each carrying all eight directions ±e_k of ℝ⁴ — 384 points for 8 axis points |
| [`dim28-joint-leech-layers`](verifications/improved/dim28-joint-leech-layers/) | 28 | joint selection of even heads and auxiliary layers over Li's odd Leech blocks |
| `[dim29-30-frame-layer](verifications/improved/dim29-30-frame-layer/)`                           | 29, 30 | the dimension-28 norm-8 frame layer is not about dimension 28: at height `√2` with `x = v/2`, `|v|² = 8`, the layer deletes nothing, one head carries the whole cross-polytope of `ℝᵏ` and all 48 vectors of a Leech frame carry every direction, so it is `96k` points in every dimension. The price is that all `⌊τ(k)/3⌋` owner classes must be type B for one frame — and in the coordinates where the frame is `8eᵢ` the type-B lines are the octad vectors of the Golay code, which the monomial group `2¹²:M₂₄` permutes, so the filter disappears and what is left is a packing solved by coordinate descent on the Golay sign words |
| `[dim31-frame-layer](verifications/improved/dim31-frame-layer/)`                             | 31     | the norm-8 frame layer at height √2 carried to k = 7: all 48 vectors of a Leech frame on the whole cross-polytope of ℝ⁷, 672 points, deleting no equator point.  It needs 42 pairwise disjoint type-B classes holding at least 10 252 of the 10 416 lines 42 full classes would give; the published attempt reached 10 183 and this one carries 10 328, found by screening candidates on the free slots their octads still have rather than by scanning uniformly — a sampled clique search needs pair density p ≥ 0.5097 to reach 42 at all, and flattening the octad coverage takes p from 0.4984 to 0.5467.  The axis is E7 rotated by 45°, 45°, 90° in the coordinate planes (0,3), (2,4), (5,6), keeping 118 of 126 against the 110 the published axis kept |
| `[dim38-leech-large-codimension](verifications/improved/dim38-leech-large-codimension/)`     | 38     | the same construction at codimension 14, where the binding constraint flips and the whole Leech shell partitions into 644 classes                                                            |
| `[dim39-ers-constant-weight](verifications/improved/dim39-ers-constant-weight/)`             | 39     | Edel–Rains–Sloane with the 2026 constant-weight codes, at n₀ = n rather than n₀ = 32                                                                                                         |
| [`dim32-35-signed-layers`](verifications/improved/dim32-35-signed-layers/) | 32–35 | signed octad replacements, dense-code coset heads and parity-coordinated root replacements; two independent exact checkers |
| [`dim36-37-coxeter-todd`](verifications/improved/dim36-37-coxeter-todd/) | 36, 37 | Coxeter–Todd phase orbits; exact short-shell and angle replay; two latitude layers |
| [`dim49-60-p48-pointwise-refill`](verifications/improved/dim49-60-p48-pointwise-refill/) | 49–60 | twelve pointwise-refill improvements, up to 58,627,711 in 60D; two exact head proofs and complete tail/pole checks |
| [`dim60-antipodal840-lift`](verifications/improved/dim60-antipodal840-lift/) | 60 | 58,543,525; 840 antipodal directions and 841 independent poles over 420 P48 classes |
| [`dim41-46-p48-heads`](verifications/improved/dim41-46-p48-heads/) | 41–46 | rational heads, explicit partial 43D families and coordinated 45D movement; inherited complete-family populations |
| `[dim49-63-p48-caps](verifications/improved/dim49-63-p48-caps/)`                             | 49–59 (60 retained) | the cap construction over P₄₈, with the 7 783-line class and 577 explicit disjoint images                                                                                             |
| `[dim62-63-ers-chain](verifications/improved/dim62-63-ers-chain/)`                         | 62, 63 | retained Edel–Rains–Sloane chain (n, 15, 2); it beats the caps restricted to the shipped bases and is superseded by the October central sections |
| `[dim68-69-gamma72-cross-sections](verifications/improved/dim68-69-gamma72-cross-sections/)` | 68, 69 | realizability of a k-point Gram settled by **exhibiting** the k-tuple in Γ₇₂ rather than through the LP's own witness lemma, which unlocks k = 3 and k = 4 |
| `[dim70-71-gamma72-cross-sections](verifications/improved/dim70-71-gamma72-cross-sections/)` | 70, 71 | central sections and the exact three-slice affine continuation on Nebe's Γ₇₂                                                                                                                                                        |
| `[dim73-95-gamma72-caps](verifications/improved/dim73-95-gamma72-caps/)`                     | 73–95  | the cap construction over Γ₇₂, whose class threshold γ = 1/4 is exactly what permits zero-sum triples                                                                                        |
| [`dim88-96-z8-golay`](verifications/improved/dim88-96-z8-golay/) | 88–96 | nonlinear published Golay parent, exact sections and sparse layers; two independent proofs |
| [`dim61-68-95-96-code-layers`](verifications/improved/dim61-68-95-96-code-layers/) | 61–63, 66–68, 95–96 | binary compression, central sections and enlarged support codes |
| [`dim64-65-signed-parity`](verifications/improved/dim64-65-signed-parity/) | 64-65 | retained 1,835-block Hamming augmentation, exact parity repairs and BCH owner caps |
| [`dim73-88-gamma72-cosets`](verifications/improved/dim73-88-gamma72-cosets/) | 73–88 | certified nearest subsets, repaired and scaled group images |


Each package has a README with the idea, a description of every file, and **the history of
the lower bound in its dimensions** — who held it, when, and by what construction.

### `verifications/recovered/` — nothing claimed


| package                                                                               | dims   | what it is                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[dim46-47-p48-cross-sections](verifications/recovered/dim46-47-p48-cross-sections/)` | 46, 47 | Venkov's theorem applied to **products** of inner products with several fixed minimal vectors. Both values were already published; re-deriving them from a method that knows nothing of either is the sharpest check the cross-section machinery has, which is why the package is kept |



## The construction ideas

**Coxeter–Todd phase codes (36–37).** Three rank-12 blocks carry a fixed
403 704-point base and 32 661 selected three-point phase orbits. Exact shell
enumeration and integer angles prove squared distance at least 12. Two latitude
layers add 3 816 points in 37D. The earlier full-cell and phase constructions are retained
in the package history; the new latitude layer uses partial three-point orbits.

**Independently normalized body families and heads (41–46).** The October constructions
use a rank-40 complement of sqrt(3) E8 in P48n, with heads in one to six dimensions.
Each family is normalized separately, so families can have different original lengths.
Joint class replacements and head motions admit additional complete or partial families;
zero-body points in the existing head span give the October 41D, 42D, 45D and 46D
extensions. The supplied 43D construction replaces its tiny layer with 3 447 explicit
body/head pairs, while coordinated 45D motion admits further complete families.
Exact checks, population assumptions and fresh
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

**The cap construction (25–27–31, 38, 49–63, 73–95).** In ℝⁿ ⊕ ℝᵏ, put a lattice shell on the
equator, delete a few lines, and put "caps" above the deleted directions. The optimum uses
**zero-sum triples** of cap directions rather than antipodal pairs, which needs class
threshold γ ≤ 1/4 and cap level t = 2/3; three is the maximum possible, because *n* points
pairwise at cosine ≤ −1/2 force *n* ≤ 3. Run over the Leech with its 248-line class this
reproduces Cohn's table in dimensions 26, 28, 29, 30, 31 **exactly** and beats it in 25 and
27 — a strong check on the model, and it is run as such by
`[verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py](verifications/improved/dim73-95-gamma72-caps/scripts/calibrate.py)`.

**Axis-layer comparison.** In the retained 73–95D cap construction, the shell rule
became the search space, and cost 221 025 spheres relative to the shipped
rotated and record axis layers (229 743 rather than 8 718 points across the
23 layers). This is a comparison within that construction, not an additional
gain to add to the current record table; `audit.py` recomputes it from the files.

**Explicit classes.** The gain is proportional to the size of a *class* — a set of lines with
pairwise |⟨u,u′⟩| ≤ μ/4 for Γ₇₂ and ≤ μ/3 for P₄₈. Caro–Wei bounds it with no coordinates at all, but real search does
about ten times better in every lattice tried: 221 → **2118** for Γ₇₂, 712 → **7783** for
P₄₈, and 32 000 pairwise-disjoint classes totalling 47 150 230 lines.

**The Edel–Rains–Sloane chain, evaluated where its authors did not (39, 62, 63).** Their count
is a sum over a chain of support sizes, `N(n) = Σ_ν A(n, n_ν, n_ν)·A(n_ν, ⌈n_ν/4⌉)`, and they
evaluated it only at 32, 36, 40, 44, 64, 80 and 128. At `n = 39` with the 2026 constant-weight
tables it gives 756 116; at `n = 62` and `n = 63` the chain `(n, 15, 2)` gives 71 310 732 and
138 419 844. Its levels above zero are **built** here rather than cited — a cyclic `[15,6,8]₄`
through the grid map, and a `[15,10,4]` — and the whole chain is checked in coordinates at
`n = 8` and `n = 16`, where it reproduces τ(8) = 240 (the E₈ root system) and τ(16) = 4320 (the
Barnes–Wall value, and the record) exactly.


## Experimental discovery workbench

The additive [discovery package](discovery/README.md) runs composable construction
programs for group representations, finite-ring codes and certified coset unions.
It includes exact known-result benchmarks, independent proof replay and resource
limits for higher-dimensional experiments. This initial campaign has found no new
record and does not change the registered results above. Earlier exploratory
scripts and receipts are preserved in `research/group_search_20261104/`.

```sh
python -m discovery operators
python -m unittest discover -s discovery/tests -t .
```


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
lattice and with every component at a stated ceiling but one; the Cohn–Li mechanism has
method-specific ceilings in dimensions 19–21; Λ₂₁, Λ₂₂ and Λ₂₃ have exact min-max cosines
√(8/29), √(3/11), √(4/15); the antipode construction is a max-weight clique problem whose
arms lie in a 30° cap, which caps it structurally; and the maximum-class problem survives
every two-point method, the three-point SDP and the subconstituent split at 4680/11 = 425.45.

`common/checkpoint.py` checks complete integrity manifests and runs the new
source verifiers in temporary copies. `source/dim61_96/verify_support_bits.py`
in the code-layer package independently checks every support pair using packed integers.
The [code-layer package](verifications/improved/dim61-68-95-96-code-layers/)
and [nearest-coset package](verifications/improved/dim73-88-gamma72-cosets/) state their
exact replay scope and inherited inputs.


## Reproducing everything

For this repository update, use `python run_all.py --skip-paper` or
`python update.py --skip-paper --repository-only`. The manuscript checks
remain tied to its earlier registry scope until the deferred paper update.
`python run_all.py --sota` verifies the current repository constructions.

```bash
python run_all.py          # 256 scripts, including exact certificate replays
python run_all.py --changed   # skip jobs whose inputs have not moved since they passed
python run_all.py --only dim31 # just the jobs whose label or command matches
python run_all.py --full   # 279 scripts: adds the all-pairs sweep, the
                           # negative controls, the class regeneration and the LP brackets
python run_all.py --list   # what would run, and roughly how long each takes
python run_all.py --only p48-heads  # October constructions: exact consistency replay
python run_all.py --sota --quick  # short replays for current SOTA-labelled constructions, with coverage gaps
python run_all.py --sota         # also run the slow coset and binary-code replays
python update.py --skip-paper --repository-only  # maintain constructions without changing the manuscript
```

[SOTA verification coverage](SOTA_VERIFY.md) maps every SOTA-labelled result to
its replay and inherited inputs. The quick profile reports deferred coset/code
work explicitly. Dimensions 88–94 use the shipped compact class certificate,
checked completely by `scripts/verify_symmetry.py` in under 10 CPU seconds.
Its `--complete` option also checks cap directions, partitions, poles and counts;
the SOTA profile runs that option. `--changed` reuses earlier passing checks
after hashing their inputs and environment.

Or individually, each from its own directory:


| script                                                                        | what it establishes                                            | time   |
| ----------------------------------------------------------------------------- | -------------------------------------------------------------- | ------ |
| `audit.py`                                                                    | all 70 claims, mutual and external consistency                 | 30 s   |
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
| `.../dim36-37-orbit-exchange/verify.py` | fresh exact baseline and orbit-exchange proof: 501 822 in 36D and 505 638 in 37D | ~1 min |
| `.../dim61-67-compatible-sign-families/verify.py` | fresh exact compatible-support, sign-family, code and layer checks for 61–67D | ~1 min |
| `.../dim68-69-gamma72-extended-moments/verify.py --enumerate` and `verify_independent.py` | two exact D7 population and geometry proofs: 471 886 208 in 68D and 746 325 051 in 69D | ~15 s |
| `.../dim36-37-coxeter-todd/verify.py` | exact 501 813-point 36D construction and 505 629-point 37D extension, two coordinate proofs | ~3 min |
| `.../dim44-d4-section/verify.py` | 44D lattice minimum, shell moments and angle checks | several minutes |
| `.../dim41-43-affine-sections/verify.py` | 41–43D shell counts and angles from integral coordinates | several minutes |
| `.../dim41-46-p48-heads/verify.py` | earlier 41–44D configurations and current 45–46D configurations with inherited populations | ~30 s |
| `…/dim46-47-p48-cross-sections/derive.py`                                     | recovers 12 309 600 and 23 766 960, claiming neither           | 6 s    |
| `…/dim49-63-p48-caps/verify.py` + `final.py`                                  |  the P₄₈ lattice, class, and fifteen numbers                    | 3 s    |
| `…/dim49-63-p48-caps/scripts/regenerate.py`                                   | rebuilds and re-verifies all 1400 disjoint classes             | 15 min |
| `…/dim49-63-p48-caps/scripts/ceiling.py`                                      | ceilings restricted to the shipped base classes, against level 0        | 1 s    |
| `…/dim62-63-ers-chain/scripts/verify.py`                                     | τ(62) ≥ 67 108 864, τ(63) ≥ 134 217 728, mechanism verified      | 1 s    |
| `…/dim62-63-ers-chain/scripts/verify_chain.py`                               | τ(62) ≥ 71 310 732, τ(63) ≥ 138 419 844, levels 1–2 built here | 20 s   |
| `…/dim68-69-gamma72-cross-sections/weyl_averaged/verify.py 3` and `verify.py 4` | exact averaged 689 452 922-point 69D and 471 884 288-point 68D existence bounds | ~3 s |
| `…/dim68-69-gamma72-cross-sections/affine68_69/verify.py 3` and `verify.py 4` | exact 689 442 146-point 69D and 471 366 656-point 68D existence bounds | ~12 s |
| `…/dim68-69-gamma72-cross-sections/scripts/leech_truth.py`                     | the k-point LP is **exact** on the Leech through k = 4          | 2 s    |
| `…/dim68-69-gamma72-cross-sections/scripts/recheck_cert.py`                    | both certificates re-derived by a second arithmetic path        | 40 s   |
| `…/dim68-69-gamma72-cross-sections/scripts/presentation_invariance.py`         | the bound does not depend on which basis presents the Gram      | 80 s   |
| `…/dim68-69-gamma72-cross-sections/scripts/k4_gram_min_det.py`                 | no four-tuple of Γ₇₂ minimal lines beats D₄                        | 2 min  |
| `…/dim70-71-gamma72-cross-sections/derive.py`                                 | τ(70) ≥ 1 360 986 112, τ(71) ≥ 2 603 658 750                   | 11 s   |
| `…/dim73-95-gamma72-caps/scripts/calibrate.py`                                | **the calibration**: reproduces Cohn's table in dims 26, 28–31 | 2 min  |
| `…/dim73-95-gamma72-caps/scripts/verify_classes.py data/disjoint_classes.npz` | the 120 Γ₇₂ classes, exactly                                   | 46 s   |
| `…/superseded/dim44-45-p48-cross-sections/brackets.py`                        | the withdrawn claims, and Sun–Wang inside the brackets         | 40 min |
| `…/closed/dim22-23-maximal-cross-sections/verify.py`                          | fixed cross-section checks for Λ₂₁, Λ₂₂, Λ₂₃                  | 6 s    |
| `…/closed/dim32-44-ers-audit/PIPELINE.py`                                     | dimensions 32–44 are exactly at the ERS value                  | 1 s    |
| `…/closed/dim17-23-cohn-li-mechanism/scripts/verify19.py`                     | the dimension-19 Cayley graph structure                        | 3 s    |
| `…/closed/dim17-24-layer-identity/layers.py`                                  | the layer identity, as a structure inside the Leech lattice     | 30 s   |
| `…/closed/dim17-24-layer-identity/flats.py`                                   | the flat layer exactly; dims 20 and 21 rebuilt from the Golay code | 80 s |
| `…/closed/dim17-24-layer-identity/equator.py`                                 | the deep-hole radius rule, and A(16,8,6) = 16                   | 20 s   |
| `…/closed/dim96-ers-takeover/ers96.py`                                        | dimension 96 belongs to Edel–Rains–Sloane, at 12 886 999 232    | 20 s   |
| `…/closed/dim96-ers-takeover/ers_sweep.py`                                    | **the level-0 floor in every claimed dimension** — read this    | 1 s    |
| `…/closed/dim96-ers-takeover/above96.py`                                     | why the table stops at 96, and not somewhere arbitrary          | 4 s    |


GPU class searches depend on the software stack used to assemble the pool. A
resumed search must match the saved pool fingerprint; the code rejects a mismatched
or missing fingerprint. See the [GPU reproduction instructions](verifications/improved/dim73-95-gamma72-caps/scripts/gpu/README.md)
for the supported workflow.

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
`closed/dim09-19-record-maximality/sweep.py`, which checks only the historical fixed-record
question. Everything else verifies from what is here.
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

*New lower bounds for kissing numbers in dimensions 18 through 96* is a 54-page
account with a state-of-the-art table in §5 and a historical appendix. The
main table reports 66 of the 70 registered claims. Dimensions 32–35 and
every improved value in the deferred checkpoint await the next manuscript
update. Every file under `paper/` remains unchanged; [the manuscript checkpoint](common/manuscript_scope.json)
records its published bounds and file hashes. The repository audit checks that
state, the four pending dimensions and every improved value.
`paper/registry_snapshot.json` separately preserves the 53-row arithmetic of
the historical September appendix.

**It is in this repository, at [`paper/`](paper/)** — source, bibliography, `.bbl`, PDF and
three checkers tied to its existing registry scope. Their synchronization is
part of the deferred paper update:

```bash
cd paper
python factcheck.py      # 779 checks: every figure in the write-up, recomputed from here
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


## Repository maintenance and archived campaigns

Below 24D, the archived campaigns have no new certified lower-bound improvement.
They include exact reproductions, nonclassical seeds and restricted-family
certificates. The engine's 1–24 target controls and held-out discovery benchmark
are measured separately in the [research index](research/README.md).

```bash
python audit.py --repository-only  # all 70 claims; check this repository without its surrounding workspace
```

`update.py` is the maintainer's entry point rather than the reviewer's: it runs every package
driver, then `audit.py`, then `audit.py --write-results`, then each package's table generator,
then the paper's tables and its checkers, then `run_all.py` -- in that order, so a number can
only travel outward from the driver that computed it. Run it after any change (`--check` to
report without rewriting anything). Nothing in it decides a value; steps 1-5 only move numbers
a driver already produced.

`audit.py` recomputes what is cheap to recompute — dimensions 47 and 71 from Venkov's
theorem, and registered counts from the shipped class, code, head and coset certificates — and checks them against
the pinned comparison data in `common/published.py`, for monotonicity and against records
in higher dimensions. These are checks against a recorded snapshot; they do not establish
that every packaged value remains the best known today. It takes a few seconds.
