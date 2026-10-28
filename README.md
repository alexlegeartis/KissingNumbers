# Kissing numbers: constructions and verification in dimensions 18 through 96

The kissing number τ(n) is the largest number of unit spheres that can touch a central unit
sphere in ℝⁿ; equivalently, the largest set of unit vectors with pairwise inner products at
most 1/2. This repository registers constructions in **70 dimensions**, each with a package
that states its verification scope, together with the negative results that guide further
searches. Some of our earlier records have since been superseded. Their constructions,
values and verification packages remain here as part of the research record.

> [!IMPORTANT]
> **28 October 2026 — further progress in 43D.**
> Continuous head placements after shared point deletions give **2 919 782**,
> improving the preceding repository certificate by **39 points**. Replacing
> 46 tiny-layer points with 85 gives 3 447 explicit tiny points. Two independent
> exact representations each replay **6 370 783** angle conditions.
> [Coordinates, method and reproduction](verifications/improved/dim41-46-p48-heads/continuation43/)
> retain the documented complete-family population inputs. The manuscript
> remains its earlier snapshot.
>
> **28 October 2026 — further progress in 49–72D: eighteen stronger bounds.**
> The **7,783-line P48 class** improves 49–60D again, by 210 to 74,090 points.
> A new eight-coordinate bridge gives **72 376 856** in 61D (+1 066 368).
> Parity-compatible sections give **100 437 548** in 62D (+66 080) and
> **182 537 988** in 63D (+130 560). Complete quadruple packings raise 66D,
> 67D and 68D by 10 240, 20 480 and 30 720 points. All gains compare with
> commit `9d38a38`; 64–65D and 69–72D retain their preceding bounds.
> [Code constructions and exact replay](verifications/improved/dim61-68-95-96-code-layers/source/continuation49_72/README.md)
> and [P48 certificates](verifications/improved/dim49-63-p48-caps/) state every
> inherited input and verification scope. The manuscript update remains deferred.
>
> **28 October 2026 — further verified progress in eight dimensions.**
> A refilled explicit 73D layer gives **6 218 206 932** (+2 856).
> Signed-tail existence proofs give **6 220 170 666**, **6 221 008 658**,
> **6 222 025 404**, **6 223 222 688** and **6 226 194 610** in 82–85D and
> 87D, gaining up to **1 217 848** points. [Exact proofs and replays](verifications/improved/dim73-88-gamma72-cosets/source/norm10/README.md#further-progress-in-73d-and-82--87d-28-october-2026).
> A further 14-support improvement gives **7 092 331 276** in 95D (+143 360)
> and **13 290 782 912** in 96D (+229 376), with all **303 675 690** support
> pairs checked in two representations. [Support continuation](verifications/improved/dim61-68-95-96-code-layers/search/continuation_20261028.json).
> Gains compare with the immediately preceding certificates; no new global
> priority survey is claimed. The manuscript retains its earlier snapshot.
>
> **28 October 2026 — verified improvements throughout 49–60D.**
> Fill/correct local replacements produce a **7,678-line P48 class**, replacing
> the 7,387-line base. Joint image selection and refilling give **52 431 358**
> in 49D (+582) through **57 894 253** in 60D (+203 582). Two complete exact
> geometric replays check all **29 472 003** base pairs, 577 isometries and
> **4 047 491** distinct mother lines. [Certificates, proof scope and reproduction](verifications/improved/dim49-63-p48-caps/).
> Gains compare with the preceding repository constructions. PackingStar's
> published fill/correct principle informed this CPU search; its larger reported
> 26–31D values are now recorded below. The manuscript retains its earlier snapshot.
>
> **28 October 2026 — verified phase refills in 36D and 37D.**
> The new certificates give **501 672** and **505 488**, improving the preceding
> repository bounds by **39** and **129** points. A shared-deletion phase search
> and finer latitude orbits pass two exact verification routes.
> [Proof, certificates and reproduction](verifications/improved/dim36-37-coxeter-todd/).
> PackingStar now reports larger bounds in 26–31D; its correction-and-refill
> strategy motivated this continuation. Those competing bounds are reported
> comparisons, not independently replayed here.
>
> **28 October 2026 — stronger exact existence bounds in 77–79D.**
> Mixed tail shells and unequal random priorities give **6 218 490 542**,
> **6 218 688 770** and **6 218 922 650**, gaining **25 526**, **2 890** and
> **10 614** points over the preceding certificates. Every tail product,
> rank and rational expectation is replayed exactly. These are existence
> proofs under the stated Gamma72 and design inputs; selected head lists
> are not enumerated. [Proof, PackingStar comparison and fast replay](verifications/improved/dim73-88-gamma72-cosets/source/norm10/README.md#mixed-tails-and-unequal-priorities-28-october-2026).
>
> **28 October 2026 — fill/correct gains in 95D and 96D.**
> The new 24,631-support certificate gives **7 092 187 916** in 95D and
> **13 290 553 536** in 96D, increasing the preceding repository bounds by
> **180 224** and **294 912** points. Coordinate and packed-integer replays
> independently check all **303 330 765** support pairs. The
> [search and certificate](verifications/improved/dim61-68-95-96-code-layers/)
> implement neutral exchanges and shared-blocker trades inspired by PackingStar's
> fill/correct principle; this is a finite search, with no optimality claim.
> [Source review](common/source_review_20261028.json) also records PackingStar's
> larger reported 26–31D bounds. Those six retained constructions are now marked
> as surpassed below. The manuscript remains its documented earlier snapshot.
>
> **27 October 2026 — stronger verified bounds in dimensions 26–31.**
> Exact residual-head completions give **201 108**, **209 532**, **216 680**,
> **226 520**, **238 808** and **246 638**, respectively. The increases over
> the preceding repository constructions are **190**, **4 896**, **11 448**,
> **9 952**, **14 528** and **4 896** points. Independent full exact replays
> passed on Python 3.9 and 3.11. [Coordinates, proofs and reproduction](verifications/improved/dim26-31-residual-completions/).
> The record table retains its named historical comparison baselines; these
> certificate checks are not a new worldwide literature survey or an optimality claim.
> The manuscript remains at its previous snapshot.
>
> **27 October 2026 — stronger bounds in 95D and 96D.**
> The new 24,613-support certificate gives **13 290 258 624** in 96D,
> improving the preceding repository bound by **851 968** points.
> Its central section gives **7 092 007 692** in 95D, an increase of **555 008**.
> Coordinate and packed integer replays each check all **302 887 578** support pairs.
> [Code-input sources](common/source_review_20261027.json) were reviewed on 27 October.
> These gains compare with our preceding certificate; the table below retains
> its historical comparison baselines. The paper update remains deferred.
>
> **26 October 2026 — four more SOTAs in dimensions 32–35.**
> Signed octad replacements give **352 800** in 32D and **370 368** in 33D.
> Dense-code coset heads give **390 982** in 34D; parity-coordinated twelve-head
> replacement gives **415 314** in 35D. These exceed the reviewed public Lin
> constructions by **480**, **128**, **2** and **6** points, respectively.
> Both independent exact verifiers passed. The repository now covers **70 dimensions**.
> The paper remains at its preceding 66-result checkpoint; its update is deferred.
>
> **26 October 2026 — new verified bounds in 25D, 31D, 38D, 39D and 40D.**
> Collective deformation gives **197 584** in 25D, one above Takhanov–Yun.
> Midpoints and dual E7 centres give **241 742** in 31D, exceeding Li by 663.
> Joint latitude families give **610 394** in 38D (+5 072 over Lin), and
> deformed sign orbits give **769 180** in 39D (+264 over Lin).
> Three affine lattice shells give **1 126 314** in 40D, exceeding
> Dorofeev–Sun–Wang's 1 092 000 by 34 314. Every new certificate is replayed
> independently; earlier constructions remain available.
>
> **25 October 2026 — new SOTAs in dimensions 26–30, 38–39 and further high-dimensional gains.**
> Shared midpoint replacements give **200 918** in 26D and **204 636** in 27D,
> exceeding Takhanov–Yun by **1 050** and **3 047** points.
> Midpoint replacements also give **216 568** in 29D and **224 280** in 30D,
> exceeding K. Li's records by **6 073** and **2 207** points.
> The retained 38D latitude extension gives **605 346** (+24 over Lin);
> the retained 39D root replacement gives **768 924** (+8 over Lin).
> Joint Leech layers give **205 232** in 28D, improving K. Li's 205 172 by 60.
> The new Coxeter–Todd phase code gives **501 633** in 36D and **505 359**
> in 37D. The supplied 43D certificate now reaches **2 919 743**, and
> coordinated head movement raises 45D to **7 655 056**. All have passed
> exact replay; complete-family P48n populations retain their documented assumptions.
> The update also includes the 7 387-line P48 family, code improvements in
> 61-68D and 95-96D, and Gamma72 nearest-coset constructions in 73-88D.
> The signed-support construction gives **332 170 112** in 64D and
> **332 238 338** in 65D; 72D is the remaining gap in 49-96D.
> In particular, 95D reaches **7 091 452 684** and 96D **13 289 406 656**.
> A compact symmetry certificate now also improves **87–94D**. Its complete
> class geometry and disjointness check takes **under 10 seconds on CPU**.

> The **25 October CPU continuation** strengthens the coset constructions in
> 76D, 78–86D and the cap constructions throughout 87–94D. Direct optimization
> of the 95D section also raises its bound to **7 091 452 684**, with a stronger 96D parent.

> The **73D non-lattice continuation** reaches **6 218 204 076**,
> exceeding the completed 73D lattice by **15 932**.
> `dim73-88-gamma72-cosets/verify73.py` independently checks its certificate.
> The new norm-ten method also improves **77–79D**, with exact whole-shell
> existence proofs checked by `verify_norm10.py`.

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

## Results and construction history

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
the 25D row compares 197,584 with Ma et al.'s historical 197,056, a gain of
528 points (+0.268%). The advance over Takhanov-Yun's subsequent 197,583-point
record is **one point**. In 31D the current 246,638-point construction gains
**8,288 points** over the table's historical 238,350-point comparison and
**5,559 points** over Li's 241,079-point record. Its separate gain over the
preceding repository construction of 241,742 is **4,896 points**. The
[25D](verifications/improved/dim25-collective-deformation/),
[historical 31D](verifications/improved/dim31-dual-centres/) and
[current 26–31D](verifications/improved/dim26-31-residual-completions/) packages
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
[record status](#record-status-and-retained-results) for the later values and sources.

The status notes use [Cohn's leaderboard](https://cohn.mit.edu/kissing-numbers/),
last retrieved **25 October 2026**, with all 48 rows captured in
[the dated comparison](common/cohn_20261025.json). Retrieval on 26 October
failed, so this update retains that snapshot. Fresh checks of the relevant
authors' repositories and the dimension-40 paper are recorded in
[the 26 October source review](common/source_review_20261026.json).
All 70 constructions exceed their recorded comparisons. Of these, 64 remain best in
the reviewed sources; PackingStar reports larger bounds in dimensions 26–31,
as recorded in [the 28 October source review](common/source_review_20261028.json). Higher-dimensional status uses the separate sources below.
The phase continuation separately records its [comparison scope](common/source_review_20261028_coxeter_todd.json).


| dim | recorded comparison | comparison method / source | **this work** | SOTA value | gain over recorded comparison | later gain | SOTA period | our method | superseding method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 18 | 7 654 | [Cohn–Li, BW16 + hexagons](https://arxiv.org/abs/2411.04916v3) | **8 358** | **SOTA** | +9.2% | — | since 2026-09-17 | [BW16](verifications/improved/dim18-bent-hexagon/) | — |
| 25 | 197 056 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **197 584** | **SOTA** | +0.268% | — |2026-08-30 to 2026-09-28; since 2026-10-26 | [Collective deformation](verifications/improved/dim25-collective-deformation/) | — |
| 26 | 199 868 | [Takhanov–Yun, deformed caps + free heads](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **201 108** | **201 134** | +0.620% | +0.01293% | 2026-09-08 to 2026-10-05; 2026-10-25 to 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 27 | 201 589 | [Takhanov–Yun, cap colouring](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **209 532** | **209 586** | +3.940% | +0.02577% | 2026-08-30 to 2026-09-28; 2026-10-25 to 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 28 | 205 172 | [K. Li, Leech lift and auxiliary layers](https://github.com/Felpix-Studios/kissing-number/blob/b823640f4e98134b9d55b27288132ee954a8affa/preprints/leech-lift-d28-d31/paper.pdf) | **216 680** | **221 664** | +5.609% | +2.30017% | 2026-09-16 to 2026-10-10; 2026-10-24 to 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 29 | 209 496 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **226 520** | **231 636** | +8.126% | +2.25852% | 2026-09-16 to 2026-10-05; 2026-10-25 to 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 30 | 220 440 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **238 808** | **247 952** | +8.332% | +3.82902% | 2026-09-16 to 2026-10-10; 2026-10-25 to 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 31 | 238 350 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **246 638** | **250 938** | +3.477% | +1.74345% |2026-09-16 to 2026-10-05; 2026-10-26 | [Residual-head completions](verifications/improved/dim26-31-residual-completions/) | [PackingStar, joint removals / insertions / auxiliary layers](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6/25D-31D/new_bounds_20261026) |
| 32 | 352 320 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-number-32) | **352 800** | **SOTA** | +0.13624% | — | since 2026-10-26 | [Signed octad replacements](verifications/improved/dim32-35-signed-layers/) | — |
| 33 | 370 240 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **370 368** | **SOTA** | +0.0345722% | — | since 2026-10-26 | [Signed octad replacements](verifications/improved/dim32-35-signed-layers/) | — |
| 34 | 390 980 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **390 982** | **SOTA** | +0.000511535% | — | since 2026-10-26 | [Dense-coset heads](verifications/improved/dim32-35-signed-layers/) | — |
| 35 | 415 308 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **415 314** | **SOTA** | +0.00144471% | — | since 2026-10-26 | [Parity-coordinated heads](verifications/improved/dim32-35-signed-layers/) | — |
| 36 | 484 760 | [Yang et al., signed replacements](https://arxiv.org/abs/2609.35051v1) | **501 672** | **SOTA** | +3.5% | — | since 2026-10-24 | [Coxeter–Todd phase refill](verifications/improved/dim36-37-coxeter-todd/) | — |
| 37 | 500 008 | [Olivier, using Lin's code](https://doi.org/10.5281/zenodo.23569153) | **505 488** | **SOTA** | +1.1% | — | since 2026-10-24 | [Phase code + latitudes](verifications/improved/dim36-37-coxeter-todd/) | — |
| 38 | 566 652 | [Edel–Rains–Sloane](https://arxiv.org/abs/math/0207291) + [Brouwer codes](https://aeb.win.tue.nl/codes/Andw.html) | **610 394** | **SOTA** | +7.72% | — | 2026-08-30 to 2026-09-28; since 2026-10-25 | [Joint latitude families](verifications/improved/dim38-joint-latitudes/) | — |
| 39 | 755 988 | [Edel–Rains–Sloane](https://arxiv.org/abs/math/0207291) + [Brouwer codes](https://aeb.win.tue.nl/codes/Andw.html) | **769 180** | **SOTA** | +1.75% | — | 2026-08-30 to 2026-09-28; since 2026-10-25 | [Deformed sign orbits](verifications/improved/dim39-deformed-orbits/) | — |
| 40 | 1 092 000 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/html/2607.20359v4) | **1 126 314** | **SOTA** | +3.14% | — | since 2026-10-26 | [Affine lattice shells](verifications/improved/dim40-affine-shells/) | — |
| 41 | 1 324 472 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/abs/2607.20359v4) | **1 460 036** | **SOTA** | +10% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 42 | 1 792 386 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/abs/2607.20359v4) | **1 991 445** | **SOTA** | +11% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 43 | 2 553 792 | [Yang et al., P48n section](https://arxiv.org/abs/2609.35051v1) | **2 919 782** | **SOTA** | +14% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 44 | 4 153 600 | Conway–Sloane P48 sections; [counts: Dorofeev–Sun–Wang](https://arxiv.org/abs/2607.20359v4) | **4 939 281** | **SOTA** | +19% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 45 | 7 380 090 | [Yang et al., QR neighbour projection](https://arxiv.org/abs/2609.35051v1) | **7 655 056** | **SOTA** | +3.7% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 46 | 12 309 600 | Conway–Sloane P48 sections; [counts: Dorofeev–Sun–Wang](https://arxiv.org/abs/2607.20359v4) | **13 387 797** | **SOTA** | +8.8% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 49 | 52 416 002 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [1D tail](common/published.py) *(derived)* | **52 431 568** | **SOTA** | +0.0297% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 50 | 52 416 006 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [2D tail](common/published.py) *(derived)* | **52 462 704** | **SOTA** | +0.08909% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 51 | 52 416 012 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [3D tail](common/published.py) *(derived)* | **52 509 404** | **SOTA** | +0.1782% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 52 | 52 416 024 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [4D tail](common/published.py) *(derived)* | **52 602 708** | **SOTA** | +0.3562% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 53 | 52 416 040 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [5D tail](common/published.py) *(derived)* | **52 726 898** | **SOTA** | +0.5931% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 54 | 52 416 072 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [6D tail](common/published.py) *(derived)* | **52 974 438** | **SOTA** | +1.065% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 55 | 52 416 126 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [7D tail](common/published.py) *(derived)* | **53 389 476** | **SOTA** | +1.857% | — | 2026-08-30 to 2026-09-28; since 2026-10-24 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 56 | 52 416 240 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [8D tail](common/published.py) *(derived)* | **54 253 958** | **SOTA** | +3.506% | — | since 2026-08-30 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 57 | 52 416 306 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [9D tail](common/published.py) *(derived)* | **54 747 156** | **SOTA** | +4.447% | — | since 2026-08-30 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 58 | 52 416 510 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [10D tail](common/published.py) *(derived)* | **56 237 690** | **SOTA** | +7.29% | — | since 2026-08-30 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 59 | 52 416 604 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [11D tail](common/published.py) *(derived)* | **56 907 532** | **SOTA** | +8.568% | — | since 2026-08-30 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 60 | 52 416 841 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [12D tail](common/published.py) *(derived)* | **57 968 343** | **SOTA** | +10.59% | — | since 2026-08-30 | [P48 caps, 7 783-line class](verifications/improved/dim49-63-p48-caps/) | — |
| 61 | 52 417 154 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [13D tail](common/published.py) *(derived)* | **72 376 856** | **SOTA** | +38.08% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 62 | 52 417 932 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [14D tail](common/published.py) *(derived)* | **100 437 548** | **SOTA** | +91.61% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 63 | 52 418 564 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [15D tail](common/published.py) *(derived)* | **182 537 988** | **SOTA** | +248.2% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 64 | 331 737 984 | [Edel-Rains-Sloane](https://arxiv.org/abs/math/0207291) | **332 170 112** | **SOTA** | +0.130% | — | since 2026-10-24 | [Signed support augmentation](verifications/improved/dim64-65-signed-parity/) | — |
| 65 | 331 737 986 | [ERS64](https://arxiv.org/abs/math/0207291) plus two poles *(derived)* | **332 238 338** | **SOTA** | +0.151% | — | since 2026-10-24 | [Signed supports and BCH owner caps](verifications/improved/dim64-65-signed-parity/) | — |
| 66 | 331 737 990 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291) ⊕ A2 *(derived)* | **338 302 340** | **SOTA** | +1.979% | — | since 2026-10-24 | [Code layers](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 67 | 331 737 996 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291) ⊕ A3 *(derived)* | **351 420 806** | **SOTA** | +5.933% | — | since 2026-10-24 | [Code layers](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 68 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **371 552 152** | **SOTA** | +12% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 69 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **627 822 180** | **SOTA** | +89% | — | since 2026-08-30 | [Gamma72 section](verifications/improved/dim68-69-gamma72-cross-sections/) | — |
| 70 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **1 249 778 250** | **SOTA** | 3.77x | — | since 2026-08-30 | [Gamma72 section](verifications/improved/dim70-71-gamma72-cross-sections/) | — |
| 71 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **2 603 658 750** | **SOTA** | 7.85x | — | since 2026-08-30 | [Gamma72 section](verifications/improved/dim70-71-gamma72-cross-sections/) | — |
| 73 | 6 218 175 602 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [1D tail](common/published.py) *(derived)* | **6 218 206 932** | **SOTA** | +0.000504% | — | since 2026-08-30 | [Gamma72 cosets plus caps](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 74 | 6 218 175 606 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [2D tail](common/published.py) *(derived)* | **6 218 211 606** | **SOTA** | +0.000579% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 75 | 6 218 175 612 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [3D tail](common/published.py) *(derived)* | **6 218 246 060** | **SOTA** | +0.00113% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 76 | 6 218 175 624 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [4D tail](common/published.py) *(derived)* | **6 218 315 288** | **SOTA** | +0.00225% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 77 | 6 218 175 640 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [5D tail](common/published.py) *(derived)* | **6 218 490 542** | **SOTA** | +0.005064% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 78 | 6 218 175 672 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [6D tail](common/published.py) *(derived)* | **6 218 688 770** | **SOTA** | +0.008252% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 79 | 6 218 175 726 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [7D tail](common/published.py) *(derived)* | **6 218 922 650** | **SOTA** | +0.012012% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 80 | 6 218 175 840 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ E8 (tabulated direct sum) | **6 219 474 300** | **SOTA** | +0.0209% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 81 | 6 218 175 906 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [9D tail](common/published.py) *(derived)* | **6 219 649 576** | **SOTA** | +0.0237% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 82 | 6 218 176 110 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [10D tail](common/published.py) *(derived)* | **6 220 170 666** | **SOTA** | +0.032076% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 83 | 6 218 176 204 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [11D tail](common/published.py) *(derived)* | **6 221 008 658** | **SOTA** | +0.045551% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 84 | 6 218 176 441 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [12D tail](common/published.py) *(derived)* | **6 222 025 404** | **SOTA** | +0.061899% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 85 | 6 218 176 754 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [13D tail](common/published.py) *(derived)* | **6 223 222 688** | **SOTA** | +0.081148% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 86 | 6 218 177 532 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [14D tail](common/published.py) *(derived)* | **6 225 165 512** | **SOTA** | +0.112% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 87 | 6 218 178 164 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [15D tail](common/published.py) *(derived)* | **6 226 194 610** | **SOTA** | +0.128920% | — | since 2026-08-30 | [Gamma72 signed-tail existence proof](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 88 | 6 218 179 920 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [16D tail](common/published.py) *(derived)* | **6 232 239 528** | **SOTA** | +0.226% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 89 | 6 218 181 330 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [17D tail](common/published.py) *(derived)* | **6 235 525 036** | **SOTA** | +0.279% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 90 | 6 218 183 254 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [18D tail](common/published.py) *(derived)* | **6 242 027 140** | **SOTA** | +0.383% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 91 | 6 218 187 548 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [19D tail](common/published.py) *(derived)* | **6 252 202 550** | **SOTA** | +0.547% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 92 | 6 218 195 048 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [20D tail](common/published.py) *(derived)* | **6 272 431 575** | **SOTA** | +0.872% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 93 | 6 218 205 368 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [21D tail](common/published.py) *(derived)* | **6 301 636 084** | **SOTA** | +1.34% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 94 | 6 218 225 496 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [22D tail](common/published.py) *(derived)* | **6 357 657 396** | **SOTA** | +2.24% | — | since 2026-08-30 | [Gamma72 caps](verifications/improved/dim73-95-gamma72-caps/) | — |
| 95 | 6 218 268 750 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [23D tail](common/published.py) *(derived)* | **7 092 331 276** | **SOTA** | +14.056365% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |
| 96 | 6 218 372 160 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [Leech](https://doi.org/10.4153/CJM-1967-017-0) *(derived)* | **13 290 782 912** | **SOTA** | +113.734118% | — | since 2026-08-30 | [Code layers / sections](verifications/improved/dim61-68-95-96-code-layers/) | — |


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
proves 501 672 points in 36D by exact short-shell enumeration and phase-orbit
angle checks, then reaches 505 488 in 37D with two latitude layers. The
[P48n family/head package](verifications/improved/dim41-46-p48-heads/) now supplies
all six 41–46D constructions. The 43D partial-family certificate checks 3 447 explicit
tiny points; the latest 45D certificate retains 2 187 complete families and 16 pure
heads. Their integer angles and exact head dimensions replay, with complete-family
populations inherited from the documented identities. The older missing 43D writeup
certificate is historical and is not needed for the current construction. No October
improvement is claimed in 47D or 48D.

<details>
<summary>Detailed construction notes</summary>

- **18D:** the base of Cohn–Li's hexagon had a free combinatorial parameter nobody had varied: which coset of RM(1,4) the three six-set families lie in.  Three bent cosets with pairwise-bent sums that do not sum to zero give the twelve tier-B slots three different 512-word spaces instead of one, and tier B rises from 256 to 960

- **25D:** 1016 heads in the lens of a minimal vector, each removing only its owner, plus one non-lattice equator point

- **26D–31D (current):** Shared midpoint replacements with axiswise residual-head completion give 201 108, 209 532, 216 680, 226 520, 238 808 and 246 638. The final 30D/31D pure tails are represented exactly over Q(√3)/Q(√2), replacing—not supplementing—the ordinary D6/D7 tails. [Coordinate assembly and independent full replays](verifications/improved/dim26-31-residual-completions/) describe all components. Earlier midpoint, joint-Leech, frame-layer and dual-centre packages remain as historical constructions.

- **36D:** Three Coxeter–Todd lattices with 32 656 explicit three-point phase orbits over a fixed 403 704-point base; exact shells and angles give 501 672 points

- **37D:** The 36D construction as an equator, plus two latitude layers of 1 908 points each

- **38D:** The Leech cap construction at codimension 14

- **39D:** [Edel–Rains–Sloane with the updated constant-weight codes](verifications/improved/dim39-ers-constant-weight/)

- **41D:** October: body families with one-dimensional heads, plus two pure-head poles

- **42D:** October: jointly replaced auxiliary layers with two-dimensional heads, plus a pure-head hexagon

- **43D:** October: 118 complete body families, 3 447 explicit tiny points selected by continuous head placement and shared-deletion trades, and five pure heads; 2 919 782 points with a supplied exact certificate

- **44D:** October: 1 034 complete body families with rational four-dimensional heads

- **45D:** October: coordinated movement of large, medium and small heads retains 2 187 complete families plus 16 pure heads; 7 655 056 points

- **46D:** October: independently normalized P48n families with six-dimensional heads, plus 27 pure-head points

- **49–61D:** [P48 cap construction with explicit disjoint classes](verifications/improved/dim49-63-p48-caps/)

- **62D:** Edel–Rains–Sloane chain (62,15,2); the recorded tables stop at 48 and resume at 64

- **63D:** Edel–Rains–Sloane chain (63,15,2); the recorded tables stop at 48 and resume at 64

- **68D:** Gamma72 four-point cross-section programme at its exact optimum, with a D4 Gram matrix

- **69D:** Gamma72 three-point cross-section programme at its exact optimum, with the Gram determinant forced by Hermite

- **70D:** Gamma72 two-point cross-section count; the recorded tables have no entry between 65 and 71

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
  improved/               all 70 claims, 30 packages; earlier constructions remain
                          available with their original proof data
  recovered/              values already published, re-derived here and claimed nowhere
  superseded/             claims this project made and then lost, and why
  closed/                 mechanisms pushed to their exact ceiling: what NOT to retry

KNOWLEDGE.md              the full working record, 158 sections, including everything
                          that failed -- with a preamble on how to read it, which
                          sections supersede which, and where each script now lives
```

Dimension-five upper-bound research is maintained separately in
[KissingUpperBounds](https://github.com/alexlegeartis/KissingUpperBounds).
The former `dim5upper/` snapshot was transferred to that repository's
`legacy/dim5upper-public-snapshot-20261023/` package, with its original files
and SHA-256 inventory, before being removed here. The unrestricted bound remains
`40 <= tau(5) <= 44`.


### `verifications/improved/` — 30 packages, 70 dimensions

Every current claim is in this tier. The historical 96D ERS evaluation remains
in `closed/dim96-ers-takeover/`; the new doubled-support construction is larger.

These directory names preserve the classification at packaging time. A construction
in `improved/` may now be superseded by a later publication; its package is retained.


| package                                                                                      | dims   | the idea                                                                                                                                                                                     |
| -------------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`dim26-31-residual-completions`](verifications/improved/dim26-31-residual-completions/) | 26–31 | current residual-head completions and exact exceptional pure tails; independent full replays |
| [`dim25-collective-deformation`](verifications/improved/dim25-collective-deformation/) | 25 | collective rational deformation of 971 Leech circles and 102 moved points; exact biquadratic comparisons |
| [`dim31-dual-centres`](verifications/improved/dim31-dual-centres/) | 31 | shared midpoint replacements and dual E7 centres; exact coordinate replay |
| [`dim38-joint-latitudes`](verifications/improved/dim38-joint-latitudes/) | 38 | 1352 joint latitude families, signed axis rotation and coordinate axes; exact integer replay |
| [`dim39-deformed-orbits`](verifications/improved/dim39-deformed-orbits/) | 39 | rationally deformed dense sign orbits and sparse root replacements; complete relative orbit enumeration |
| [`dim40-affine-shells`](verifications/improved/dim40-affine-shells/) | 40 | three affine lattice shells with independent exact minimum enumeration and complete bounded integer cross-products |
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
| [`dim41-46-p48-heads`](verifications/improved/dim41-46-p48-heads/) | 41–46 | rational heads, explicit partial 43D families and coordinated 45D movement; inherited complete-family populations |
| `[dim49-63-p48-caps](verifications/improved/dim49-63-p48-caps/)`                             | 49–60  | the cap construction over P₄₈, with the 7 387-line class and 577 explicit disjoint images                                                                                             |
| `[dim62-63-ers-chain](verifications/improved/dim62-63-ers-chain/)`                         | 62, 63 | retained Edel–Rains–Sloane chain (n, 15, 2); it beats the caps restricted to the shipped bases and is superseded by the October central sections |
| `[dim68-69-gamma72-cross-sections](verifications/improved/dim68-69-gamma72-cross-sections/)` | 68, 69 | realizability of a k-point Gram settled by **exhibiting** the k-tuple in Γ₇₂ rather than through the LP's own witness lemma, which unlocks k = 3 and k = 4 |
| `[dim70-71-gamma72-cross-sections](verifications/improved/dim70-71-gamma72-cross-sections/)` | 70, 71 | the same k-point method on Nebe's Γ₇₂                                                                                                                                                        |
| `[dim73-95-gamma72-caps](verifications/improved/dim73-95-gamma72-caps/)`                     | 73–95  | the cap construction over Γ₇₂, whose class threshold γ = 1/4 is exactly what permits zero-sum triples                                                                                        |
| [`dim61-68-95-96-code-layers`](verifications/improved/dim61-68-95-96-code-layers/) | 61–63, 66–68, 95–96 | binary compression, central sections and enlarged support codes |
| [`dim64-65-signed-parity`](verifications/improved/dim64-65-signed-parity/) | 64-65 | 211 added Hamming sign blocks, exact parity repairs and BCH owner caps |
| [`dim73-88-gamma72-cosets`](verifications/improved/dim73-88-gamma72-cosets/) | 73–88 | certified nearest subsets, repaired and scaled group images |


Each package has a README with the idea, a description of every file, and **the history of
the lower bound in its dimensions** — who held it, when, and by what construction.

### `verifications/recovered/` — nothing claimed


| package                                                                               | dims   | what it is                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[dim46-47-p48-cross-sections](verifications/recovered/dim46-47-p48-cross-sections/)` | 46, 47 | Venkov's theorem applied to **products** of inner products with several fixed minimal vectors. Both values were already published; re-deriving them from a method that knows nothing of either is the sharpest check the cross-section machinery has, which is why the package is kept |


## Record status and retained results

**28 October review.** PackingStar reports 201 134, 209 586, 221 664,
231 636, 247 952 and 250 938 in dimensions 26–31. All six exceed the
retained repository constructions, including the 27 October completions.
Their README describes joint point removals, midpoint insertions and auxiliary
layers; it supplies integer coordinate files. The external values above are
source-reported, not independently reverified by this update. The pinned
[review](common/source_review_20261028.json) separates their October claims
from the construction and algorithm described in the June paper.

The 25 October midpoint constructions regain 26D and 27D with **200 918** and
**204 636**, following Takhanov–Yun's 199 868 and 201 589. Their paper is linked
in the comparison-source field; our earlier coset-triangle packages remain available.

The October joint Leech construction regains 28D with 205 232 points,
after K. Li's 205 172 superseded our 204 896 on 10 October.
The October 7 387-line class family regains 49–55D after their September
supersession by the 7 077-line lift. The two separate record periods are
shown above; earlier constructions and their values remain in the package history.


**Matching 18D construction.** [K. Li's hyperplane-section paper](https://github.com/Felpix-Studios/kissing-number/blob/b823640f4e98134b9d55b27288132ee954a8affa/preprints/hyperplane-section-d18/paper.pdf)
also reaches 8 358. Cohn's existing 18D table entry credits this repository;
Li's later construction is a separately documented match.

**Earlier supersessions.** Our earlier values remain valid lower bounds and stay in
this repository. The later values below are from
[Cohn's leaderboard](https://cohn.mit.edu/kissing-numbers/) as retrieved on
**25 October 2026**, and, for the earlier 49–55D supersessions, from Yang et al.,
[*Large-Scale Autonomous Discovery of Kissing Number Constructions*, Table 1](https://arxiv.org/html/2609.35051v1#S1.T1)
(28 September 2026). Cohn's table covers dimensions through 48; the higher-dimensional
comparison is given separately by its source.

| dim | earlier result | intervening result | source | current result |
|---|---:|---:|---|---:|
| 25 | 197 579 | 197 583 | [Takhanov–Yun](https://github.com/k-nic/kissing_d25-31) | 197 584 |
| 29 | 209 968 | 210 495 | [K. Li](https://github.com/Felpix-Studios/kissing-number/tree/main/preprints/leech-lift-d28-d31) | 226 520 |
| 30 | 221 012 | 222 073 | [K. Li](https://github.com/Felpix-Studios/kissing-number/tree/main/preprints/leech-lift-d28-d31) | 238 808 |
| 31 | 238 662 | 241 079 | [K. Li](https://github.com/Felpix-Studios/kissing-number/tree/main/preprints/leech-lift-d28-d31) | 246 638 |
| 38 | 591 612 | 605 322 | [Q. Lin](https://github.com/QH-Lin/kissing-numbers-33-39) | 610 394 |
| 39 | 756 116 | 768 916 | [Q. Lin](https://github.com/QH-Lin/kissing-numbers-33-39) | 769 180 |

The retained 29–30D and 38–39D certificates exceed those intervening records;
PackingStar subsequently reports larger 26–31D values in its 26 October update.
The 26 October collective and dual-centre constructions regain 25D and 31D. The table above
preserves the supersession history.

The linked sources give the construction behind each displayed later result:
Takhanov–Yun use new heads with coordinated bulk motion in 25D, a deformed second
cap layer with free heads in 26D, and graph colouring for a larger second cap layer
in 27D. K. Li uses lifted 528-vector Leech blocks with auxiliary layer points in
28–31D. Q. Lin adds 6855 equatorial lines to our 38D base and uses a larger
weight-eight support code in the 39D ERS construction. Yang et al. enlarge the
P48 mother class from 7069 to 7077 lines and select its isometric images for
49–55D. These are the methods for the latest values above; the earlier historical
periods end at the first superseding publication.

These dated sources determine the supersession notes; a new global record count
is not inferred from the repository's older comparison cache.
The October 41–46D constructions are presented alongside the earlier results above,
with the supplied partial-family 43D and coordinated-motion 45D certificates.

**Historical recovered cross-sections in dimensions 46 and 47.** τ(47) ≥ 23 766 960 is due to
Boyvalenkov–Cherkashin, *Results in Mathematics* **80** (2025), Paper No. 3, and is
equation (3) of the preprint [arXiv:2312.05121](https://arxiv.org/abs/2312.05121)
(December 2023), which observes *in the same paragraph* that those A₀ vectors form a
47-dimensional kissing configuration.
τ(46) ≥ 12 309 600 is Table 3 of Sun–Wang,
[arXiv:2607.20359v3](https://arxiv.org/abs/2607.20359) (18 August 2026), and its ingredients
were in Ozeki's 2016 Siegel theta tables a decade earlier, though not stated as a kissing
number. Both were obtained here independently, and both are **recovered rather than
claimed**: they are not in `RESULTS.md`, not counted among the 70, and not presented as
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

**Historical calibration.** The 19D construction studied here was bounded by
11 948, matching Ho's record at that time. K. Li's later construction now gives
12 460 in dimension 19 (and 31 560 in dimension 21); the earlier ceiling concerns
only the old construction family. Dimension 12's record likewise landed elsewhere
while this work was in progress.

**Registered comparison scope.** Dimensions 18, 25–46, 49–71, 73–96.
These are the registered construction results compared with their recorded baselines;
the October doubled-support construction raises 96D above the retained ERS value. The list includes the superseded
results above. An earlier improvement remains a constructive result after a larger
configuration is found. The literature status of the 70D and 71D consequences of the
one-point distribution of Γ₇₂'s minimal shell still merits specialist review.

## The construction ideas

**Coxeter–Todd phase codes (36–37).** Three rank-12 blocks carry a fixed
403 704-point base and 32 656 selected three-point phase orbits. Exact shell
enumeration and integer angles prove squared distance at least 12. Two latitude
layers add 3 816 points in 37D. The earlier full-cell and phase constructions are retained
in the package history; the new latitude layer uses partial three-point orbits.

**Independently normalized body families and heads (41–46).** The October constructions
use a rank-40 complement of sqrt(3) E8 in P48n, with heads in one to six dimensions.
Each family is normalized separately, so families can have different original lengths.
Joint class replacements and head motions admit additional complete or partial families;
zero-body points in the existing head span give the latest 41D, 42D, 45D and 46D
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

**The cap construction (25–27, 29–31, 38, 49–63, 73–95).** In ℝⁿ ⊕ ℝᵏ, put a lattice shell on the
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
python run_all.py          # 88 scripts; earlier suite about 40 minutes, plus the October certificate replays, one verdict
python run_all.py --changed   # skip jobs whose inputs have not moved since they passed
python run_all.py --only dim31 # just the jobs whose label or command matches
python run_all.py --full   # 107 scripts, budget about 10 hours: adds the all-pairs sweep, the
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
| `.../dim36-37-coxeter-todd/verify.py` | exact 501 672-point 36D construction and 505 488-point 37D extension | ~15 s |
| `.../dim41-46-p48-heads/verify.py` | current 41–46D totals, partial 43D points and coordinated 45D; inherited populations | ~30 s |
| `…/dim46-47-p48-cross-sections/derive.py`                                     | recovers 12 309 600 and 23 766 960, claiming neither           | 6 s    |
| `…/dim49-63-p48-caps/verify.py` + `final.py`                                  |  the P₄₈ lattice, class, and fifteen numbers                    | 3 s    |
| `…/dim49-63-p48-caps/scripts/regenerate.py`                                   | rebuilds and re-verifies all 1400 disjoint classes             | 15 min |
| `…/dim49-63-p48-caps/scripts/ceiling.py`                                      | ceilings restricted to the shipped base classes, against level 0        | 1 s    |
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
main table reports 66 of the 70 registered claims. Dimensions 32–35 await the
next manuscript update, as do the improved values in 26–31D, 36–37D, 43D, 95D and 96D. Every file
under `paper/` remains unchanged; [the manuscript checkpoint](common/manuscript_scope.json)
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

