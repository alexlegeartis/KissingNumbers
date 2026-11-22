# Kissing numbers: constructions and verification in dimensions 18 through 96

[Henry Cohn’s kissing-number table](https://cohn.mit.edu/kissing-numbers/) credited this repository in **23 dimensions, 18 and 25–46**, in the 17 November 2026 snapshot. [Dated comparison and attribution](common/cohn_20261117.json).

Verified construction lower bounds in **70 dimensions**, with exact proof packages and retained research history. A kissing configuration consists of unit vectors whose pairwise inner products are at most 1/2.

> [!IMPORTANT]
> These are lower bounds, not proofs of optimality. Each verification package identifies its exact checks and inherited lattice, population or design premises; some large constructions are existence proofs rather than complete coordinate lists. Comparisons are limited to the named, dated sources, not a worldwide-priority guarantee. The [manuscript update remains deferred](RESEARCH_GUIDE.md#the-write-up).

[Weekly history](#recent-updates) · [Construction and verification guide](RESEARCH_GUIDE.md) · [Replay coverage](SOTA_VERIFY.md) · [Sphere-packing companion](https://github.com/alexlegeartis/SpherePackings)

The twelve-fiber continuation gives **12,884,901,888 in 88D** (+34.34659% over the preceding repository result), also valid in 89D by embedding. Mixed latitudes then raise 89D to **13,958,643,714** (about +8.33333% over that embedding).

[Re-marked axial Pless and Steiner squares](verifications/improved/dim86-87-remarked-square/) yield **8,008,545,838 in 86D** and **8,160,981,972 in 87D**. A [two-coordinate shortening in 88D](verifications/improved/dim88-axial-pless/) gives **13,357,769,652** (+2.52%).

Two asymmetric latitude pairs and [353 new Golay support patterns](verifications/improved/dim88-89-golay-variable-punctures/) raise **89D to 15,177,390,082** (+1,125,793,792 over the prior repository construction). The same new patterns also improve **88D to 13,029,906,432**.

Larger [shifted dodecacode lists](verifications/improved/dim92-96-shifted-dodecacode-outer-extension/) and [BCH and complementary parity flat layers](verifications/improved/dim94-96-hexary-parity-packing/) are retained in the [adapted dodecad sign kernels](verifications/improved/dim96-adapted-dodecad-fibres/), and the [mixed-tetrad axis and sign fibres](verifications/improved/dim96-mixed-tetrad-axis-sign-fibres/) and [row-alphabet repacking](verifications/improved/dim96-mixed-tetrad-row-alphabet-repack/) bring the construction to **209,299,907,072 in 96D**. The separate [108-head E7 construction](verifications/improved/dim95-e7-108-heads/) with [coordinate-aligned dodecad fibres](verifications/improved/dim94-95-coordinate-aligned-dodecad-fibres/) and [mixed-tetrad fibres](verifications/improved/dim94-95-mixed-tetrad-axis-sign-fibres/) reaches **117,097,130,108 in 95D**.

Individually normalized Z8 fibers and the [enlarged flat layer](verifications/improved/dim91-93-flat-parity-extensions/) reach **25,917,421,568 in 91D**, **34,507,356,160 in 92D**, and **47,392,258,048 in 93D**.

The [shortening and larger pattern lists](verifications/improved/dim88-95-golay-shortening/) extend the Golay–dodecacode method to **88D and 89D**, and further strengthen **92–95D**.

## Results and construction history

**Recorded comparison** is the historical, attributed baseline, not necessarily the latest competing result. **Source years** date the cited results or the ingredients of a derived bound, not a standalone publication of the derived total. Earlier dates are retained when a later linked revision contains the same value. **Gain** is measured against that baseline; later improvements and supersessions are kept in the [comparison history](HISTORY.md#record-status-and-retained-results). **SOTA** means best in the reviewed sources, not proven optimal. Cohn’s live table was checked on 17 November 2026; its 1–48D coverage does not establish the status of higher-dimensional results, which use the separately named sources. [Full column definitions and comparison scope](RESEARCH_GUIDE.md#reading-the-results-table).

| dim | recorded comparison | comparison method / source | **this work** | SOTA value | gain over recorded comparison | later gain | SOTA period | our method | superseding method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 18 | 7 654 | [Cohn–Li (2024; v3 2026), BW16 + hexagons](https://arxiv.org/abs/2411.04916v3) | **8 358** | **SOTA** | +9.2% | — | since 2026-09-17 | [BW16](verifications/improved/dim18-bent-hexagon/) | — |
| 25 | 197 056 | [Ma et al., PackingStar (2025; v4 2026)](https://arxiv.org/abs/2511.13391v4) | **197 586** | **SOTA** | +0.269% | — |2026-08-30 to 2026-09-28; since 2026-10-26 | [Critical-cluster replacement and exact contacts](verifications/improved/dim25-critical-cluster/) | — |
| 26 | 199 868 | [Takhanov–Yun, deformed caps + free heads (2026)](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **201 204** | **SOTA** | +0.66844% | — | updated 2026-11-20 | [27D section and rational kernel](verifications/improved/dim26-section-rational-kernel/) | — |
| 27 | 201 589 | [Takhanov–Yun, cap colouring (2026)](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **209 884** | 210 108 | +4.1148% | +0.10673% | superseded 2026-11-17 | [E8 head-frame replacement](verifications/improved/dim26-27-radical-orbits/) | [Takhanov’s Leech midpoint construction (2026)](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1) |
| 28 | 205 172 | [K. Li, Leech lift and auxiliary layers (2026)](https://github.com/Felpix-Studios/kissing-number/blob/b823640f4e98134b9d55b27288132ee954a8affa/preprints/leech-lift-d28-d31/paper.pdf) | **223 528** | **SOTA** | +8.94664% | — | updated 2026-11-18 | [E8 root and auxiliary exchange](verifications/improved/dim28-takhanov-root-exchange/) | — |
| 29 | 209 496 | [Ma et al., PackingStar (2025; v4 2026)](https://arxiv.org/abs/2511.13391v4) | **232 157** | **SOTA** | +10.81691% | — | updated 2026-11-19 | [14-for-11 even-grid exchange](verifications/improved/dim29-even-grid-exchange/) | — |
| 30 | 220 440 | [Ma et al., PackingStar (2025; v5 2026)](https://arxiv.org/abs/2511.13391v5) | **248 506** | **SOTA** | +12.73181% | — | updated 2026-11-22 | [Lifted Leech contact midpoints](verifications/improved/dim30-row-midpoint-packing/) | — |
| 31 | 238 350 | [Ma et al., PackingStar (2026; v5)](https://arxiv.org/abs/2511.13391v5) | **251 524** | **SOTA** | +5.52717% | — | updated 2026-11-18 | [Joint midpoint exchange](verifications/improved/dim30-31-midpoint-exchanges/) | — |
| 32 | 352 320 | [Lin, support code + ERS layers (2026)](https://github.com/QH-Lin/kissing-number-32/tree/f367a9c5a051c0d6821ce71324b9dfa1ead76c2e) | **352 992** | **SOTA** | +0.190736% | — | since 2026-11-07 | [Non-antipodal signed octad continuation](verifications/improved/dim32-35-signed-layers/nonantipodal32/) | [Comparison scope](verifications/improved/dim32-35-signed-layers/nonantipodal32/evidence/public_comparison.json) |
| 33 | 370 240 | [Lin, support code + ERS layers (2026)](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **370 528** | **SOTA** | +0.077787% | — | updated 2026-11-22 | [Joint signed refill](verifications/improved/dim33-joint-signed-refill/) | [Comparison scope](verifications/improved/dim32-35-signed-layers/nonantipodal33/evidence/public_comparison.json) |
| 34 | 390 980 | [Lin, support code + ERS layers (2026)](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **391 004** | **SOTA** | +0.00614% | — | updated 2026-11-20 | [Octad swaps and signed stars](verifications/improved/dim34-octad-switches/) | — |
| 35 | 415 308 | [Lin, support code + ERS layers (2026)](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **415 360** | **SOTA** | +0.01252% | — | updated 2026-11-21 | [Joint signed-star and root repacking](verifications/improved/dim35-joint-star-root/) | — |
| 36 | 484 760 | [Yang et al., signed replacements (2026)](https://arxiv.org/abs/2609.35051v1) | **502 902** | **SOTA** | +3.7425% | — | updated 2026-11-22 | [Coxeter–Todd triangle switches and phase refill](verifications/improved/dim36-ct-triangle-refill/) | — |
| 37 | 500 008 | [Olivier, using Lin's code (2026)](https://doi.org/10.5281/zenodo.23569153) | **544 666** | **SOTA** | +8.9315% | — | updated 2026-11-22 | [QR 12-flat phase repacking](verifications/improved/dim37-qr-flat-repack/) | — |
| 38 | 605 322 | [Qihong Lin, Leech equatorial augmentation (2026)](https://arxiv.org/abs/2611.03682) | **652 076** | **SOTA** | +7.7238% | — | updated 2026-11-20 | [QR offsets, octad switches and signed stars](verifications/improved/dim38-qr-octad-switches/) | — |
| 39 | 755 988 | [Edel–Rains–Sloane (1998 method)](https://arxiv.org/abs/math/0207291) + [Brouwer’s archived 39D bound](https://aeb.win.tue.nl/codes/Andw.html-old5#kiss) | **783 468** | **SOTA** | +3.63% | — | updated 2026-11-20 | [Three-coordinate QR heads](verifications/improved/dim39-qr-three-coordinate-heads/) | — |
| 40 | 1 092 000 | [Dorofeev–Sun–Wang, P48 sections (2026)](https://arxiv.org/html/2607.20359v4) | **1 126 351** | **SOTA** | +3.15% | — | updated 2026-11-22 | [P48q affine refill through empty difference classes](verifications/improved/dim40-p48q-affine-refill/) | — |
| 41 | 1 324 472 | [Dorofeev–Sun–Wang, P48 sections (2026)](https://arxiv.org/abs/2607.20359v4) | **1 460 036** | **SOTA** | +10% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 42 | 1 792 386 | [Dorofeev–Sun–Wang, P48 sections (2026)](https://arxiv.org/abs/2607.20359v4) | **1 991 445** | **SOTA** | +11% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 43 | 2 553 792 | [Yang et al., P48n section (2026)](https://arxiv.org/abs/2609.35051v1) | **2 919 803** | **SOTA** | +14% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 44 | 4 153 600 | [Conway–Sloane P48 sections (1982)](https://doi.org/10.2307/2007025); [counts: Dorofeev–Sun–Wang (2026)](https://arxiv.org/abs/2607.20359v4) | **4 939 281** | **SOTA** | +19% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 45 | 7 380 090 | [Yang et al., QR neighbour projection (2026)](https://arxiv.org/abs/2609.35051v1) | **7 655 056** | **SOTA** | +3.7% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 46 | 12 309 600 | [Conway–Sloane P48 sections (1982)](https://doi.org/10.2307/2007025); [counts: Dorofeev–Sun–Wang (2026)](https://arxiv.org/abs/2607.20359v4) | **13 387 797** | **SOTA** | +8.8% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 49 | 52 416 002 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [1D tail](common/published.py) *(derived)* | **52 431 624** | **SOTA** | +0.0298% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 50 | 52 416 006 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [2D tail](common/published.py) *(derived)* | **52 462 872** | **SOTA** | +0.0894% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 51 | 52 416 012 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [3D tail](common/published.py) *(derived)* | **52 509 740** | **SOTA** | +0.1788% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 52 | 52 416 024 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [4D tail](common/published.py) *(derived)* | **52 603 402** | **SOTA** | +0.3575% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 53 | 52 416 040 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [5D tail](common/published.py) *(derived)* | **52 728 054** | **SOTA** | +0.5953% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 54 | 52 416 072 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [6D tail](common/published.py) *(derived)* | **52 976 623** | **SOTA** | +1.0694% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 55 | 52 416 126 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [7D tail](common/published.py) *(derived)* | **53 393 546** | **SOTA** | +1.8647% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 56 | 52 416 240 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [8D tail](common/published.py) *(derived)* | **54 263 398** | **SOTA** | +3.5240% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 57 | 52 416 306 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [9D tail: Leech–Sloane 1971](https://doi.org/10.4153/CJM-1971-081-3) *(derived)* | **54 760 693** | **SOTA** | +4.4726% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 58 | 52 416 510 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [10D tail: Ganzhinov 2022](https://arxiv.org/abs/2207.08266) *(derived)* | **56 270 076** | **SOTA** | +7.3518% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 59 | 52 416 604 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [11D tail: Bianchi et al. 2026](https://arxiv.org/abs/2606.10402) *(derived)* | **56 952 414** | **SOTA** | +8.6534% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 60 | 52 416 841 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [12D tail: Takhanov et al. 2026](https://arxiv.org/abs/2606.18984) *(derived)* | **58 631 571** | **SOTA** | +11.8564% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 61 | 52 417 154 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [13D tail: Zinoviev–Ericson 1999](https://www.mathnet.ru/eng/ppi457) *(derived)* | **91 038 762** | **SOTA** | +73.681238% | — | updated 2026-11-22 | [Omitted-fibre weighted spikes](verifications/improved/dim61-62-omitted-fibre-spikes/) | — |
| 62 | 52 417 932 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [14D tail: Ganzhinov 2022](https://arxiv.org/abs/2207.08266) *(derived)* | **126 571 579** | **SOTA** | +141.466181% | — | updated 2026-11-22 | [Partial nonary signed refill](verifications/improved/dim62-partial-nonary-refill/) | — |
| 63 | 52 418 564 | [Leech–Sloane P48 (1971)](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [15D tail: Leech–Sloane 1971](https://doi.org/10.4153/CJM-1971-081-3) *(derived)* | **184 428 848** | **SOTA** | +251.838803% | — | updated 2026-11-22 | [Mixed-depth dense owner clones](verifications/improved/dim63-64-mixed-depth-clones/) | — |
| 64 | 331 737 984 | [Edel-Rains-Sloane (1998)](https://arxiv.org/abs/math/0207291) | **335 798 944** | **SOTA** | +1.224147% | — | updated 2026-11-22 | [Mixed-depth dense owner clones](verifications/improved/dim63-64-mixed-depth-clones/) | — |
| 65 | 331 737 986 | [ERS64 (1998)](https://arxiv.org/abs/math/0207291) plus two poles *(derived)* | **336 258 432** | **SOTA** | +1.362656% | — | updated 2026-11-22 | [Mixed-depth clone transfer](verifications/improved/dim65-mixed-depth-transfer/) | — |
| 66 | 331 737 990 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291) ⊕ A2 *(derived)* | **342 390 951** | **SOTA** | +3.211257% | — | updated 2026-11-22 | [Mixed-depth clone transfer](verifications/improved/dim66-67-mixed-depth-transfer/) | — |
| 67 | 331 737 996 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291) ⊕ A3 *(derived)* | **355 203 873** | **SOTA** | +7.073618% | — | updated 2026-11-22 | [Mixed-depth clone transfer](verifications/improved/dim66-67-mixed-depth-transfer/) | — |
| 68 | 331 737 984 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **472 598 888** | **SOTA** | +42.461494% | — | updated 2026-11-17 | [Exact Gamma72 D4 fiber rescaling](verifications/improved/dim68-91-lattice-code-continuation/) | — |
| 69 | 331 737 984 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **780 388 888** | **SOTA** | +135.242549% | — | updated 2026-11-17 | [Exact Gamma72 D4 fiber rescaling](verifications/improved/dim68-91-lattice-code-continuation/) | — |
| 70 | 331 737 984 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **1 361 151 808** | **SOTA** | 4.10x | — | updated 2026-11-22 | [Mixed Γ72 A2 normal embedding](verifications/improved/dim70-gamma72-mixed-a2-normals/) | — |
| 71 | 331 737 984 | [Edel–Rains–Sloane 64D (1998)](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **2 623 656 288** | **SOTA** | +690.882086% | — | updated 2026-11-22 | [Single-block Γ72 refill](verifications/improved/dim71-gamma72-single-block-refill/) | — |
| 73 | 6 218 175 602 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [1D tail](common/published.py) *(derived)* | **6 218 449 010** | **SOTA** | +0.004396% | — | updated 2026-11-22 | [Depth-six Γ72 lattice shell](verifications/improved/dim73-82-gamma72-depth6-lattice/) | — |
| 74 | 6 218 175 606 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [2D tail](common/published.py) *(derived)* | **6 218 995 830** | **SOTA** | +0.013190% | — | updated 2026-11-22 | [Depth-six Γ72 lattice shell](verifications/improved/dim73-82-gamma72-depth6-lattice/) | — |
| 75 | 6 218 175 612 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [3D tail](common/published.py) *(derived)* | **6 219 816 924** | **SOTA** | +0.026395% | — | updated 2026-11-22 | [Depth-six Γ72 lattice shell](verifications/improved/dim73-82-gamma72-depth6-lattice/) | — |
| 76 | 6 218 175 624 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [4D tail](common/published.py) *(derived)* | **6 221 459 976** | **SOTA** | +0.052818% | — | updated 2026-11-22 | [Depth-six Γ72 lattice shell](verifications/improved/dim73-82-gamma72-depth6-lattice/) | — |
| 77 | 6 218 175 640 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [5D tail](common/published.py) *(derived)* | **6 222 013 132** | **SOTA** | +0.061714% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 78 | 6 218 175 672 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [6D tail](common/published.py) *(derived)* | **6 223 119 444** | **SOTA** | +0.079505% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 79 | 6 218 175 726 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [7D tail](common/published.py) *(derived)* | **6 224 780 712** | **SOTA** | +0.106220% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 80 | 6 218 175 840 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ E8 (tabulated direct sum) | **6 228 104 952** | **SOTA** | +0.159678% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 81 | 6 218 175 906 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [9D tail: Leech–Sloane 1971](https://doi.org/10.4153/CJM-1971-081-3) *(derived)* | **6 228 378 586** | **SOTA** | +0.164078% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 82 | 6 218 176 110 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [10D tail: Ganzhinov 2022](https://arxiv.org/abs/2207.08266) *(derived)* | **6 228 925 854** | **SOTA** | +0.172876% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 83 | 6 218 176 204 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [11D tail: Bianchi et al. 2026](https://arxiv.org/abs/2606.10402) *(derived)* | **6 229 747 428** | **SOTA** | +0.186087% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 84 | 6 218 176 441 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [12D tail: Takhanov et al. 2026](https://arxiv.org/abs/2606.18984) *(derived)* | **6 231 391 584** | **SOTA** | +0.212524% | — | updated 2026-11-22 | [Half-root Γ72 lattice shell](verifications/improved/dim77-84-gamma72-half-root-shells/) | — |
| 85 | 6 218 176 754 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [13D tail: Zinoviev–Ericson 1999](https://www.mathnet.ru/eng/ppi457) *(derived)* | **6 235 993 850** | **SOTA** | +0.286532% | — | updated 2026-11-22 | [Finer exact Γ72 recurrence](verifications/improved/dim74-85-gamma72-fine-recurrence/) | — |
| 86 | 6 218 177 532 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [14D tail: Ganzhinov 2022](https://arxiv.org/abs/2207.08266) *(derived)* | **8 008 552 282** | **SOTA** | +28.792596% | — | updated 2026-11-22 | [Rotated plane trades](verifications/improved/dim86-87-90-rotated-plane-trades/) | — |
| 87 | 6 218 178 164 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [15D tail: Leech–Sloane 1971](https://doi.org/10.4153/CJM-1971-081-3) *(derived)* | **8 160 988 416** | **SOTA** | +31.244043% | — | updated 2026-11-22 | [Rotated plane trades](verifications/improved/dim86-87-90-rotated-plane-trades/) | — |
| 88 | 6 218 179 920 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [16D tail: Barnes–Wall 1959](https://doi.org/10.1017/S1446788700025064) *(derived)* | **13 358 363 692** | **SOTA** | +114.827552% | — | updated 2026-11-22 | [Full Pless census and new chord marks](verifications/improved/dim88-full-pless-chord-marks/) | — |
| 89 | 6 218 181 330 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [17D tail: Cohn–Li 2024](https://arxiv.org/abs/2411.04916v1) *(derived)* | **15 185 856 402** | **SOTA** | +144.217008% | — | updated 2026-11-22 | [Restricted Reed–Solomon axis code](verifications/improved/dim89-96-rs-axis-code-repack/) | — |
| 90 | 6 218 183 254 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [18D tail: Cohn–Li 2024](https://arxiv.org/abs/2411.04916v1) *(derived)* | **21 944 111 386** | **SOTA** | +252.902295% | — | updated 2026-11-22 | [Rotated plane trades](verifications/improved/dim86-87-90-rotated-plane-trades/) | — |
| 91 | 6 218 187 548 | [Nebe Γ72 (2010 preprint; 2012 journal)](https://arxiv.org/abs/1008.2862v3) ⊕ [19D tail: Ho 2026](https://arxiv.org/abs/2603.10425) *(derived)* | **25 923 338 000** | **SOTA** | +316.895402% | — | updated 2026-11-22 | [Non-Steiner triality append](verifications/improved/dim91-93-nonsteiner-triality-append/) | — |
| 92 | 8 589 934 592 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code (2001)](https://doi.org/10.1109/18.923742) *(shortened)* | **34 513 272 592** | **SOTA** | +301.787374% | — | updated 2026-11-22 | [Non-Steiner triality append](verifications/improved/dim91-93-nonsteiner-triality-append/) | — |
| 93 | 17 179 869 184 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code (2001)](https://doi.org/10.1109/18.923742) *(shortened)* | **49 545 659 928** | **SOTA** | +188.393697% | — | updated 2026-11-22 | [Bipartite D5 heads and 46 Z8 fibres](verifications/improved/dim93-z8-split-d5-heads/) | — |
| 94 | 34 359 738 368 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code (2001)](https://doi.org/10.1109/18.923742) *(shortened)* | **74 147 456 772** | **SOTA** | +115.797501% | — | updated 2026-11-22 | [Mixed-tetrad axis and sign fibres](verifications/improved/dim94-95-mixed-tetrad-axis-sign-fibres/) | — |
| 95 | 68 719 476 736 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code (2001)](https://doi.org/10.1109/18.923742) *(shortened)* | **117 097 130 108** | **SOTA** | +70.398751% | — | updated 2026-11-22 | [Mixed-tetrad axis and sign fibres](verifications/improved/dim94-95-mixed-tetrad-axis-sign-fibres/) | — |
| 96 | 137 438 953 472 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code (2001)](https://doi.org/10.1109/18.923742) *(derived by sign embedding)* | **209 299 907 072** | **SOTA** | +52.285725% | — | updated 2026-11-22 | [Mixed-tetrad row-alphabet repack](verifications/improved/dim96-mixed-tetrad-row-alphabet-repack/) | — |

## Recent updates

**22 November 2026 — 62D partial nonary refill.** An exact sign exchange adds 1,914 and retires 971 old rays, producing **126,571,579** (+943). [Proof and replay](verifications/improved/dim62-partial-nonary-refill/).

**22 November 2026 — 71D single-block refill.** Keeping the entire norm-ten slice, a Leech-block exchange adds 81,952 rays and removes 4,600, giving **2,623,656,288** (+77,352). [Exact proof and replay](verifications/improved/dim71-gamma72-single-block-refill/).

**22 November 2026 — Half-root Γ72 shells.** Exact lattice shells improve **eight bounds in 77–84D**; the largest gains over the preceding certificates are **3,314,392** each in 81D and 82D. [Proof and replay](verifications/improved/dim77-84-gamma72-half-root-shells/).

Expand a week for the main milestones; follow its history link for the full dated sequence, intermediate values and proof notes. Current bounds are in the table above.

Empty short-difference classes in the \(P_{48q}\) complement take **40D to 1,126,351** (+24). [Exact replay](verifications/improved/dim40-p48q-affine-refill/).

Individual Coxeter–Todd triangle switches and a three-block phase refill raise **36D to 502,902** (+81). [Exact replays](verifications/improved/dim36-ct-triangle-refill/).

<details>
<summary>16–22 November 2026 (through 19 November)</summary>

A larger radical support packing takes **37D to 544,474** (+10,610 over the original QR result). [37D proof](verifications/improved/dim37-qr-radical-families/). The [536,440-point intermediate](verifications/improved/dim37-qr-shared-stars/) is retained.

Joint repacking raises **33D to 370,528** (+112). [Exact 33D proof](verifications/improved/dim33-joint-signed-refill/).

Four octad swaps raise **34D to 391,004** (+12). [34D proof](verifications/improved/dim34-octad-switches/).

Joint signed-star and root repacking raises **35D to 415,360** (+2). [35D proof](verifications/improved/dim35-joint-star-root/).

A coordinated pair of octad switches raises **38D to 652,076** (+4). A nonlocal Coxeter–Todd phase rebuild raises **36D to 502,821** (+54). [38D proof](verifications/improved/dim38-qr-octad-switches/) · [36D proof](verifications/improved/dim36-ct-nonlocal/).

QR coordinate compression gives **533,864 in 37D** (+27,155) and **643,836 in 38D** (+32,038). [Exact proofs](verifications/improved/dim37-38-qr-merge/).

Rebuilding the Coxeter–Todd cell layer raises **36D to 502,767** (+945) and **37D to 506,709** (+1,071). [Exact proof](verifications/improved/dim36-37-ct-cell-rebuild/).

Triangular Hamming latitudes raise **38D to 611,798**, a gain of 160. [Exact proof](verifications/improved/dim38-triangle-hamming/).

A rebuilt 10D section of a 27D Leech-based code gives [**201,204 in 26D**](verifications/improved/dim26-section-rational-kernel/), 12 above [Takhanov’s 201,192](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1). The earlier alternate construction tied Takhanov at 201,192. Fixed-projection obstructions cover the 26D and 27D comparison witnesses. [26D and 27D research](research/26_27_native_projection_obstructions_20261118/)

The 29D projected-grid search gives **232,157**, 3 above the preceding repository witness and 4 above Takhanov’s pinned release. [14-for-11 even-grid exchange](verifications/improved/dim29-even-grid-exchange/).

Exact exchanges give **223,528 in 28D**, **248,506 in 30D**, and **251,524 in 31D**. [28D proof](verifications/improved/dim28-takhanov-root-exchange/) · [30D proof](verifications/improved/dim30-row-midpoint-packing/) · [31D proof](verifications/improved/dim30-31-midpoint-exchanges/).

A Gram-dominated E6 assignment strengthens 94D to 73,016,281,804 points. Direct high-dimensional layers strengthen 79–96D, while planar-head compression improves 70D. Mixed-shell sections and code constructions also advance 68–69D and 71D. Cohn’s current table explicitly credits this repository in 18D and 25–46D.

[Full weekly history and verification links](HISTORY.md#week-2026-11-16).

</details>

<details>
<summary>9–15 November 2026</summary>

Exact exchanges and certificate unions advance 26–31D; affine sign families and Gamma72 sections strengthen 61–69D. Further 25D, 36–37D and 49–60D refinements retain their independent proof packages.

[Full weekly history and verification links](HISTORY.md#week-2026-11-09).

</details>

<details>
<summary>2–8 November 2026</summary>

Algebraic completions, latitude exchanges, Coxeter–Todd switches and pointwise P48 refills improve several lower- and higher-dimensional constructions.

[Full weekly history and verification links](HISTORY.md#week-2026-11-02).

</details>

<details>
<summary>26 October–1 November 2026</summary>

Joint replacements, shared parity repairs and nonlinear-code sections broaden the construction families. The nonlinear parent attribution is corrected explicitly; its published contribution is separated from our added layers.

[Full weekly history and verification links](HISTORY.md#week-2026-10-26).

</details>

<details>
<summary>19–25 October 2026</summary>

Midpoint replacements, joint Leech layers, phase codes and nearest-coset constructions establish the earlier checkpoints. Their superseded values remain part of the research record.

[Full weekly history and verification links](HISTORY.md#week-2026-10-19).

</details>

## Reproducing the results

```bash
python audit.py --repository-only  # registry and repository consistency
python run_all.py --sota           # current construction replay suite
```

These commands have different verification scopes. See [coverage and inherited inputs](SOTA_VERIFY.md), [full reproduction instructions](RESEARCH_GUIDE.md#reproducing-everything), and [requirements](RESEARCH_GUIDE.md#reproducing-everything). Use `--skip-paper` for the documented deferred-manuscript workflow; do not interpret a quick audit as a complete geometric replay.

## Licence, citation, and how this was made

MIT (see `[LICENSE](LICENSE)`). `[CITATION.cff](CITATION.cff)` carries the metadata GitHub
uses for its "Cite this repository" button.

The constructions, the verification code and the write-ups in this repository were developed
by the author working with **Claude Code** (Anthropic), **ChatGPT Pro (GPT 6)**, and **Astra & Sol 6.1 in Codex**. The repository is deliberately
structured so that someone else can pick it up the same way: every package states its idea in
prose before its code, every negative result says what was tried and why it failed, and
`KNOWLEDGE.md` is the full working record — 160 sections, most of them about things that
did not work.

Alexey Kravatskiy, MIRIAI (Moscow Independent Research Institute of Artificial Intelligence),
Moscow, Russia. [kravatskii.a@miriai.org](mailto:kravatskii.a@miriai.org)
