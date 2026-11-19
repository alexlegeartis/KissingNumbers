# Kissing numbers: constructions and verification in dimensions 18 through 96

[Henry Cohn’s kissing-number table](https://cohn.mit.edu/kissing-numbers/) credited this repository in **23 dimensions, 18 and 25–46**, in the 17 November 2026 snapshot. [Dated comparison and attribution](common/cohn_20261117.json).

Verified construction lower bounds in **70 dimensions**, with exact proof packages and retained research history. A kissing configuration consists of unit vectors whose pairwise inner products are at most 1/2.

> [!IMPORTANT]
> These are lower bounds, not proofs of optimality. Each verification package identifies its exact checks and inherited lattice, population or design premises; some large constructions are existence proofs rather than complete coordinate lists. Comparisons are limited to the named, dated sources, not a worldwide-priority guarantee. The [manuscript update remains deferred](RESEARCH_GUIDE.md#the-write-up).

[Weekly history](#recent-updates) · [Construction and verification guide](RESEARCH_GUIDE.md) · [Replay coverage](SOTA_VERIFY.md) · [Sphere-packing companion](https://github.com/alexlegeartis/SpherePackings)

The twelve-fiber continuation gives **12,884,901,888 in 88D** (+34.34659% over the preceding repository result), also valid in 89D by embedding. Mixed latitudes then raise 89D to **13,958,643,714** (about +8.33333% over that embedding).

## Results and construction history

**Recorded comparison** is the historical, attributed baseline, not necessarily the latest competing result. **Gain** is measured against that baseline; later improvements and supersessions are kept in the [comparison history](HISTORY.md#record-status-and-retained-results). **SOTA** means best in the reviewed sources, not proven optimal. Cohn’s live table was checked on 17 November 2026; its 1–48D coverage does not establish the status of higher-dimensional results, which use the separately named sources. [Full column definitions and comparison scope](RESEARCH_GUIDE.md#reading-the-results-table).

| dim | recorded comparison | comparison method / source | **this work** | SOTA value | gain over recorded comparison | later gain | SOTA period | our method | superseding method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 18 | 7 654 | [Cohn–Li, BW16 + hexagons](https://arxiv.org/abs/2411.04916v3) | **8 358** | **SOTA** | +9.2% | — | since 2026-09-17 | [BW16](verifications/improved/dim18-bent-hexagon/) | — |
| 25 | 197 056 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **197 586** | **SOTA** | +0.269% | — |2026-08-30 to 2026-09-28; since 2026-10-26 | [Critical-cluster replacement and exact contacts](verifications/improved/dim25-critical-cluster/) | — |
| 26 | 199 868 | [Takhanov–Yun, deformed caps + free heads](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **201 188** | 201 192 | +0.6604% | +0.00199% | superseded 2026-11-17 | [Radical-orbit continuation](verifications/improved/dim26-27-radical-orbits/) | [Takhanov’s Leech midpoint construction](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1) |
| 27 | 201 589 | [Takhanov–Yun, cap colouring](https://github.com/k-nic/kissing_d25-31/blob/b3768310ba566da5c8947546266c103fbcdce181/paper.pdf) | **209 884** | 210 108 | +4.1148% | +0.10673% | superseded 2026-11-17 | [E8 head-frame replacement](verifications/improved/dim26-27-radical-orbits/) | [Takhanov’s Leech midpoint construction](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1) |
| 28 | 205 172 | [K. Li, Leech lift and auxiliary layers](https://github.com/Felpix-Studios/kissing-number/blob/b823640f4e98134b9d55b27288132ee954a8affa/preprints/leech-lift-d28-d31/paper.pdf) | **223 528** | **SOTA** | +8.94664% | — | updated 2026-11-18 | [E8 root and auxiliary exchange](verifications/improved/dim28-takhanov-root-exchange/) | — |
| 29 | 209 496 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v4) | **232 157** | **SOTA** | +10.81691% | — | updated 2026-11-19 | [14-for-11 even-grid exchange](verifications/improved/dim29-even-grid-exchange/) | — |
| 30 | 220 440 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v5) | **248 416** | **SOTA** | +12.69098% | — | updated 2026-11-18 | [Reflected certificate unions](verifications/improved/dim30-reflected-certificate-unions/) | — |
| 31 | 238 350 | [Ma et al., PackingStar](https://arxiv.org/abs/2511.13391v5) | **251 524** | **SOTA** | +5.52717% | — | updated 2026-11-18 | [Joint midpoint exchange](verifications/improved/dim30-31-midpoint-exchanges/) | — |
| 32 | 352 320 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-number-32/tree/f367a9c5a051c0d6821ce71324b9dfa1ead76c2e) | **352 992** | **SOTA** | +0.190736% | — | since 2026-11-07 | [Non-antipodal signed octad continuation](verifications/improved/dim32-35-signed-layers/nonantipodal32/) | [Comparison scope](verifications/improved/dim32-35-signed-layers/nonantipodal32/evidence/public_comparison.json) |
| 33 | 370 240 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **370 408** | **SOTA** | +0.045376% | — | since 2026-11-07 | [Non-antipodal octads and individual root-star trades](verifications/improved/dim32-35-signed-layers/nonantipodal33/) | [Comparison scope](verifications/improved/dim32-35-signed-layers/nonantipodal33/evidence/public_comparison.json) |
| 34 | 390 980 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **390 992** | **SOTA** | +0.00306921% | — | since 2026-11-06 | [Rational heavy-head rotation](verifications/improved/dim32-35-signed-layers/rotated34/) | — |
| 35 | 415 308 | [Lin, support code + ERS layers](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) | **415 336** | **SOTA** | +0.00674198% | — | updated 2026-11-19 | [Shared root deletions across heads and stars](verifications/improved/dim35-shared-root-heads/) | — |
| 36 | 484 760 | [Yang et al., signed replacements](https://arxiv.org/abs/2609.35051v1) | **502 821** | **SOTA** | +3.7258% | — | updated 2026-11-19 | [Nonlocal Coxeter–Todd phase rebuild](verifications/improved/dim36-ct-nonlocal/) | — |
| 37 | 500 008 | [Olivier, using Lin's code](https://doi.org/10.5281/zenodo.23569153) | **544 280** | **SOTA** | +8.8543% | — | updated 2026-11-19 | [QR compression and radical support packing](verifications/improved/dim37-qr-radical-families/) | — |
| 38 | 605 322 | [Qihong Lin, Leech equatorial augmentation](https://arxiv.org/abs/2611.03682) | **652 028** | **SOTA** | +7.7159% | — | updated 2026-11-19 | [New QR offsets after compression](verifications/improved/dim38-qr-free-offsets/) | — |
| 39 | 755 988 | [Edel–Rains–Sloane](https://arxiv.org/abs/math/0207291) + [Brouwer codes](https://aeb.win.tue.nl/codes/Andw.html) | **769 852** | **SOTA** | +1.83% | — | 2026-08-30 to 2026-09-28; since 2026-10-25 | [Reduced-symmetry sign orbits](verifications/improved/dim39-deformed-orbits/) | — |
| 40 | 1 092 000 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/html/2607.20359v4) | **1 126 314** | **SOTA** | +3.14% | — | since 2026-10-26 | [Affine lattice shells](verifications/improved/dim40-affine-shells/) | — |
| 41 | 1 324 472 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/abs/2607.20359v4) | **1 460 036** | **SOTA** | +10% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 42 | 1 792 386 | [Dorofeev–Sun–Wang, P48 sections](https://arxiv.org/abs/2607.20359v4) | **1 991 445** | **SOTA** | +11% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 43 | 2 553 792 | [Yang et al., P48n section](https://arxiv.org/abs/2609.35051v1) | **2 919 803** | **SOTA** | +14% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 44 | 4 153 600 | Conway–Sloane P48 sections; [counts: Dorofeev–Sun–Wang](https://arxiv.org/abs/2607.20359v4) | **4 939 281** | **SOTA** | +19% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 45 | 7 380 090 | [Yang et al., QR neighbour projection](https://arxiv.org/abs/2609.35051v1) | **7 655 056** | **SOTA** | +3.7% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 46 | 12 309 600 | Conway–Sloane P48 sections; [counts: Dorofeev–Sun–Wang](https://arxiv.org/abs/2607.20359v4) | **13 387 797** | **SOTA** | +8.8% | — | since 2026-10-23 | [P48n heads](verifications/improved/dim41-46-p48-heads/) | — |
| 49 | 52 416 002 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [1D tail](common/published.py) *(derived)* | **52 431 624** | **SOTA** | +0.0298% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 50 | 52 416 006 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [2D tail](common/published.py) *(derived)* | **52 462 872** | **SOTA** | +0.0894% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 51 | 52 416 012 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [3D tail](common/published.py) *(derived)* | **52 509 740** | **SOTA** | +0.1788% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 52 | 52 416 024 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [4D tail](common/published.py) *(derived)* | **52 603 402** | **SOTA** | +0.3575% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 53 | 52 416 040 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [5D tail](common/published.py) *(derived)* | **52 728 054** | **SOTA** | +0.5953% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 54 | 52 416 072 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [6D tail](common/published.py) *(derived)* | **52 976 623** | **SOTA** | +1.0694% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 55 | 52 416 126 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [7D tail](common/published.py) *(derived)* | **53 393 546** | **SOTA** | +1.8647% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 56 | 52 416 240 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [8D tail](common/published.py) *(derived)* | **54 263 398** | **SOTA** | +3.5240% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 57 | 52 416 306 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [9D tail](common/published.py) *(derived)* | **54 760 693** | **SOTA** | +4.4726% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 58 | 52 416 510 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [10D tail](common/published.py) *(derived)* | **56 270 076** | **SOTA** | +7.3518% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 59 | 52 416 604 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [11D tail](common/published.py) *(derived)* | **56 952 414** | **SOTA** | +8.6534% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 60 | 52 416 841 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [12D tail](common/published.py) *(derived)* | **58 631 571** | **SOTA** | +11.8564% | — | updated 2026-11-13 | [Signed P48 classes](verifications/improved/dim49-60-p48-signed-refill/) | — |
| 61 | 52 417 154 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [13D tail](common/published.py) *(derived)* | **88 979 562** | **SOTA** | +69.7528% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim61-62-66-affine-reconstruction/) | — |
| 62 | 52 417 932 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [14D tail](common/published.py) *(derived)* | **122 543 436** | **SOTA** | +133.7815% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim61-62-66-affine-reconstruction/) | — |
| 63 | 52 418 564 | [Leech–Sloane P48](https://doi.org/10.4153/CJM-1971-081-3) ⊕ [15D tail](common/published.py) *(derived)* | **184 415 700** | **SOTA** | +251.8137% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim63-65-affine-reconstruction/) | — |
| 64 | 331 737 984 | [Edel-Rains-Sloane](https://arxiv.org/abs/math/0207291) | **335 749 696** | **SOTA** | +1.2093% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim63-65-affine-reconstruction/) | — |
| 65 | 331 737 986 | [ERS64](https://arxiv.org/abs/math/0207291) plus two poles *(derived)* | **335 817 922** | **SOTA** | +1.2299% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim63-65-affine-reconstruction/) | — |
| 66 | 331 737 990 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291) ⊕ A2 *(derived)* | **342 291 652** | **SOTA** | +3.1813% | — | updated 2026-11-13 | [Construction and exact proof](verifications/improved/dim61-62-66-affine-reconstruction/) | — |
| 67 | 331 737 996 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291) ⊕ A3 *(derived)* | **355 095 446** | **SOTA** | +7.0409% | — | since 2026-11-13 | [Compatible affine sign classes](verifications/improved/dim61-67-compatible-sign-families/) | — |
| 68 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **472 598 888** | **SOTA** | +42.461494% | — | updated 2026-11-17 | [Exact Gamma72 D4 fiber rescaling](verifications/improved/dim68-91-lattice-code-continuation/) | — |
| 69 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **780 388 888** | **SOTA** | +135.242549% | — | updated 2026-11-17 | [Exact Gamma72 D4 fiber rescaling](verifications/improved/dim68-91-lattice-code-continuation/) | — |
| 70 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **1 361 134 720** | **SOTA** | 4.10x | — | updated 2026-11-17 | [D4-profile compression with anisotropic planar heads and independently replayed Gamma72 populations](verifications/improved/dim70-gamma72-d4-planar-heads/) | — |
| 71 | 331 737 984 | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* | **2 623 578 936** | **SOTA** | +690.858769% | — | updated 2026-11-17 | [Mixed-shell Gamma72 sections](verifications/improved/dim68-91-lattice-code-continuation/) | — |
| 73 | 6 218 175 602 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [1D tail](common/published.py) *(derived)* | **6 218 206 932** | **SOTA** | +0.000504% | — | since 2026-08-30 | [Gamma72 cosets plus caps](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 74 | 6 218 175 606 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [2D tail](common/published.py) *(derived)* | **6 218 211 606** | **SOTA** | +0.000579% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 75 | 6 218 175 612 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [3D tail](common/published.py) *(derived)* | **6 218 246 060** | **SOTA** | +0.00113% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 76 | 6 218 175 624 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [4D tail](common/published.py) *(derived)* | **6 218 315 288** | **SOTA** | +0.00225% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 77 | 6 218 175 640 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [5D tail](common/published.py) *(derived)* | **6 218 490 542** | **SOTA** | +0.005064% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 78 | 6 218 175 672 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [6D tail](common/published.py) *(derived)* | **6 218 688 770** | **SOTA** | +0.008252% | — | since 2026-08-30 | [Gamma72 nearest cosets](verifications/improved/dim73-88-gamma72-cosets/) | — |
| 79 | 6 218 175 726 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [7D tail](common/published.py) *(derived)* | **6 219 067 454** | **SOTA** | +0.014341% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 80 | 6 218 175 840 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ E8 (tabulated direct sum) | **6 219 714 034** | **SOTA** | +0.024737% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 81 | 6 218 175 906 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [9D tail](common/published.py) *(derived)* | **6 220 589 740** | **SOTA** | +0.038819% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 82 | 6 218 176 110 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [10D tail](common/published.py) *(derived)* | **6 221 706 936** | **SOTA** | +0.056782% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 83 | 6 218 176 204 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [11D tail](common/published.py) *(derived)* | **6 223 067 174** | **SOTA** | +0.078656% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 84 | 6 218 176 441 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [12D tail](common/published.py) *(derived)* | **6 224 663 766** | **SOTA** | +0.104328% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 85 | 6 218 176 754 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [13D tail](common/published.py) *(derived)* | **6 226 484 034** | **SOTA** | +0.133597% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 86 | 6 218 177 532 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [14D tail](common/published.py) *(derived)* | **6 228 511 388** | **SOTA** | +0.166188% | — | updated 2026-11-17 | [Gamma72 norm-ten shells with exact joint-moment and marginal duals, local triangle bounds and blocked independence recurrence](verifications/improved/dim79-86-norm10-triangle-layers/) | — |
| 87 | 6 218 178 164 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [15D tail](common/published.py) *(derived)* | **7 266 268 830** | **SOTA** | +16.855269% | — | updated 2026-11-17 | [Compressed ternary fibers, nonlinear phases and pure-head roots](verifications/improved/dim87-89-hexagonal-phase-unions/) | — |
| 88 | 6 218 179 920 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [16D tail](common/published.py) *(derived)* | **12 884 901 888** | **SOTA** | +107.213398% | — | updated 2026-11-17 | [Twelve punctured Z8 fibers](verifications/improved/dim88-89-z8-prefix-anticode/) | — |
| 89 | 6 218 181 330 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [17D tail](common/published.py) *(derived)* | **13 958 643 714** | **SOTA** | +124.481130% | — | updated 2026-11-17 | [Thirteen mixed-latitude Z8 fibers](verifications/improved/dim88-89-z8-prefix-anticode/mixed89/) | — |
| 90 | 6 218 183 254 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [18D tail](common/published.py) *(derived)* | **21 943 640 934** | **SOTA** | +252.894729% | — | updated 2026-11-17 | [Nonlinear hexagonal phase layers](verifications/improved/dim90-92-hexagonal-phase-layers/) | — |
| 91 | 6 218 187 548 | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [19D tail](common/published.py) *(derived)* | **22 359 345 896** | **SOTA** | +259.579793% | — | updated 2026-11-17 | [Direct linear-phase latitude layers](verifications/improved/dim91-92-direct-latitude-layers/) | — |
| 92 | 8 589 934 592 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742) *(shortened)* | **29 797 911 144** | **SOTA** | +246.893342% | — | updated 2026-11-17 | [Four-row partition heads with RS/root augmentation](verifications/improved/dim92-four-row-partition-heads/) | — |
| 93 | 17 179 869 184 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742) *(shortened)* | **42 951 379 608** | **SOTA** | +150.009934% | — | updated 2026-11-17 | [Three alphabets, E7 heads and RS/root layers](verifications/improved/dim93-96-e7-rs-root-augmentation/) | — |
| 94 | 34 359 738 368 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742) *(shortened)* | **73 016 281 804** | **SOTA** | +112.505349% | — | updated 2026-11-17 | [E6 and mixed E7 head exchanges with RS/root layers](verifications/improved/dim94-95-e6-e7-head-exchange/) | — |
| 95 | 68 719 476 736 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742) *(shortened)* | **109 523 504 196** | **SOTA** | +59.377675% | — | updated 2026-11-17 | [E6 and mixed E7 head exchanges with RS/root layers](verifications/improved/dim94-95-e6-e7-head-exchange/) | — |
| 96 | 137 438 953 472 | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742) | **206 160 530 240** | **SOTA** | +50.001528% | — | updated 2026-11-17 | [Three alphabets, E7 heads and RS/root layers](verifications/improved/dim93-96-e7-rs-root-augmentation/) | — |

## Recent updates

Expand a week for the main milestones; follow its history link for the full dated sequence, intermediate values and proof notes. Current bounds are in the table above.

<details>
<summary>16–22 November 2026 (through 19 November)</summary>

A larger radical support packing takes **37D to 544,280** (+10,416 over the original QR result). [37D proof](verifications/improved/dim37-qr-radical-families/). The [536,440-point intermediate](verifications/improved/dim37-qr-shared-stars/) is retained.

New offsets raise **38D to 652,028** (+8,192). A nonlocal Coxeter–Todd phase rebuild raises **36D to 502,821** (+54). [38D proof](verifications/improved/dim38-qr-free-offsets/) · [36D proof](verifications/improved/dim36-ct-nonlocal/).

QR coordinate compression gives **533,864 in 37D** (+27,155) and **643,836 in 38D** (+32,038). [Exact proofs](verifications/improved/dim37-38-qr-merge/).

Rebuilding the Coxeter–Todd cell layer raises **36D to 502,767** (+945) and **37D to 506,709** (+1,071). [Exact proof](verifications/improved/dim36-37-ct-cell-rebuild/).

Triangular Hamming latitudes raise **38D to 611,798**, a gain of 160. [Exact proof](verifications/improved/dim38-triangle-hamming/).

The alternate 26D construction reaches **201,192**, tying Takhanov. Fixed-projection obstructions cover the 26D and 27D comparison witnesses. [26D and 27D research](research/26_27_native_projection_obstructions_20261118/)

The 29D projected-grid search gives **232,157**, 3 above the preceding repository witness and 4 above Takhanov’s pinned release. [14-for-11 even-grid exchange](verifications/improved/dim29-even-grid-exchange/).

Exact exchanges give **223,528 in 28D**, **248,416 in 30D**, and **251,524 in 31D**. [28D proof](verifications/improved/dim28-takhanov-root-exchange/) · [30D proof](verifications/improved/dim30-reflected-certificate-unions/) · [31D proof](verifications/improved/dim30-31-midpoint-exchanges/).

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
