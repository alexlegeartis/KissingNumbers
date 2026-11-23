# Verification of the retained SOTA results

23 November 2026 — Replacing the remaining pure poles by a 94-point antipodal code improves 59D to **56,969,157**. [Exact replay](verifications/improved/dim59-p48-pure94/).

23 November 2026 — The complete 126-hyperplane atlas gives **481,013,656 in 68D**. An independent seven-chart Γ72 construction raises the bounds in 77–79D, including **6,226,104,558 in 79D**. [Atlas replay](verifications/improved/dim65-68-complete-hyperplane-atlas/) · [Seven-chart replay](verifications/improved/dim77-79-gamma72-seven-charts/). Two complete trine fibres give **502,929 in 36D**. [36D replay](verifications/improved/dim36-ct-double-trine/).

23 November 2026 — Signed sparse shelves and a 24-cell head on 94 affine hyperplanes raise the verified kissing lower bounds to **448,201,606 in 67D** and **479,699,848 in 68D**, gains of 1.33% and 1.50% over the preceding repository constructions. [Exact replay](verifications/improved/dim63-68-affine-hyperplane-heads/).

23 November 2026 — Support-specific Hamming codes allow signed supports to overlap in 9 or 10 coordinates. The new bounds in 63–67D include **216,972,548** in 63D and **442,334,086** in 67D. [Exact replay](verifications/improved/dim63-67-signed-overlap/).

23 November 2026 — Internal weight-12 tags strengthen the 63D, 64D, and 66D shells. A three-cell BCH polar layer and the new equator bring 65D to **413,985,666** (+196,608). [Exact replay](verifications/improved/dim63-66-internal-tags-bch/).

23 November 2026 — A diagonally weighted Pless-frame cube improves 59D by 1,520 points to **56,969,091**. [Exact replay](verifications/improved/dim59-p48-frame-cube/).

23 November 2026 — Weight-12 bodies with axial tags add 65,536 kissing points in 66D and 24,576 in 67D. [Exact replay](verifications/improved/dim66-67-rm2-axial-weight12/).

23 November 2026 — Two alternating-checksum polar shelves add 262,144 kissing points in 65D, reaching **413,789,058**. [Exact replay](verifications/improved/dim65-rm2-polar-lift/).

23 November 2026 — A common-class 48-frame and 32-color tail give **58,665,998** in 60D (+1,814), with literal signed deletion accounting. [Exact replay](verifications/improved/dim60-p48-frame48-tail/).

23 November 2026 — One binary code supplies both dense and sparse layers, improving all six 62–67D lower bounds again. The 64D bound is **413,526,912** (+10,706,944), and three tagged weight-14 shells raise 67D to **437,513,094** (+15,308,784). [Exact replay](verifications/improved/dim62-67-rm2-selfsparse/). Run `python verifications/improved/dim62-67-rm2-selfsparse/verify.py`.

22 November 2026 — Eighty cubic cosets of RM(2,6) give 335,544,320 dense signs and improve all dimensions 62–67. The 64D lower bound is **402,819,968** (+67,021,024), and the 67D lower bound is **422,204,310**. [Exact replay](verifications/improved/dim62-67-rm2-cubic-cosets/).

22 November 2026 — The signed-owner replay corrects the 60D bound to **58,664,184** (+462 over 58,663,722), with 231 compatible ray pairs and no blockers. [Exact replay and correction](verifications/improved/dim60-p48-n12-blockers/).

22 November 2026 — Global owner migration and a second refill improve 50–60D. The 60D construction reaches 58,663,722, combining exact owner changes with the 16-head sign cube and nested tetrahedral palettes. [Exact replay](verifications/improved/dim50-60-p48-global-owner-refill/).

22 November 2026 — New signed-head exchanges and complementary tail directions improve every dimension from 49D through 60D. [Transport proof](verifications/improved/dim50-59-p48-joint-transport/). [Tail proof and replay](verifications/improved/dim54-60-p48-tail-antipodes/) · [49D proof](verifications/improved/dim49-joint-signed-exchange/).

The [62D partial sign refill](verifications/improved/dim62-partial-nonary-refill/) replays with `python verifications/improved/dim62-partial-nonary-refill/verify.py`. It checks the pinned parent, all sign overlaps, and the removal union.

The [71D single-block Γ72 refill](verifications/improved/dim71-gamma72-single-block-refill/) replays with `python verifications/improved/dim71-gamma72-single-block-refill/verify.py`. It checks the block shell and the exact replacement count.

The [77–84D half-root Γ72 proof](verifications/improved/dim77-84-gamma72-half-root-shells/) replays with `python verifications/improved/dim77-84-gamma72-half-root-shells/verify.py`. It verifies all shorter-length gates and convolves exact vertical parity shells with Leech contact counts.

The [73–82D depth-six Γ72 lattice proof](verifications/improved/dim73-82-gamma72-depth6-lattice/) replays with `python verifications/improved/dim73-82-gamma72-depth6-lattice/verify.py` (NumPy and SymPy). It reconstructs the catalogue gluing and short vertical shells, then uses exact theta coefficients and Leech norm classes to count all contacts.

The [92D variable-radius fiber proof](verifications/improved/dim92-z8-variable-radius-heads/) replays with `python verifications/improved/dim92-z8-variable-radius-heads/verify_flat.py` (Node.js for the sparse ancestor).

The [93D variable-radius fiber proof](verifications/improved/dim93-z8-variable-radius-heads/) replays with `python verifications/improved/dim93-z8-variable-radius-heads/verify_flat.py` (Node.js for the sparse ancestor).

The [Marked caps and nine-support refill](verifications/improved/dim61-62-marked-caps/) replays with `python verifications/improved/dim61-62-marked-caps/verify.py`.

The [Balanced signed nine-support refill](verifications/improved/dim63-balanced-nine-refill/) replays with `python verifications/improved/dim63-balanced-nine-refill/verify.py`.

The [62D coordinate-weighted distance upgrade](verifications/improved/dim62-coordinate-weighted-distance/) replays with `python verifications/improved/dim62-coordinate-weighted-distance/verify.py`. It reads the sixty explicit rank-16 kernels and checks every weighted difference and translated word assignment in exact arithmetic.

The [61–62D omitted-fibre construction](verifications/improved/dim61-62-omitted-fibre-spikes/) replays with `python verifications/improved/dim61-62-omitted-fibre-spikes/verify.py` and checks the exact shortened-code distances and weighted geometry.

The [66–67D mixed-depth transfer](verifications/improved/dim66-67-mixed-depth-transfer/) replays with `python verifications/improved/dim66-67-mixed-depth-transfer/verify.py`, including the Kerdock and two-tail sparse interfaces.

The [65D mixed-depth transfer](verifications/improved/dim65-mixed-depth-transfer/) replays with `python verifications/improved/dim65-mixed-depth-transfer/verify.py`, checking BCH coordinate alignment, discarded dense owners and retained GF(8) interfaces.

The [63–64D mixed-depth clones](verifications/improved/dim63-64-mixed-depth-clones/) replay with `python verifications/improved/dim63-64-mixed-depth-clones/verify.py`. This trades against the same unmodified parent stock; gains are not additive to the old clones.

The [63–64D dense owner clones](verifications/improved/dim63-64-dense-owner-clones/) replay with `python verifications/improved/dim63-64-dense-owner-clones/verify.py` and validate all altered signs and pair products. These are alternative witnesses, not an append to the GF(8) layers.

The [GF(8) norm pole replacement](verifications/improved/dim63-67-gf8-norm-poles/) improves 63–67D. Run `node verifications/improved/dim63-67-gf8-norm-poles/verify.cjs` for exact checks of the exterior cosets, the 63D section, and the new vectors.

The [signed Netto trades](verifications/improved/dim66-67-netto-signed-trades/) refine the 66D and 67D results to **342,389,574** and **355,195,404**, with exact nine-support phase and signed-triple checks.

The [exterior BCH cap bodies](verifications/improved/dim65-67-bch-exterior-caps/) raise 65D to **336,256,796**, 66D to **342,389,558**, and 67D to **355,195,391**. Replay `python verifications/improved/dim65-67-bch-exterior-caps/verify.py`; the optional Node checker establishes the scoped fixed-cap no-append result in 65D.

The [Kerdock-axis owner lifts](verifications/improved/dim66-67-kerdock-axis-owners/) replay with `python verifications/improved/dim66-67-kerdock-axis-owners/verify.py`.

The [Complementary parity six-row packings](verifications/improved/dim94-96-hexary-parity-packing/) replay with `python verifications/improved/dim94-96-hexary-parity-packing/verify.py`.

The [Delsarte–Goethals dense owners](verifications/improved/dim65-dg-dense-owners/) replays with `python verifications/improved/dim65-dg-dense-owners/verify.py`.

The [Signed tail-support Netto exchange](verifications/improved/dim66-tail-signed-nine-refill/) replays with `python verifications/improved/dim66-tail-signed-nine-refill/verify.py`.

The [64–65D partial signed fibers](verifications/improved/dim64-65-partial-nonary-signs/) replay with `python verifications/improved/dim64-65-partial-nonary-signs/verify.py`.

The [67D sign refill](verifications/improved/dim67-signed-nine-refill/) replays with `python verifications/improved/dim67-signed-nine-refill/verify.py`; its independent scalar check is `python verifications/improved/dim67-signed-nine-refill/audit.py`.

The [A2 normal fibers of Gamma72](verifications/improved/dim70-gamma72-a2-fibers/) replay with `python verifications/improved/dim70-gamma72-a2-fibers/verify.py`.

The [alternate Γ72 A2 normal embedding](verifications/improved/dim70-gamma72-alternate-a2-normals/) replays with `OPENBLAS_NUM_THREADS=1 python verifications/improved/dim70-gamma72-alternate-a2-normals/verify.py` (NumPy, SymPy, g++), counting the three pairing classes directly.
The [mixed Γ72 A2 normal embedding](verifications/improved/dim70-gamma72-mixed-a2-normals/) replays with `OPENBLAS_NUM_THREADS=1 python verifications/improved/dim70-gamma72-mixed-a2-normals/verify.py` (NumPy, SymPy, g++ and Node.js), using an exact link enumeration and harmonic identities.

The [Three-marginal Gamma72 proof](verifications/improved/dim77-85-delta-marginal/) replays with `python verifications/improved/dim77-85-delta-marginal/verify.py`.

The [88D full Pless census](verifications/improved/dim88-full-pless-chord-marks/) replays with `python verifications/improved/dim88-full-pless-chord-marks/verify.py` and a C++17 compiler. It certifies the six-column QR projection rank and enumerates the sparse weight-15 words from their unique lighter information half.

The [rooted exchange in 74–85D](verifications/improved/dim74-85-gamma72-rooted-exchange/) replays with `python verifications/improved/dim74-85-gamma72-rooted-exchange/verify.py`. It uses exact double counts at each graph root, rational triangle duals and the inherited integer recurrence. These are existence bounds.

The [finer exact Γ72 recurrence](verifications/improved/dim74-85-gamma72-fine-recurrence/) replays with `python verifications/improved/dim74-85-gamma72-fine-recurrence/verify.py`, verifies the rooted parent, and uses up to 32,768 blocks to certify the higher integer bounds in dimensions 74–85.

The [harmonic-theta 74–85D latitude bounds](verifications/improved/dim74-85-gamma72-harmonic-latitudes/) replay with `python verifications/improved/dim74-85-gamma72-harmonic-latitudes/verify.py`. The verifier rebuilds the harmonic shell identities, exact linear dual bounds and downward-rounded graph recurrence. These are existence bounds.

The [Shaped Gamma72 norm-ten latitudes](verifications/improved/dim76-85-mixed-shape-norm10/) replays with `python verifications/improved/dim76-85-mixed-shape-norm10/verify.py`.

The [Re-marked axial Pless and Steiner squares](verifications/improved/dim86-87-remarked-square/) replays with `python verifications/improved/dim86-87-remarked-square/verify.py`.

The [Oriented polyhedral norm-ten latitude](verifications/improved/dim75-76-oriented-norm10-latitudes/) replays with `python verifications/improved/dim75-76-oriented-norm10-latitudes/verify.py`.

The [latest 89D mixed matching-row construction](verifications/improved/dim89-mixed-matching-rows/) gives **15,184,866,834**; run `python verifications/improved/dim89-mixed-matching-rows/verify.py` for exact mixed pattern and inner-code checks.

The [new 88–89D signed plane and row designs](verifications/improved/dim88-89-signed-row-designs/) are verified by `python verifications/improved/dim88-89-signed-row-designs/verify.py`. This checks the support and signed-row geometry against the pinned parents.

The [94–95D mixed triality append](verifications/improved/dim94-95-mixed-triality/) replays with `python verifications/improved/dim94-95-mixed-triality/verify.py` with the true physical RS grouping.

The [94–95D adapted dodecad sign fibres](verifications/improved/dim94-95-adapted-dodecad-fibres/) replay with `python verifications/improved/dim94-95-adapted-dodecad-fibres/verify.py`, including the latest mixed and partial triality and matching layers in the compatibility checks.

The [94–95D coordinate-aligned dodecad fibres](verifications/improved/dim94-95-coordinate-aligned-dodecad-fibres/) replay with `python verifications/improved/dim94-95-coordinate-aligned-dodecad-fibres/verify.py`, checking the chosen axes, adapted sign kernels, and partial sign-fibre selections against the retained stock.

The [91–93D complementary parity Steiner append](verifications/improved/dim91-93-complementary-parity-steiner-append/) replays with `python verifications/improved/dim91-93-complementary-parity-steiner-append/verify.py`, verifying the mixed frame geometry and retained stock.

The [91–93D non-Steiner triality append](verifications/improved/dim91-93-nonsteiner-triality-append/) replays with `python verifications/improved/dim91-93-nonsteiner-triality-append/verify.py`, checking all affected retained modes and the fourteen new ones.

The [93D bipartite D5 head construction](verifications/improved/dim93-z8-split-d5-heads/) replays with `python verifications/improved/dim93-z8-split-d5-heads/verify.py`. It checks 46 full Z8 fibres, the interfaces to the retained body layers, and the 1,800 new roots.

The [94–95D mixed-tetrad axis and sign fibres](verifications/improved/dim94-95-mixed-tetrad-axis-sign-fibres/) replay with `python verifications/improved/dim94-95-mixed-tetrad-axis-sign-fibres/verify.py` and a C++17 compiler. The replay checks the retained post-refinement monochrome geometry, matching modes, and actual shifted cross-frame interfaces.

The [96D adapted sign kernels](verifications/improved/dim96-adapted-dodecad-fibres/) replay with `python verifications/improved/dim96-adapted-dodecad-fibres/verify.py`. Each kernel is reconstructed from its active projection equations, checked at minimum Hamming distance three, and checked against the unchanged post-RS stock.

The [96D mixed-tetrad axis and sign fibres](verifications/improved/dim96-mixed-tetrad-axis-sign-fibres/) replay with `python verifications/improved/dim96-mixed-tetrad-axis-sign-fibres/verify.py` and a C++17 compiler. The replay checks same-frame and off-frame signs against the retained mixed layers and verifies geometric compatibility with the adapted monochrome stock.

The [96D row-alphabet repack](verifications/improved/dim96-mixed-tetrad-row-alphabet-repack/) replays with `python verifications/improved/dim96-mixed-tetrad-row-alphabet-repack/verify.py` and a C++17 compiler, checking the prior mixed layer and the revised axis alphabets.

The [96D joint dodecad fibres](verifications/improved/dim96-joint-dodecad-fibres/) replay with `python verifications/improved/dim96-joint-dodecad-fibres/verify.py`, checking individual axis words and their Hamming-syndrome sign fibres against all active projections.

The [post-RS partial-fibre packings](verifications/improved/dim94-96-post-rs-partial-fibres/) replay with `python verifications/improved/dim94-96-post-rs-partial-fibres/verify.py`, checking each explicit set of axis words directly against the repacked and untouched interfaces.

The [164-word axis-code repack](verifications/improved/dim89-96-rs-axis-code-repack/) replays with `python verifications/improved/dim89-96-rs-axis-code-repack/verify.py`, regenerating the restricted Reed–Solomon code and checking every replacement eligibility and count in 89D and 91–96D.

The [96D partial triality fibre append](verifications/improved/dim96-partial-triality-fibres/) replays with `python verifications/improved/dim96-partial-triality-fibres/verify.py`. It directly checks actual retained axis words and all old and new interfaces.

The [96D affine-phase append](verifications/improved/dim96-affine-phase-projections/) replays with `python verifications/improved/dim96-affine-phase-projections/verify.py`, verifying projected affine codes across all interfaces with the earlier mixed and hexary patterns.

The [96D mixed triality append](verifications/improved/dim96-mixed-triality-six-rows/) replays with `python verifications/improved/dim96-mixed-triality-six-rows/verify.py` and checks the pinned grouped roots, flat and Reed–Solomon interfaces.

The [94–95D matching-row append](verifications/improved/dim94-95-matching-row-append/) replays with `python verifications/improved/dim94-95-matching-row-append/verify.py` and checks the pinned Golay, flat and Reed–Solomon interfaces.

The [91–93D triality and matching-row transfer](verifications/improved/dim91-93-triality-matching-transfer/) replays with `python verifications/improved/dim91-93-triality-matching-transfer/verify.py` against pinned sparse and axis code inputs.

The [86/87/90D rotated plane trades](verifications/improved/dim86-87-90-rotated-plane-trades/) replay with `python verifications/improved/dim86-87-90-rotated-plane-trades/verify.py`.

The [88D axial Pless construction](verifications/improved/dim88-axial-pless/) replays with `python verifications/improved/dim88-axial-pless/verify.py`.

The [86–87D axial Pless construction](verifications/improved/dim86-87-axial-pless/) replays with `python verifications/improved/dim86-87-axial-pless/verify.py`.

The [86–87D QR48 chord sections, Pless phases and E8 poles proof](verifications/improved/dim86-87-ternary-chord-sections/) replays with `python verifications/improved/dim86-87-ternary-chord-sections/verify.py`.

The [91–93D flat six-/seven-row parity layers proof](verifications/improved/dim91-93-flat-parity-extensions/) replays with `python verifications/improved/dim91-93-flat-parity-extensions/verify.py`.

The [94–96D BCH flat syndrome gluing proof](verifications/improved/dim94-96-bch-flat-phase/) replays with `python verifications/improved/dim94-96-bch-flat-phase/verify.py`.

The [94–96D paired-group BCH append](verifications/improved/dim94-96-bch-paired-group-append/) replays with `python verifications/improved/dim94-96-bch-paired-group-append/verify.py`.

The [91D asymmetric head-plane proof](verifications/improved/dim91-z8-asymmetric-planes/) runs with `python verifications/improved/dim91-z8-asymmetric-planes/verify_flat.py` (Node.js for the inherited sparse layer). Its analytic alternative head assignments in 91–93D are checked by `heads_analytic.py`.

The [expanded 92D variable-radius proof](verifications/improved/dim92-z8-variable-radius-heads-continuation/) replays with `python verifications/improved/dim92-z8-variable-radius-heads-continuation/verify_flat.py` (Node.js for the sparse ancestor).

The [expanded 93D variable-radius proof](verifications/improved/dim93-z8-variable-radius-heads-continuation/) replays with `python verifications/improved/dim93-z8-variable-radius-heads-continuation/verify_flat.py` (Node.js for the sparse ancestor).

The [95D 108-head proof](verifications/improved/dim95-e7-108-heads/) checks all 5,778 head pairs and the unchanged sparse-layer bounds. Run `python verifications/improved/dim95-e7-108-heads/verify.py`.

The [signed nine-support layers](verifications/improved/dim61-66-nonary-affine-caps/) replay with `python verifications/improved/dim61-66-nonary-affine-caps/verify.py`; the original body stocks are inherited unchanged.

The [variable-puncture Golay supports proof](verifications/improved/dim88-89-golay-variable-punctures/) replays with `python verifications/improved/dim88-89-golay-variable-punctures/verify.py`.

The [flat Reed–Solomon layers proof](verifications/improved/dim92-96-flat-rs-layers/) replays with `python verifications/improved/dim92-96-flat-rs-layers/verify_flat.py`.

The [two-color 89D latitude replay](verifications/improved/dim89-two-color-latitudes/) checks all 91 cross-fiber pairs and the unchanged sparse layer. Run `python verifications/improved/dim89-two-color-latitudes/verify.py` (with Node.js installed).

The [enlarged 92–96D outer lists](verifications/improved/dim92-96-shifted-dodecacode-outer-extension/) have two exact replays: `node verifications/improved/dim92-96-shifted-dodecacode-outer-extension/verify.js` and `python verifications/improved/dim92-96-shifted-dodecacode-outer-extension/verify_shifted.py`.

The [shifted-dodecacode replay](verifications/improved/dim92-96-shifted-dodecacode/) checks all three quotient maps, the 637 rescued outer pairs across three lists, and all remaining pairs. Run `node verifications/improved/dim92-96-shifted-dodecacode/verify.js`.

The [90D/91D rotated-square replay](verifications/improved/dim90-91-a2-rotated-square/) checks all 15,180 triples of the 46-point Steiner system and the 3,465 surviving support blocks, 128 even signs per support, symbolic geometry, and root-layer counts. Run `python3 verifications/improved/dim90-91-a2-rotated-square/verify.py`.

The [shortened Golay–dodecacode replay](verifications/improved/dim88-95-golay-shortening/) verifies the 629-pattern 88D/89D construction and larger 92–95D pattern lists, and checks a shortened dodecacode and all distances in the 144 binary signs. Run `node verifications/improved/dim88-95-golay-shortening/verify.js`.

The [92–96D Golay–dodecacode augmentation](verifications/improved/dim92-96-golay-dodecacode/) verifies three outer pattern lists (309, 938 and 2,583 patterns) and their two exact inner codes. Run `node verifications/improved/dim92-96-golay-dodecacode/verify.js`; the proof checks compatibility with the pinned dense, Reed–Solomon and root layers without enumerating billions of points.

Latest 37D: **544,666** (+192 over the radical-family proof). The [12-flat repacking](verifications/improved/dim37-qr-flat-repack/) replaces 384 octad phases by 576 points, checking every new pair and all inherited-layer interfaces exactly.

Preceding 37D: **544,474**, 10,610 above the original QR construction. The [radical-layer verification](verifications/improved/dim37-qr-radical-families/) checks all 57,797,376 new-layer pairs and exact cross-layer bounds.

Intermediate 37D: **536,440** (+2,576), with 720 stars replacing 576 roots and a further 2,432 seven-site points. Independent exact [sparse-layer checks](verifications/improved/dim37-qr-shared-stars/) cover all affected pairs and the signed-root ledger.

Latest 38D: **652,076** (+4). The [octad-switch verification](verifications/improved/dim38-qr-octad-switches/) checks the ten replacements and all final support pairs, and independently regenerates the signed stars and root deletions. The 31 dense fibres are unchanged and replayed.

The latest 36D construction gives **502,902** (+81) by switching individual Coxeter–Todd triangles and refilling the three-block layer. The two replays check the quotient geometry and 792,834,848 explicit pairs involving changed points. [36D triangle switch proof](verifications/improved/dim36-ct-triangle-refill/).

The preceding 36D result of **502,821** (+54) replays all 13,002,450 outer coset pairs. [36D phase proof](verifications/improved/dim36-ct-nonlocal/).

The preceding 38D result of **652,044** (+8,208) replays all 465 inter-fibre comparisons. The signed star replay checks every added star and the 128 root deletions. [38D offset proof](verifications/improved/dim38-qr-free-offsets/) · [38D star trade](verifications/improved/dim38-qr-star-trades/).

Preceding 37D and 38D: **533,864** and **643,836**, gains of 27,155 and 32,038. The [QR compression proof](verifications/improved/dim37-38-qr-merge/) checks the complete dense and sparse layers, all product classes in two formulations, and an optional full dense-pair expansion.

Earlier 36D and 37D: **502,767** and **506,709**, gains of 945 and 1,071. The [Coxeter–Todd cell proof](verifications/improved/dim36-37-ct-cell-rebuild/) checks short-shell enumerations, all 12,992,253 equator coset pairs, and all 1,140,333 relevant latitude coset pairs in two exact replays.

Preceding 38D: **611,798**, a gain of 160. The [triangular Hamming proof](verifications/improved/dim38-triangle-hamming/) gives two exact coordinate reconstructions and changed-pair checks, inheriting the unchanged Leech minimum-norm argument.

The latest 33D witness has **370,528** points, up 112. The complete signed geometry checks the dense code, the full octad layer and every explicit point. [33D refill proof](verifications/improved/dim33-joint-signed-refill/).

The preceding 33D witness has **370,416** (+8). Two [exact replays](verifications/improved/dim33-octad-switches/) check the octad switches, all exceptional interactions and the signed-root ledger.

Latest 34D: **391,004** (+12). The [octad-switch proof](verifications/improved/dim34-octad-switches/) checks the support and root ledgers, integer angles and all 391,004 distinct coordinates of exact rank 34.

Latest 35D: **415,360** (+2). The [exact replay](verifications/improved/dim35-joint-star-root/) checks all five layers and the chosen star/root replacement in exact integer arithmetic.

Preceding 35D: **415,358** (+22). The [octad-switch proof](verifications/improved/dim35-octad-switches/) has two independent integer replays of the changed supports, all sparse layer interactions and the dense code.

Preceding 35D: **415,336**, two above its predecessor. The [shared-root proof](verifications/improved/dim35-shared-root-heads/) checks the complete code and support layers and a separate expansion verifies all 67,258,188 changed/anchor pairs.

Latest 28D: **223,528**, 32 above Takhanov’s pinned construction. The [E8 root exchange](verifications/improved/dim28-takhanov-root-exchange/) checks all **24,982,271,628** unordered final pairs and separately verifies the 320-for-288 exchange in signed int64.

The latest 29D witness has **232,157** points, 3 above the previous certificate and 4 above Takhanov’s pinned release. The direct integer replay checks all 3,250,156 donor–incumbent products and 91 donor pairs; `--full` checks all **26,948,320,246** final unordered pairs. [14-for-11 even-grid exchange](verifications/improved/dim29-even-grid-exchange/).

The latest 30D witness has **248,506** points. The verifier checks all **30,877,491,765** unordered pairs and independently checks every pair involving one of the 340 insertions. [30D midpoint packing](verifications/improved/dim30-row-midpoint-packing/)

The preceding 30D witness has **248,416** points, a gain of 65. The verifier checks all **30,855,130,320** unordered pairs and rebuilds the exact 456-for-391 exchange. [30D reflected unions](verifications/improved/dim30-reflected-certificate-unions/)

The preceding 30D witness has **248,351** points, a gain of 62 over 248,289. The verifier checks all **30,838,985,425** unordered pairs and rebuilds the 188-for-126 exchange. [30D contact-sum proof](verifications/improved/dim30-contact-sum-exchanges/).

The joint midpoint package verifies **251,524** in 31D and the preceding 30D witness of **248,289**. [Joint midpoint proofs](verifications/improved/dim30-31-midpoint-exchanges/).

Latest 39D: **783,468**. The [three-coordinate QR certificate](verifications/improved/dim39-qr-three-coordinate-heads/) verifies all exact quotient distances, projected sparse sign gates, and the full root/star layer.

Earlier 39D: **774,972** (up 5,120). The [five-head QR verification](verifications/improved/dim39-qr-pentapolar/) checks the shortened kernel, all fibre distance minima, all sparse supports, and exact sparse–dense angular inequalities.

The 44D verifier checks a 48D Gram matrix, an extremal lattice minimum and a complete 1,656-point section. Spherical moments and exact class distances certify all 4,940,961 projected points. [Witness and replay](verifications/improved/dim44-d4-section/).

The 41D–43D affine-shell verifier reconstructs every family from its Gram matrix and residue. Short-difference enumerations and direct integer products certify the angles. [Witnesses and full replay](verifications/improved/dim41-43-affine-sections/).

The latest 40D result gives **1,126,365** with a joint short-layer exchange in the \(P_{48q}\) complement. The replay proves the lattice minimum, enumerates the 1,084,590-point affine core and checks all remaining unordered pairs. [40D P48q proof](verifications/improved/dim40-p48q-affine-refill/).

The \(P_{48n}\) construction gave **1,126,327** with an alternate affine centre and thirteen individual short-shell points. Its two replays check the gluing geometry and enumerate the affine core. [40D P48n proof](verifications/improved/dim40-affine-centre-refill/).

The preceding 40D result gave **1,126,322** with eight minimum-shell points. Exact verification checks all 9,010,512 new–old products and all 28 new–new products. [40D affine shell proof](verifications/improved/dim40-affine-shells/).

Latest 26D: **201,204** (12 above Takhanov). The [section rebuild](verifications/improved/dim26-section-rational-kernel/) restores the Leech roots, verifies the new 340-point kernel in exact rational arithmetic, and offers a complete exact exterior-pair check with `--full`.

Retained 26D and 27D certificates give **201,188** and **209,884**, respectively **8 and 22 above** pinned PackingStar. Takhanov’s later release gives **201,192** and **210,108**. The [radical-orbit and E8 head-frame proofs](verifications/improved/dim26-27-radical-orbits/) reconstruct the 160-for-152 and 238-for-216 exchanges and check all **20,238,205,078** and **22,025,541,786** final unordered pairs exactly. Run `python verifications/improved/dim26-27-radical-orbits/verify.py` and `python verifications/improved/dim26-27-radical-orbits/negative_controls.py`. [Takhanov’s release](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1)

The alternate **201,192**-point 26D witness passes complete integer and rational verification. Fixed-projection verifiers cover the two Takhanov witnesses and the 26D alternate. [26D and 27D research](research/26_27_native_projection_obstructions_20261118/)

Retained 30D: **248,212**, +32 over the retained native design and +220 over pinned PackingStar. The [exact certificate-union proof](verifications/improved/dim30-certificate-unions/) reconstructs 52 additions and 20 removals and checks all **30,804,474,366** final unordered pairs. Run `python verifications/improved/dim30-certificate-unions/verify.py` and `python verifications/improved/dim30-certificate-unions/negative_controls.py`. The [native-design package](verifications/improved/dim29-30-native-partial-design/) retains the 231,814-point 29D repository certificate and the 248,180-point 30D certificate.

Latest 61–66D: [61/62/66 affine subcodes](verifications/improved/dim61-62-66-affine-reconstruction/) and [63–65 affine subcodes and explicit section](verifications/improved/dim63-65-affine-reconstruction/) extend the retained parent proof. The SOTA runner executes both new proofs and the unchanged parent geometry.

Retained 31D: **251,420**, with [annealed-design proof](verifications/improved/dim31-annealed-design/). Run `python verifications/improved/dim31-annealed-design/verify.py`. The historical checkpoints below remain valid.

13 November 2026 continuation: [31D joint design](verifications/improved/dim31-joint-design/) gives **251,396**, +1 over the retained 251,395 certificate and +142 over the pinned PackingStar release. Full exact final geometry covers **31,599,848,710** unordered pairs, and direct int64 exchange checks verify all added-point interactions. Run `python verifications/improved/dim31-joint-design/verify.py`.

13 November 2026: the [30D Leech–Sloane section](verifications/improved/dim30-leech-sloane-section/) and [31D Leech–Sloane module](verifications/improved/dim31-leech-sloane-module/) give **248,084** and **251,395**, 92 and 141 above PackingStar's pinned November 11 release. Complete exact final geometry covers **30,772,711,486** and **31,599,597,315** unordered pairs. Separate integer checkers rebuild the full 2,564-point recipe, verify its section and exact orientations, reconstruct the pinned sources and check every inserted-point interaction. Four malformed-input controls and complete manifests are replayed. The earlier [248,020-point Barnes–Wall package](verifications/improved/dim30-bw14-module/) remains available. The [fresh comparison](common/source_review_20261113_26_31.json) records the remaining 26D/27D gaps.

The 13 November continuation has three new verification routes. The [Coxeter–Todd orbit exchange](verifications/improved/dim36-37-orbit-exchange/) gives **501,822 in 36D** and **505,638 in 37D**. Its fresh replay reconstructs the source proof and checks every interaction involving the added phase orbits; malformed witnesses, including a correctly counted but geometrically invalid restored orbit, are rejected.

The [extended Gamma72 moments](verifications/improved/dim68-69-gamma72-extended-moments/) give **471,886,208 in 68D** and **746,325,051 in 69D**. Two independent standard-library proofs enumerate the complete 162,969-profile D7 domain, rebuild the nested D4 orientation objective and exact moment rows, check rational duals and all geometric pair types, and validate the two additional 69D poles. Thirteen production rejection controls pass. The published Gamma72 minimum, shell population and spherical 11-design remain explicit inputs.

The [compatible-sign-family package](verifications/improved/dim61-67-compatible-sign-families/) supplies the new 61–67D bounds. The default replay freshly checks finite code/support/sign conditions and the layer formulas. Its independent 67D checker constructs the 49 added affine classes by filtering all 65,536 sign words, validates 357 reordered parent classes, and checks 115,998,720 exact sign-pair products; the new count is **355,095,446**. In 66D and 67D, direct checks of every retained parent-support pair replace reliance on the full 73-bit parent-code minimum. The historical billion-word replay remains available in its retained package; the new construction is checked through its actual retained supports.

The 13 November [signed P48 refill](verifications/improved/dim49-60-p48-signed-refill/) improves all 49–60D bounds, reaching **58,631,571 in 60D**. The class model permits individual signed vectors without an antipodality requirement. Verification checks complete finite classes, all changed and retained interactions, full-coordinate ownership and tail geometry; the reusable new class also has an independent scalar-integer all-pairs replay. Published P48 minimum and shell population remain inputs.

The 13 November [69D slice-compression proof](verifications/improved/dim69-gamma72-head-compression/) establishes **743,522,304** points, **54,069,382 (+7.8424%)** above the preceding certificate. Two standard-library exact checkers independently reconstruct all 1,681 profiles, 16 averaging orbits, 10 moment rows and 55 geometric pair types. Published Gamma72 minimum, shell population and design facts are explicit premises; the certificate proves a favorable orientation exists.

The 12 November [28D/29D sextet continuation](verifications/improved/dim28-29-november-sextets/) gives **222,120** and **231,776**, respectively **16** and **2** above PackingStar’s November 11 files. Complete exact integer geometry covers every unordered pair; direct-int64 changed-pair checks and four rejection controls are included. At that date 26D, 27D, 30D and 31D remained below that release; the November 13 campaign above subsequently reclaimed 30D and 31D.

The 8 November [31D latitude-contact exchange](verifications/improved/dim31-latitude-contact-exchange/) gives **251,200** points, four above the preceding repository certificate. Both complete exact pair proofs cover all 31,550,594,400 unordered pairs. One uses direct integer arithmetic without antipodal compression; reconstruction, all changed-pair checks and semantic rejection controls are included. This proof assumes no lattice minimum or unenumerated population. The source-comparison scope is recorded separately.

The 8 November [pointwise P48 refill](verifications/improved/dim49-60-p48-pointwise-refill/) gives twelve improved 49–60D bounds, including **58,627,711 in 60D**. Its portable replay uses two exact head implementations, every added/retained and added/added pair, full-coordinate ownership, complete direction/pole checks and an algebraic lift proof. Twelve rejection controls test malformed inputs. Published P48 minimum six and shell population 52,416,000 remain inherited; the full lifted code is not expanded.

The 6 November [840-direction lift](verifications/improved/dim60-antipodal840-lift/) gives **58,543,525 in 60D**. Its replay checks all direction/pole/base-class pairs, 420 isometries and all selected mother-line disjointness. The 52,416,000-point P48 shell population and minimum six remain published inputs; the full lifted code is not expanded.

The 29 October [joint-refill certificates](verifications/improved/dim27-31-joint-refills/) give 209,807, 221,988, 231,697, 247,972, 251,188 in 27–31D, respectively, exceeding PackingStar by 47, 124, 21, 5, 19. Complete exact all-pairs checks and direct-int64 inserted-point checks pass. Binary32 and binary64 replays, pinned-source rebuilds and six rejection controls are supplied. Corrected sextet symmetries reclaim 28D and strengthen three further dimensions.

The retained synchronized 28 October 26–31D results use the new
[asymmetric exchanges](verifications/improved/dim26-30-asymmetric-exchanges/),
[27D quadratic module](verifications/improved/dim27-quadratic-normal-module/)
and [31D normal-module exchange](verifications/improved/dim31-normal-module-exchange/).
Their full coordinate checks establish every norm, distinctness condition and
point-pair inequality, including unpaired points. These replace the current
routes for the retained residual-completion certificates.

The further 28 October 33–37D continuation adds shared root-star deletions in
33–35D and five phase orbits in 36–37D, preserving every latitude body. Both
proof representations and deliberate rejection controls pass. These gains
are 8, 8, 20, 15 and 15 points over the preceding repository certificates.

The 28 October 38D continuation gives **611,522** points (+1,128). The
[small constructive recipe](verifications/improved/dim38-joint-latitudes/refill/README.md)
replaces the latitude layer and deletes 120 caps. Both exact proof
representations replay all retained layers; the independent coordinate
representation checks all 22,147,840 distinct latitude pairs. Reconstruction
is deterministic and both verifiers reject four certificate corruptions.

The 28 October 43D continuation reaches 2,919,782 (+39) by continuous head placement. Two independent exact representations check the complete tiny layer; the complete-family population inputs remain documented in the [43D continuation](verifications/improved/dim41-46-p48-heads/continuation43/README.md).

The further 28 October continuation strengthens 49–63D and 66–68D. The P48 replay now checks the 7,783-line base. `dim61-68-95-96-code-layers/verify_continuation.py` replays the weight-eight bridge, parity-compatible sections and complete quadruple designs, and the package default replay includes all of them. The 64D parent was also freshly enumerated and checked.

The 28 October P48 refilling certificate improves 49–60D. `verify_october.py` also runs `verify_refilled.py`: full int64 base/isometry checks and an independent exact binary64/full-coordinate replay. `scripts/refill_family.py` reconstructs every family array from the exhibited base.

The 27 October [26–31D residual completion package](verifications/improved/dim26-31-residual-completions/) supersedes the earlier 26–31D packages. Its default replay checks all supplied coordinates independently, including a self-contained shell lemma, and uses exact quadratic-field arithmetic for the final 30D/31D pure tails. Historical packages remain reproducible. No new worldwide record survey is asserted.

The repository registers 70 constructions. The README labels all 70 as best in its
dated reviewed sources; the November 14 certificates exceed pinned PackingStar in 26D and 27D. Every construction has a verification route, but a short run does
not freshly establish all of them.
The commands below check constructions; they do not perform a literature search.

From this repository's root:

```bash
python run_all.py --sota --quick          # short construction replays and registry audit
python run_all.py --sota --quick --changed # reuse passing checks only while their inputs match
python run_all.py --sota --list           # show the longer replays and coverage before running
python run_all.py --sota                  # also replay cosets, code enumeration and cap geometry
```

Before the 29–30D, 38–39D and norm-ten additions, the quick profile completed in about 3.5 minutes;
the unchanged-input `--changed` repeat completed in under 3 seconds. These are
local timings, and the long coset/code work is intentionally outside that short
profile.

The ordinary `run_all.py` modes remain available for historical constructions,
negative results and manuscript checks. `--sota` selects the current construction
packages; it always uses the audit's `--repository-only` mode. `--quick` selects
the short replays without weakening their angle checks. The per-dimension report
states which checks were deferred. `--changed` means reuse of a local passing run,
not a new mathematical verification. All commands require assertions enabled.

| Dimensions | Available verification | Remaining inputs or longer work |
|---|---|---|
| 18 | Reconstruct all 8,358 points and check every pair exactly | None beyond the stated mathematical construction |
| 25 | Independent Leech rebuild, all 971 circles and 102 moved points; exact biquadratic comparisons | Standard Leech shell kissing property |
| 27–31 | Complete integer coordinates, all unordered point pairs and direct-int64 inserted-point checks; current 30D/31D modules also have independent exact reconstruction | NumPy; no baseline validity assumption. The retained earlier 31D package additionally provides an integer-only C++ full-pair replay |
| 40 | Independent lattice minimum enumeration and all auxiliary/core pair types | NumPy and C++17 g++ with __int128; exact bounded integer BLAS |
| 26 | 27D section, exact rational kernel, and every unordered exterior pair with `--full` | NumPy |
| 38 | Two independent exact proofs rebuild all cap classes, equatorial lines and axes; check 120 cap deletions and all 6,656 expanded latitude points | Golay lattice congruences and weight enumerator prove the parent minimum |
| 39 | Replay sparse supports and root replacements; check every relative pair of 1,281 dense sign orbits | Bounded integer arithmetic and complete Walsh enumeration |
| 32–35 | Current 32D 352992 and 33D 370408 configurations reconstructed with independent exact proofs; the 33D audit checks all 2129016 exceptional pairs in two arithmetic representations, every star category, exact source/deltas, rank 33, 500 missing antipodes and 16 controls for each of its final and unpublished intermediate witnesses; retained 32D controls, all historical proofs and exact 34–35D head/root replays, including the 415336-point 35D construction | Python 3.10+ and NumPy; no search or unverified population input |
| 36–37 | Fresh exact source replay and independent reconstruction of all retained base switches; every added/retained phase and latitude interaction, plus malformed-proof controls | Python 3.9+, NumPy and SymPy; no unverified parent population |
| 41–43 | Fresh counts, head ranks and integer angle checks | Gram matrices and explicit head coordinates supplied |
| 44 | Exact parent minimum, class distances and shell counts by spherical moments | Gram matrix, normals and a complete pair section supplied |
| 45–46 | Exact head ranks and angles, including the added 45D medium classes | Complete-family population identities and tables are inherited; the latest large shells have not all been freshly enumerated |
| 49–60 | All 7,387 base lines, 577 integral isometries, 3,901,694 distinct lines and cap directions | Published parent lattice minimum/shell size; the 12D Coxeter–Todd direction count is a literature input |
| 61–67 | Fresh exact compatible-support, sign-family and code replays; direct retained parent-support pairs, exact section/compression geometry and controls | NumPy and C++17 g++; the unchanged auxiliary layer constructions are replayed locally |
| 68–69 | Two independent exact D7 population and D4 geometry proofs; complete 162,969-profile domain, rational duals, integer ceiling and two 69D poles | Published Gamma72 minimum, shell population and spherical 11-design; no ambient shell enumeration |
| 70–71 | Exact moment/LP certificates and exhibited section Grams | Published Gamma72 minimum and design theorem |
| 73 | Completed lattice certificate, int64 cap surgery, 7663 norm-ten heads in two exact arithmetics and non-lattice witness (`verify73.py`) | Published Gamma72 parent facts; `verify73.py --rebuild` regenerates the block enumerations |
| 77–79, 82–85, 87 | Exact theta coefficients, shell moments, mixed-tail ranks/products and rational weighted-priority bound (`verify_norm10.py`; `verify_weighted_norm10.py` covers only 77–79D) | Published Gamma72 minimum and Venkov all-shell 11-design theorem; existence proof without a full selected head list |
| 74–76, 80–81, 86 | Nearest-subset checks, integral isometries, tail ranks and all retained cross-head inequalities | Deferred in quick mode; about an hour for the full replay. The 21 completed block lists and completion metadata are inputs; the original enumerations are not rerun |
| 88–94 | Complete CPU replay of the compact class certificate, direction partitions, axis geometry and counts | Published Gamma72 minimum and shell size |

The SOTA profile is validated against both README.md and RESULTS.md. A new SOTA
dimension or changed package without a matching verification route causes an
error. PASS means the named replay succeeded within its documented scope;
DEFERRED and GAP do not mean verified. Failure or reuse of an audit or shared
cap-algebra check also affects the dependent coverage rows.

## Fast certificates and stronger bounds in 87–94D

On 25 October 2026, the local CPU continuation increased the seed from
2,454 to 2,465 lines and selected its images by maximum marginal coverage.
The shipped 24.4 MB certificate contains that seed, 16,632 explicit integral
isometries and survivor masks, specifying 34,857,986 distinct retained lines.
It improves the registered cap bounds throughout 87–94D. From the cap package:

```bash
python scripts/verify_symmetry.py
python scripts/verify_symmetry.py --complete
```

The first command passed in 6.64 seconds on local CPU. It checks every seed
angle once, then proves every matrix preserves the Gram. Exact factored hashes
and a global sort establish disjointness without materialising every image.
Equal or antipodal vectors always have equal keys, so hash collisions cause
rejection rather than a false PASS. Integer partial-sum bounds justify every
binary32 and binary64 multiplication. The full seed check uses int64.

The second command checks every auxiliary inequality as well. Large pole
products use a rigorous coarse bound, followed by exact Python-integer checks
for every undecided pair. Exact span checks use a small basis and its
annihilator. These shortcuts preserve the original geometric conditions.
The quick SOTA profile runs this complete command and no longer reports a
GPU-family GAP for these eight dimensions.

The older 32,000-class greedy family remains construction history. Its
`scripts/gpu/check_cert.py` checks a run receipt and statistics, not every
angle. Rebuilding that entire discovery run is unnecessary for the stronger
registered bounds: their complete coordinates are specified by the new
certificate, and the checker uses only NumPy, SciPy/SymPy for the auxiliary
checks, and the shipped files. The published Gamma72 minimum and shell size
remain parent inputs.

## Where the long coset replay spends time

The installed 88D report records 11,640 independent product families,
1,059,240 coset pairs and 105,000,212,525 representative products. Its completed
replay took about 34 minutes on the source machine. The package's full replay
also checks the 4-seed and 8-seed constructions and their common nearestness
inputs. These are coordinate products, not a reconstruction of the discovery
search. The verifier already uses isometries to reuse product families and
proves an absolute-term bound before using exact binary32 arithmetic.

Checking the displayed counts in `final.py` is fast, but does not replace these
geometric checks. For repeated local review, `--changed` avoids repeating a
successful replay when the measured inputs and environment are unchanged.

## Retained mixed-pair support continuation in 95--96D

The code-layer package checks all 338,559,231 pairs of its 26,022 physical
supports. `python verifications/improved/dim61-68-95-96-code-layers/verify95_96.py`
replays the coordinate and packed-integer certificates, sign distances,
95D section counts and GF(4) dual trace obstructions. The latter certify
2,073 supports that occupy thirteen or fourteen rows. Both fresh receipts
must match the pinned files. The retained predecessor bounds are 13,313,343,680 in 96D
(+22,560,768) and 7,101,541,132 in 95D (+9,209,856) over the preceding
repository certificate. The default package verifier includes both replays.

## Norm-ten continuation in 73D and 77--79D

`python verifications/improved/dim73-88-gamma72-cosets/verify_norm10.py`
checks the explicit 73D addition and the exact whole-shell existence bounds
in 77--79D. The latter derive the theta coefficients and design moments,
check all finite tail products, and evaluate the weighted-priority refinement
of Caro--Wei over rational numbers. Mixed shells give 760 tails in 77D and
2232 in 78D; 79D reweights its existing 4158 tails. The fast command
`python verifications/improved/dim73-88-gamma72-cosets/verify_weighted_norm10.py`
replays the new certificates alone, including all 22,848,388 ordered products.
They do not enumerate the selected full-shell head sets. The 73D parent
certificates are replayed as part of this command.


## Further 28 October 73–96D certificates

`verify_norm10.py` in the Gamma72 coset package now checks the 7,663-head
explicit 73D layer in int64 and provably exact binary32 (120,048,558
products per representation), and the signed-orbit existence proofs in
82–85D and 87D. `signed_norm10_certificate.py` derives complete orbit
profiles combinatorially; `verify_signed_profiles.py` re-derives them from
polynomial coefficients and checks all pure-root products and ranks.
The latter proofs retain the published Gamma72/all-shell design inputs
and do not enumerate selected norm-ten heads.

`verify95_96.py` checks the new 24,645-support family: all 303,675,690
support pairs in physical coordinates and again in packed integer form,
plus the exact sign and section counts.

## Nonlinear Golay parent and sections in 88–96D

`python verifications/improved/dim88-96-z8-golay/verify.py` checks both exact
compressed-construction proofs and malformed-certificate controls. Each proof
checks all 16,777,216 Z4 words, 4,166,151 octad-coset exclusions and 338,559,231
support pairs. It derives every dense projection and sparse section count.
The (96,2^37,24) parent is published work of Duursma–Greferath–Litsyn–Schmidt
(2001), previously missed by the linear-code comparison. The 92–96D floors
now use that code and its shortenings; the standard full-chain evaluation is
recorded separately as a derived comparison. No 2^37-word enumeration is claimed.
Both the SOTA profile and CI run this complete package replay.

**17 November 2026: twelve punctured Z8 fibers.** **12,884,901,888 in 88D**, +3,294,110,112 (+34.34659%) over the preceding repository bound. The same bound improves 89D by zero-coordinate embedding. [Complete exact parent, population, pair-type and injectivity proof](verifications/improved/dim88-89-z8-prefix-anticode/); no exhaustive worldwide-priority claim.

Run `python verifications/improved/dim88-89-z8-prefix-anticode/verify.py`.

**17 November 2026: mixed-latitude Z8 fibers in 89D.** **13,958,643,714**, +1,073,741,826 (about +8.33333%) over the twelve-fiber embedding. Thirteen fibers with two opposite nonzero heights plus two pure poles, exact pointwise normalization and full inherited parent replay. [Independent complete proof](verifications/improved/dim88-89-z8-prefix-anticode/mixed89/AUDIT.md); no exhaustive worldwide-priority claim.

