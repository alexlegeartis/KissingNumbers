# Verification of the retained SOTA results

Latest 61–66D: [61/62/66 affine subcodes](verifications/improved/dim61-62-66-affine-reconstruction/) and [63–65 affine subcodes and explicit section](verifications/improved/dim63-65-affine-reconstruction/) extend the retained parent proof. The SOTA runner executes both new proofs and the unchanged parent geometry.

Latest 31D: **251,420**, with [annealed-design proof](verifications/improved/dim31-annealed-design/). Run `python verifications/improved/dim31-annealed-design/verify.py`. The historical checkpoints below remain valid.

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

The repository registers 70 constructions. The README labels 68 as best in its
dated reviewed sources; PackingStar has superseded the repository in 26D and 27D. Every construction has a verification route, but a short run does
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
| 26 | Complete asymmetric integer coordinates; every unordered point pair | NumPy; retained 28–30D certificates remain available |
| 38 | Two independent exact proofs rebuild all cap classes, equatorial lines and axes; check 120 cap deletions and all 6,656 expanded latitude points | Golay lattice congruences and weight enumerator prove the parent minimum |
| 39 | Replay sparse supports and root replacements; check every relative pair of 1,281 dense sign orbits | Bounded integer arithmetic and complete Walsh enumeration |
| 32–35 | Current 32D 352992 and 33D 370408 configurations reconstructed with independent exact proofs; the 33D audit checks all 2129016 exceptional pairs in two arithmetic representations, every star category, exact source/deltas, rank 33, 500 missing antipodes and 16 controls for each of its final and unpublished intermediate witnesses; retained 32D controls, all historical proofs and exact 34–35D head/root replays | Python 3.10+ and NumPy; no search or unverified population input |
| 36–37 | Fresh exact source replay and independent reconstruction of all retained base switches; every added/retained phase and latitude interaction, plus malformed-proof controls | Python 3.9+, NumPy and SymPy; no unverified parent population |
| 41–46 | Exact head ranks and angles; two independent 43D tiny-layer replays; explicit tiny bodies and added medium classes | Complete-family population identities and tables are inherited; the latest large shells have not all been freshly enumerated |
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
