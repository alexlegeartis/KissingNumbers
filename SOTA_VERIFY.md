# Verification of the retained SOTA results

The 29 October [joint-refill certificates](verifications/improved/dim27-31-joint-refills/) give 209,775, 231,684, 247,972 and 251,180 in 27D and 29–31D. Complete exact all-pairs checks and direct-int64 inserted-point checks pass; both binary32 and binary64 replays are exact under their explicit coefficient bounds. Rebuilds and six rejection controls are supplied. The 28D result remains below PackingStar’s 221,864.

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

The README labels 69 results SOTA against its dated sources. Every labelled result
has a verification route, but a short run does not freshly establish all of them.
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
| 27, 29–31 | Complete integer coordinates, all unordered point pairs and direct-int64 inserted-point checks | NumPy; exact bounded products, no baseline validity assumption |
| 40 | Independent lattice minimum enumeration and all auxiliary/core pair types | NumPy and C++17 g++ with __int128; exact bounded integer BLAS |
| 26 | Complete asymmetric integer coordinates; every unordered point pair | NumPy; retained 28–30D certificates remain available |
| 38 | Two independent exact proofs rebuild all cap classes, equatorial lines and axes; check 120 cap deletions and all 6,656 expanded latitude points | Golay lattice congruences and weight enumerator prove the parent minimum |
| 39 | Replay sparse supports and root replacements; check every relative pair of 1,281 dense sign orbits | Bounded integer arithmetic and complete Walsh enumeration |
| 32–35 | Two independent exact representations: full dense-code enumeration, signed-support checks, coordinated heads, shared root deletions and every added-star coordinate product | Python 3.9+ and NumPy; no search or unverified population input |
| 36–37 | Two exact routes: fresh short-shell enumeration, phase refill, partial latitude orbits, 97,983 expanded phase points and 1,908 latitude bodies | Python 3.9+ and NumPy; no unverified parent population |
| 41–46 | Exact head ranks and angles; two independent 43D tiny-layer replays; explicit tiny bodies and added medium classes | Complete-family population identities and tables are inherited; the latest large shells have not all been freshly enumerated |
| 49–60 | All 7,387 base lines, 577 integral isometries, 3,901,694 distinct lines and cap directions | Published parent lattice minimum/shell size; the 12D Coxeter–Todd direction count is a literature input |
| 61–63, 95–96 | Rebuild explicit codes, sparse layers, sections and counts | Published parent code inputs stated by the source proofs |
| 64–65 | Fresh sign-code enumeration, support intersections and parity repairs | Parent code inputs stated in the package |
| 66–68 | Layer geometry and counts in quick mode | `--sota` also enumerates all 68,719,476,736 words of the [73,36,16] code; quick mode reads its recorded distance certificate |
| 69–71 | Exact moment/LP certificates and exhibited section Grams | Published Gamma72 minimum and design theorem |
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
