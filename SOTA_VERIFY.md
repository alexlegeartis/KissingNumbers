# Verification of the retained SOTA results

The 28 October P48 refilling certificate improves 49–60D. `verify_october.py` also runs `verify_refilled.py`: full int64 base/isometry checks and an independent exact binary64/full-coordinate replay. `scripts/refill_family.py` reconstructs every family array from the exhibited base.

The 27 October [26–31D residual completion package](verifications/improved/dim26-31-residual-completions/) supersedes the earlier 26–31D packages. Its default replay checks all supplied coordinates independently, including a self-contained shell lemma, and uses exact quadratic-field arithmetic for the final 30D/31D pure tails. Historical packages remain reproducible. No new worldwide record survey is asserted.

The README labels 64 results SOTA against its dated sources. The six retained
26–31D constructions are surpassed by PackingStar's reported October bounds
([28 October review](common/source_review_20261028.json)). Every labelled result
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
| 31 | Independent Leech rebuild, all midpoint/centre/pure-tail pairs and expanded-coordinate comparison | Standard Leech shell kissing property; signed exact radical bounds |
| 40 | Independent lattice minimum enumeration and all auxiliary/core pair types | NumPy and C++17 g++ with __int128; exact bounded integer BLAS |
| 26–27 | Regenerate the Leech shell, check deleted endpoints and midpoint identities, and check every inserted/inserted pair exactly | Standard Leech shell kissing property; midpoint convexity and pure-tail Cauchy–Schwarz bounds |
| 28 | Exact Leech partition, head angles and auxiliary layers | Standard kissing property of the generated Leech shell |
| 29–30 | Regenerate the Leech shell, check all endpoint/midpoint identities and every inserted pair exactly | Standard Leech shell kissing property; exact small integer BLAS with no tolerance; pure-tail Cauchy–Schwarz bound |
| 38 | Rebuild all cap classes and equatorial lines; check all 1,352 joint latitude families and retained/coordinate axes with integers | Golay lattice congruences and weight enumerator prove the parent minimum |
| 39 | Replay sparse supports and root replacements; check every relative pair of 1,281 dense sign orbits | Bounded integer arithmetic and complete Walsh enumeration |
| 32–35 | Two independent exact representations: full dense-code enumeration, signed-support checks, normalized head products and explicit coordinate checks | Python 3.9+ and NumPy; no search or unverified population input |
| 36–37 | Two exact routes: fresh short-shell enumeration, phase refill, partial latitude orbits, 97,968 expanded phase points and 1,908 latitude bodies | Python 3.9+ and NumPy; no unverified parent population |
| 41–46 | Exact head ranks and angles; explicit tiny bodies and added medium classes | Complete-family population identities and tables are inherited; the latest large shells have not all been freshly enumerated |
| 49–60 | All 7,387 base lines, 577 integral isometries, 3,901,694 distinct lines and cap directions | Published parent lattice minimum/shell size; the 12D Coxeter–Todd direction count is a literature input |
| 61–63, 95–96 | Rebuild explicit codes, sparse layers, sections and counts | Published parent code inputs stated by the source proofs |
| 64–65 | Fresh sign-code enumeration, support intersections and parity repairs | Parent code inputs stated in the package |
| 66–68 | Layer geometry and counts in quick mode | `--sota` also enumerates all 68,719,476,736 words of the [73,36,16] code; quick mode reads its recorded distance certificate |
| 69–71 | Exact moment/LP certificates and exhibited section Grams | Published Gamma72 minimum and design theorem |
| 73 | Completed lattice certificate, int64 cap surgery, 6235 norm-ten additions and non-lattice witness (`verify73.py`) | Published Gamma72 parent facts; `verify73.py --rebuild` regenerates the block enumerations |
| 77–79 | Exact theta coefficients, shell moments, mixed-tail ranks/products and rational weighted-priority bound (`verify_norm10.py`, or fast `verify_weighted_norm10.py`) | Published Gamma72 minimum and Venkov all-shell 11-design theorem; existence proof without a full selected head list |
| 74–76, 80–86 | Nearest-subset checks, integral isometries, tail ranks and all retained cross-head inequalities | Deferred in quick mode; about an hour for the full replay. The 21 completed block lists and completion metadata are inputs; the original enumerations are not rerun |
| 87–94 | Complete CPU replay of the compact class certificate, direction partitions, axis geometry and counts | Published Gamma72 minimum and shell size |

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

## Mixed-pair support continuation in 95--96D

The code-layer package now checks all 303,330,765 pairs of its 24,631
physical supports. `python verifications/improved/dim61-68-95-96-code-layers/verify95_96.py`
replays this certificate in isolation, including sign distances and the
95D central section. It also checks every pair independently using packed
integer masks and matches both fresh receipts to the pinned certificate.
The current bounds are 13,290,553,536 in 96D (+294,912 over the preceding
repository value) and 7,092,187,916 in 95D (+180,224). The default package
verifier includes both replays.

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

