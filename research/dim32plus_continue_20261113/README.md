# Continuation in dimensions 32–48 and 61 upward

**13 November 2026.** This campaign continues the published checkpoint
`d595c660a8672e0ee436935b74790dadb5d9a838`. It produces eleven further verified
lower-bound improvements. The geometric proofs and count drivers live in the
verification packages; this directory preserves discovery and review history.

| Dimension | Starting bound | New bound | Gain |
|---:|---:|---:|---:|
| 36 | 501,813 | **501,822** | 9 |
| 37 | 505,629 | **505,638** | 9 |
| 61 | 88,260,394 | **88,975,146** | 714,752 |
| 62 | 121,824,268 | **122,539,020** | 714,752 |
| 63 | 184,286,468 | **184,393,988** | 107,520 |
| 64 | 335,496,064 | **335,711,104** | 215,040 |
| 65 | 335,564,290 | **335,779,330** | 215,040 |
| 66 | 341,994,884 | **342,271,364** | 276,480 |
| 67 | 355,090,822 | **355,095,446** | 4,624 |
| 68 | 471,884,288 | **471,886,208** | 1,920 |
| 69 | 743,522,304 | **746,325,051** | 2,802,747 |

These are improvements over the named preceding certificates. The dated
bounded literature check is in `source_review.json`; it does not establish
comprehensive worldwide priority. None of the searches proves the unrestricted
kissing number optimal.

## Three proof packages

1. [36D/37D orbit exchange](../../verifications/improved/dim36-37-orbit-exchange/)
   replaces three phase orbits by six, adding nine points. The inherited
   short shells are freshly enumerated; all 9,067,122 affected products are
   checked exactly, together with the baseline construction. All existing
   latitude orbits remain valid.
2. [61D–67D compatible sign families](../../verifications/improved/dim61-67-compatible-sign-families/)
   optimizes after shortening, permits parity repairs between newly added
   blocks, and adds smaller affine sign classes in 67D. Three complete
   full-class certificates and the 49-class affine extension are explicit.
   Scalar and vector replays check the support geometry, parent words,
   sign kernels, balanced section, owner lift and extra-tail interactions.
3. [68D/69D extended moments](../../verifications/improved/dim68-69-gamma72-extended-moments/)
   embeds the actual D4 normals into an actual D7 normal system. Additional
   mixed moments sharpen the orientation-averaged population bound while
   retaining the existing slice geometry. Two poles complete the 69D count.
   Independent exact checkers enumerate all 162,969 profiles and verify all
   18 moment rows, rational duals and geometric inequalities. The published
   Gamma72 minimum, population and design properties remain explicit inputs.

## What the successful searches changed

**Optimize the retained problem.** The 61D/62D construction previously inherited
the shortening of a family optimized in 64D. Searching after shortening removes
unnecessary conflicts with absent supports. The count of added supports rises
from 588 to 937, each carrying 2,048 points.

**Repair new/new conflicts.** A ten-coordinate overlap can be valid when the
two sign codes impose opposite parities on a common eight-coordinate check.
The original old/new checks leave one Hamming-dual completion free. Using that
freedom for new/new pairs yields 1,940 added supports in 64D: exactly 32,768
full support classes, or 67,108,864 sparse points. The final global certificate
contains 303 such new/new repairs.

**Allow smaller sign classes.** The 67D full-class search did not improve its
1,792-block seed. A separate scan allows proper affine subcodes whose minimum
distance is at least four. A joint 49-class witness contributes 4,624 points
while 357 affected parent classes retain their full 2,048 signs. The independent
review enumerates every sign word and checks 115,998,720 exact sign-pair products
on all critical overlaps. It finds maximum 32 at squared norm 64.

**Use larger actual normal systems.** The D7 moment extension strengthens the
69D ten-slice mean to `6,716,925,440/9`. Its ceiling, plus two compatible poles,
is 746,325,051. This improves the first campaign's 743,522,304 without changing
the ten slice heights. The larger normal system constrains possible shell
populations; it is not assumed to be an ambient lattice symmetry.

## Lanes and finite negative outcomes

| Directory | Work and scope |
|---|---|
| [`codes/`](codes/) | Full-class searches, affine-plane candidate generation, exact export and scalar/vector proof prototypes; 61D–66D improvements. |
| [`gamma/`](gamma/) | D5/D6/D7 normal extensions, exact moment certificates, D7×A1 boundary check and alternate rank-two searches. |
| [`lowdim/`](lowdim/) | 36D/37D orbit exchange, bounded further orbit/pointwise searches, independent D7 and 67D reviews. |
| [`middim/`](middim/) | P48 head/split searches in 40D–48D, independent code review and the successful 67D affine-subcode extension. |

The P48 searches in 40D–48D produced no further record. The 43D body-dependent
scan checked 3,000 bodies and found only two alternative splits of the same
source body, both requiring the same two deletions; they cannot combine for
positive gain. These are finite search outcomes.

The final 37D latitude search rebuilt 168,344 candidates against the changed
equator and found no immediate addition. Its bounded 891-orbit low-blocker
search retained the existing 636 classes. Several higher-shell body families
were excluded by exact base products. The READMEs state which exclusions are
exhaustive for their specified families and which searches merely stopped at
their time limits.

For Gamma72, D7×A1 introduces 296 profile orbits and 43 moment rows but gives
the same tested bounds after exact boundary-membership checks. Alternate
rank-two Grams and the tested affine/head models in 70D/71D do not beat the
existing certificates. A larger numerical candidate or an optimizer status
was never registered without an exact rational population and geometry proof.

## Verification and preservation

The new packages include independent replays, explicit input scopes, corruption
controls and complete SHA-256 manifests. Search binaries and large reproducible
candidate pools are ignored. Compact explicit certificates, source code,
bounded search logs and mathematical negative results are retained.

The repository registry, current tables, status notes, verification routes and
manuscript-pending metadata are updated together. The manuscript itself remains
the separately preserved snapshot. Run `python audit.py --repository-only` for
repository consistency; each proof package's `verify.py` performs its geometric
and arithmetic replay.
