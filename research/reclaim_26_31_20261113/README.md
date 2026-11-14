# Reclaiming PackingStar dimensions 26–31: 13 November 2026

This campaign establishes **248,084 in 30D** and **251,395 in 31D**, respectively
**92 and 141 above PackingStar's November 11 release**. Together with the
preceding 28D and 29D certificates, the repository now exceeds the pinned
release in four of the six target dimensions. **26D and 27D remain open.**

## Current comparison

The PackingStar source is commit
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
The repository started at commit `30b72dab72b69bcb2760a2cce0ec5edfa0ab56c8`.
The [machine-readable comparison](../../common/source_review_20261113_26_31.json)
records the source identities and scope. This is a comparison with the checked
release, not an exhaustive claim of worldwide priority.

| Dimension | PackingStar | Pre-campaign repository | Current repository | Margin over PackingStar |
|---:|---:|---:|---:|---:|
| 26 | 201,180 | 201,136 | 201,136 | −44 |
| 27 | 209,862 | 209,825 | 209,825 | −37 |
| 28 | 222,104 | 222,120 | 222,120 | +16 |
| 29 | 231,774 | 231,776 | 231,776 | +2 |
| 30 | 247,992 | 247,972 | **248,084** | **+92** |
| 31 | 251,254 | 251,200 | **251,395** | **+141** |

An earlier conversation reported a 251,256-point 31D checkpoint, but no
coordinate certificate or published commit was recovered in the checked main
tree. It was not used as an established result. The new 251,395-point package
has its own complete explicit certificate and proof.

## Successful constructions

### 30D: partial fibers in a 1,356-point section

The final construction replaces 932 normal-module heads by 1,056 heads and
deletes 32 exterior blockers:

**247,992 − 932 − 32 + 1,056 = 248,084.**

The 1,356-point module is a coordinate section of the reassembled 2,564-point
Leech–Sloane configuration. It shares the incumbent's 240-point E8 equator and
60-point D6 tail set. An exact R/4 orientation aligns it with the surrounding
layers. The minimal exchange has 964 deletions and 1,056 insertions.

The 1,056 heads use 192 fibers of size two, 48 of size ten, and 12 of size
sixteen. This bypasses a proved 960-head ceiling for complete sixteen-point
fibers in the same integer tail alphabet. The gain comes from allowing partial
fibers, not from irrational coordinates.

The [portable proof](../../verifications/improved/dim30-leech-sloane-section/)
checks all **30,772,711,486 final unordered pairs**, independently reconstructs
the parent and section, and checks every inserted-point comparison directly
with integers. The final configuration contains 1,216 unpaired points.

Earlier in the same campaign, a 1,260-point Barnes–Wall section with no exterior
deletions gave **248,020**. Its
[complete certificate](../../verifications/improved/dim30-bw14-module/) is
retained. A matching certificate proves that the literal union of the two
fixed 30D witnesses cannot improve 248,084.

### 31D: the larger module and an individual-point union

The normal 15D sector grows from 2,282 points to a 2,564-point candidate module:

| Component | Incumbent sector | New candidate module |
|---|---:|---:|
| E8 equator | 240 | 240 |
| Mixed heads | 1,958 | 2,240 |
| D7 pure tails | 84 | 84 |
| Total | 2,282 | 2,564 |

The key was to rebuild the known 15D recipe with a coordinate E8 equator.
Seventy transversal quadruples plus twenty-one pair-union quadruples form a
Steiner quadruple system on fourteen coordinates. Two explicit affine binary
codes supply the glue layers. The complete module is independently reconstructed
and every pair is checked; no solver feasibility claim substitutes for geometry.

The first complete-module replacement gave 251,358. A bounded search of 3,530
orientations, evaluating each two-copy union exactly, improved this to
**251,395**. It retains 2,559 module points and removes 136 exterior points:

**251,254 − 2,282 − 136 + 2,559 = 251,395.**

After cancelling coincidences, this is 2,219 additions and 2,078 deletions.
The exact orientation is R/2. Five otherwise new heads are omitted, and the
final configuration has 61 unpaired points. The
[portable proof](../../verifications/improved/dim31-leech-sloane-module/)
checks all **31,599,597,315 final unordered pairs**, all 3,285,766 module pairs,
and 557,843,286 inserted-point comparisons by an independent integer method.

## PackingStar's actual restriction

The general method is not restricted to rational coordinates or rational Gram
matrices. The latest reviewed paper is
[arXiv:2511.13391v5](https://arxiv.org/html/2511.13391v5). Its rational-Gram
restriction concerns a separate 13D experiment. Methods equation 18 uses a
finite basis-product alphabet C1, while equation 19 permits a continuous C2.
Structural membership restrictions can additionally limit the admitted family.
The current 26D archive itself includes sixteen points over Q(√3).

An irrational global rotation changes neither the Gram matrix nor the kissing
geometry. The useful attack changes the admitted module, partial-fiber structure,
relative orientations, or actual Gram values. Both successful constructions in
this campaign have integer coordinates before unit normalization.

## Other completed experiments and their limits

| Work | Outcome | Scope |
|---|---|---|
| All 759 Golay-octad transplants in 26D | No improving tested two-copy union; 2,004 cheap candidates retained | One deterministic ordering per octad; combinations remain open |
| Joint 2,004-candidate 26D optimization | All source blockers replayed exactly; solver process hit its 160-second wall limit without returning a selection | No exclusion or optimum follows |
| Exact 438-point 11D module with E8 equator in 27D | Valid local geometry; 256 axis orientations and 1,024 rational rotations gave no improving tested union | Exterior costs were at least 8,000 or 8,384 in the respective pilots |
| Continuous 26D, rational tail and latitude pilots | No feasible improvement found | Numerical failures and numerical finite-pool optima are not global upper bounds |
| Fourteen direct BW14-to-BW15 extensions in 31D | No improving tested union | Those fourteen extensions only |
| 31D missing-axis catalogue | 8,008 exact new points generated; 3,960 with at most eight blockers | Entire admitted union has exact gain at most zero by capacity and matching |
| 31D three-axis rational rotations | 478,560 tested states; best four-for-six bundle had 90 external blockers | Does not exclude coordinated changes of several fibers |
| Irrational seven-cube tail family | All 276,480 candidates priced; every one has at least 21 blockers | A separate geometric proof limits a full replacement to 1,920 heads; mixed unions are not excluded |
| 30D full-head catalogue and partial-fiber heuristics | No gain over the initial BW14 seed in the tested searches | Later explicit Leech–Sloane section demonstrates that the unrestricted partial-fiber model was not excluded |
| Eight sextet frames, eight core reflections and twelve mixed image trials in 31D | No improving tested union | Only those sampled images and blocker cutoffs |

The initial 26D latitude pilot contained an int16 intermediate-overflow bug.
Its invalid outputs were discarded, arithmetic widened before multiplication,
and norm checks added. Only the corrected outputs are retained. This correction
did not affect either successful 30D/31D construction or their independent proofs.

## Directory map

| Path | Contents |
|---|---|
| [packingstar_audit/](packingstar_audit/) | Source audit, initial BW14 improvement, direct extensions, full-head searches, exact restricted 960-head ceiling |
| [algebraic/](algebraic/) | Final 31D construction, portable recipe, exact-cover and orientation searches, first witness, cube and rotation analyses |
| [continuous26/](continuous26/) | Final 30D section, lower-dimensional pilots, 26D catalogue, 27D module and source pricing |
| [bulk_fibers31/](bulk_fibers31/) | Exact missing-axis pool, interrupted initial solve, independent capacity/matching exclusion |
| [joint26_pool/](joint26_pool/) | Exact blocker replay and bounded numerical selection attempt |
| [SEARCH_26_27.md](SEARCH_26_27.md) | New component-wise 26D optimizer and full-alphabet 27D partial-fiber search |
| [STRATEGY.md](STRATEGY.md) | Prioritized remaining attack with explicit objectives and stopping rules |

The scripts at this directory's root generate the missing-axis, joint-26D,
core-reflection and mixed-image experiments. The source is normalized from its
declared squared norm; 31D uses stored norm128 coordinates, most other rational
files use norm32, and 26D has a separate quadratic format.

All heavy jobs respected the available CPU affinity, reserving at least one
logical CPU. No search jobs remain active at the end of the campaign. The
research packages preserve their seeds, exact inputs, search scopes and proof
receipts; only the portable verification packages register lower bounds.
