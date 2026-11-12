# Routes to a substantial improvement

This is a research plan grounded in the November 12 experiments. It does not
assert that any unexecuted optimization will succeed. Here a useful milestone
is a net gain of at least 1,000 over the pinned PackingStar source in one
dimension, followed by a construction that transfers to neighboring dimensions.
The completed gains of 16 and two do not reach that milestone.

## 1. Joint selection of individual points from several radii

This has the strongest immediate connection to the measured obstruction.
Whole fibers can be very expensive even when a few sign choices are cheap.
Different radii may delete the same incumbent points, so solving each radius
alone can miss profitable combinations.

Let X be the current complete configuration, C a supplied set of distinct
new norm-32 candidates, and B(c) the old points with c·x>16. Use binary
variables z_c for selected candidates and y_x for deleted old points:

```math
\max \quad \sum_{c\in C}z_c-\sum_{x\in X}y_x,
\qquad z_c\le y_x\quad(x\in B(c)),
\qquad z_c+z_{c'}\le1\quad(c\cdot c'>16).
```

All old layers remain eligible for deletion. Do not impose antipodality or
require a full sign fiber. An old deletion is charged once, regardless of
how many selected candidates require it. If old points are compressed into
weighted profiles, the profile equivalence must preserve their products
with **every admitted candidate family**; a compression valid for one
eight-coordinate normal space cannot be reused blindly after changing it.

Start with a merged pool containing the following families:

- Individual signed octads with three core and five tail coordinates, two
  core and six tails, and one core and seven tails, using all 24 core positions.
- Norm-24 midpoint bodies paired with new norm-8 tail shapes
  `(±2,±1,±1,±1,±1)`, where the dimension permits them.
- Norm-20 core bodies with norm-12 tetrahedral tails, including actual
  quarter-angle Leech midpoints rather than just signed five-support bodies.
- New norm-8/norm-24 points as complementary variables, despite the standalone
  30D obstruction: its deletion-incidence proof does not automatically exclude
  a profitable union with other families that pay for the same old deletions.
- Leech restorations and ordinary midpoint/reflection candidates exposed by
  each tentative changed configuration.

Price complete small support families, then select neighborhoods by common
blockers and mutual compatibility. Selecting only the individually cheapest
points can throw away a useful high-cost group with nearly identical blockers.
Use candidate overlap with the current deletion set as a second score.

A practical implementation alternates a large-neighborhood heuristic and
exact optimization on the affected incidence component. Maintain the union
of old blockers explicitly; remove a cluster of chosen new points, release
deletions that are no longer needed, and solve or greedily refill a broader
cluster. Preserve all genuine positive witnesses before enlarging the pool.
Use conflict cliques and finite fiber-capacity inequalities to strengthen
the component model. Generate candidate/candidate conflicts on demand instead
of materializing an unrestricted quadratic matrix.

**Acceptance gates.** First seek a completely verified positive exchange,
then a reproducible gain of at least 100, and finally 1,000. A candidate pool
whose merged old/new graph has a matching saturating every candidate is closed
without solving its new/new graph. Certificates for separate pools do not
combine unless their assigned old blockers remain jointly distinct. A
numerical solver bound is a scheduling signal; claim an exact exclusion only
after reconstructing a checkable integer or rational certificate.

This route is best motivated in 29–31D, where larger tail alphabets exist.
The normal-octad whole-fiber pilot is only a starting restriction. Its
unresolved status is not evidence that a larger point-level model succeeds.

**Broader point pilot already executed.** The campaign subsequently priced
726,528 individual points on 42 mixed triples, all 276 core pairs and all
24 singletons. The adaptive 19,110-point pool has 302,518 old/new incidences
and 4,105,565 new/new conflicts. Minimum costs fall to four and two in the
3+5 and 2+6 families, but three bounded annealing/continuation schedules found
no positive exchange. Continuation of the deletion penalty from 0.02 toward
one reached 262 selected points with 1,667 deletions as a trajectory diagnostic,
far from the desired trade. The saved best nonempty exchange still has gain
minus one. Thus simply removing the whole-fiber restriction or heating this
same pool has already been tried. The remaining proposed enlargement combines
the complementary radii and exact incidence-component repairs rather than
repeating that short heuristic unchanged.

The unrestricted three-family signed-octad universe has 11,381,760 candidates.
The source-dependent signed-subset count gives 473,381,256 old/new incidences,
about 1.89 GB for one uint32 adjacency direction before offsets. Enumeration
and pricing are plausible on the current host, but an unrestricted quadratic
new/new conflict matrix and a monolithic MILP are not. These counts and the
executed pilot are recorded in
[`source_and_layers/reports/partial_sign_pilot.json`](source_and_layers/reports/partial_sign_pilot.json).

For the 31D three-family alphabet, exact same-tail capacities are four,
two and one. Summing over the signed tail vectors gives a total capacity
of 3,712 before conflicts between tails. Even if this full capacity could be
attained, a net gain of 1,000 would require deleting at most 2,712 incumbent
points. This is a useful quantitative target for the shared-deletion search.
For any admitted pool C, a bipartite maximum matching of size ν between new
candidates and old blockers gives the further exact bound `gain ≤ |C|−ν`.
This uses a vertex-cover argument on the old/new graph; a large matching
deficiency merely permits the target and does not establish a compatible
selection.

## 2. Use module transplants as changes of geometry

The Barnes–Wall experiment offers a concrete reason to allow temporarily
smaller configurations. Its best forced seed has 251,244 points, ten below
PackingStar, but its 1,996 new heads change a substantial part of the contact
geometry. Retaining it as an independent search seed is reasonable even
though it is not a registered improvement.

The completed free Leech and classical refills, and the new-head contact
closure at blocker cutoffs two and four, gave no gain. Repeating those same
closures is not the next experiment. Broaden the module problem by selecting
individual heads from multiple differently oriented BW modules together with
non-BW auxiliary points, sharing external deletions. Permit old module heads
to survive where compatible; do not force an entire module when partial
replacement can be better.

The measured full-module objective is transparent:

```math
251254-1958-68+2016=251244.
```

Reducing the external cost below 58 would produce a record in this restricted
full-module approach. However, even reducing it to zero gives only 58 points
for one module. A thousand-point target therefore requires multiple jointly
compatible subspace replacements or a larger replacement construction. The
single-module near miss should not be sold as a thousand-point route by itself.

## 3. Replace the support system, with a quantitative target

The whole signed-octad construction makes the target especially clear.
Retain the 98,304 odd Leech vectors, include all 2n(n−1) coordinate roots,
and use 128 even sign patterns on each admissible eight-support. Supports
must intersect pairwise in at most four positions. A support is admissible
if it is a core Golay octad or contains at least two tail coordinates.
Then

```math
N(n)=98304+2n(n-1)+128M.
```

In 31D, 1,181 supports give 251,332, 1,189 give 252,356, and 1,250 give
260,164. These would respectively exceed PackingStar by 78, 1,102 and 8,910.
They are conditional targets. The executed support search reached only 928
supports and 218,948 points, after enumerating all 4,731,285 admissible supports.
The gap to the first useful target is large.

A useful new search would exploit nonlinear support trades or compatible
partial sign families, instead of another greedy whole-support refill.
Partial fibers turn the intersection-four rule into actual signed-product
constraints and can recover points lost by the conservative whole-support
condition. Their count is no longer the simple 128M formula, so publish an
explicit full witness or a complete finite counting proof. The previous
32D parent-shortening experiment also failed to reach the current 31D bound;
it does not provide a hidden shortcut to 1,189 supports.

## 4. For 26–28D, change the latitude geometry

The large-tail octads require at least five tail coordinates and therefore
do not directly address 26–28D. A serious all-dimensions program needs a
second mechanism rather than claiming that a 31D construction automatically
projects without losses.

Two exact local replacements supply useful starting identities. For
norm-32 Leech endpoints with product sixteen, their norm-24 midpoint can
be shrunk and lifted on an equilateral tail triangle, producing three
points after potentially deleting two parents. For endpoint product eight,
the midpoint has norm twenty and supports a four-point norm-twelve tetrahedron
when at least three tail coordinates are available. Other Leech points clear
the midpoint or its shrink by convexity. These identities do not settle the
cross products between different fibers or with existing auxiliary layers.

The fixed-plane triangle pilot in 31D was expensive, and its fixed-phase
code has an exact 280-head bound. Thus a new experiment must genuinely vary
the plane/phase assignments, couple their parent deletions, and allow changes
to the incumbent layers. For 26D, an angular interval arrangement on the
two-dimensional tail circle is a concrete finite subproblem. For 27D and
28D, several tetrahedral or triangle orientations can be admitted jointly.
Price these against each dimension's actual source; the 31D blocker counts
do not establish their costs in lower dimensions.

This direction is more speculative than the point-level mixed-radius
optimization. It is listed because it changes the geometry available in
the dimensions where the large-tail families do not exist. No achieved gain
or unrestricted bound follows from the local three-for-two or four-for-two
identities alone.

## What should not be repeated unchanged

The six literal previous-source unions, free classical completions and
56-point-exchange tail amplifications now have exact finite exclusions or
complete scans. The tested 30D latitude has an explicit five-versus-four
incidence obstruction. The displayed 21 whole-fiber orientation has an
integer flow certificate excluding gain from every subset of its fibers.
These results are useful stopping rules for those exact models. They leave
new frames, partial signs, merged families, shared deletions and different
radii open, with no guarantee that those enlargements will improve a record.
