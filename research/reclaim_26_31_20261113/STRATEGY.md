# Next attack after the verified 30D and 31D gains

The immediate unsolved targets are **201,181 in 26D** and **209,863 in 27D**:
one point above the pinned PackingStar configurations. Use those stronger
external configurations as working seeds. The smaller published repository
values are historical constructions, not the starting point for search.

The successful mechanism is now concrete: a larger local module, an embedding
that preserves the actual E8 equator, and individual-point selection against
the surrounding layers. The next work should exploit that mechanism while
changing the construction where the measurements show a barrier.

## 1. Optimize the new design against its exterior

The new 15D module uses one solution of a small exact-cover problem. Selecting
the design jointly with exterior deletion cost is more targeted than choosing
another arbitrary solution and searching only its orientation.

For a fixed pairing and orientation, let z_b select one of 560 transversal
quadruples. Each of the 280 triples using three distinct pairs must belong to
exactly one selected quadruple. These constraints choose seventy quadruples;
the other twenty-one pair-union quadruples remain fixed. Each selected support
contributes sixteen signed points. The fixed affine glue codes, axes and outer
points complete the construction.

Price each block against the actual exterior. A binary deletion variable for
an old exterior point is charged once even when several selected blocks
conflict with it. Require its deletion whenever a selected block conflicts.
Include every deletion forced by the fixed part of the module. The first identity placement required 178 exterior deletions. The final
orientation has 142 blockers for the complete module, reduced to 136 by
omitting five heads. When varying the design, recompute the complete blocker
union across all candidate supports rather than restricting deletion variables
to the blockers of one selected design. This gives a design-and-deletion model
before allowing finer point selection.

Every chosen exact cover still needs reconstruction and an independent full
module pair check. Then solve the old/new union at point level: the final
31D witness improved by omitting five heads, so complete-block or antipodal
selection must not be imposed on the final exchange.

The same design controls its 1,356-point 30D sections. A joint objective can
favor both dimensions, or maintain a small set of nondominated designs rather
than optimize one scalar score. This is an unrun continuation, not an asserted
improvement or an optimality claim for the present design.

## 2. Attack 26D through shared deletions and coupled deformations

**14 November update:** the original 2,004-point transplant pool now has an
[independently checked exact no-gain certificate](../26_27_exchange_exclusions_20261114/).
The continuation also certifies zero possible gain for a merged 4,662-point
26D pool and a merged 5,230-point 27D pool. All source points remain eligible
for deletion. These are finite-pool exclusions; the broader deformation
and candidate-generation questions below remain open. The interrupted solve
described next is the historical November 13 run.

The complete 759-octad sweep supplies 2,004 distinct rational candidate points
with two to four old blockers each. Their old-point blockers were replayed
against the full quadratic source. The joint model has 4,094 distinct old
blockers and 8,724 candidate/candidate conflicts. The first bounded MILP did
not return before its external wall limit. It did not prove the pool exhausted.

The objective must count shared deletions:

**gain = number of admitted new points − size of the union of their blockers.**

A useful continuation decomposes the incidence graph where possible, uses
integer matching or capacity bounds to close components, and searches the
remaining components with explicit feasible incumbents. An empty exchange is
always feasible. Save every positive exact exchange immediately; a solver's
time-limit message must not erase a witness or become an exclusion claim.

The normal 10D sector itself contains 380 points. Sixteen radical points have
few local contacts, but attempts to replace them by seventeen and broader
381-point numerical searches did not find a valid improvement. Their apparent
freedom alone is insufficient evidence that another point fits.

The fixed exterior compresses to 2,300 nonzero norm16 profiles of multiplicity
32, 240 norm8 profiles of multiplicity512, and one zero profile of multiplicity
4,320. Violating a single active exterior profile can therefore be expensive.
Build deformation cavities that include the points inducing the binding
profiles, rather than move a small normal-sector cluster against a completely
frozen exterior. Preserve many inherited contacts by exact rational conic or
orthogonal parametrizations and optimize the remaining variables numerically.
Only an exact reconstructed witness should advance the registered bound.

## 3. In 27D, expand the coupled region before enlarging the module

An explicit 438-point 11D code with the required E8 equator is already available.
Its nominal sixteen-point advantage over the incumbent 422-point sector is
overwhelmed by exterior conflicts in every tested orientation. Whole-module
transplantation with that frozen exterior is therefore a poor next use of
compute unless a new alignment argument substantially reduces those costs.

Identify which exterior profile classes block the most useful partial
submodules. Optimize several fibers and their shared old-point deletions
together. Keep the old source available point by point; requiring either the
entire old sector or the entire new sector excludes the kind of asymmetric
gain that succeeded in 31D.

The lower coordinate sections of the new 15D module have only 308 points in
the 26D setting and 348 in the 27D setting, below the current normal sectors.
Repeating those coordinate sections cannot help. A different section, norm
split, body embedding or larger coupled region is needed.

## 4. Use the exact exclusions to avoid repeating closed models

Three restrictions have actual proofs in this campaign:

- Complete sixteen-point core fibers over the fixed 30D integer tail alphabet
  contain at most 960 heads. The successful 1,056-head section escapes by using
  fibers of sizes two, ten and sixteen.
- The particular seven-cube tail alphabet in 31D supports at most 1,920 heads
  in a full replacement, even if the core coordinates are arbitrary real
  numbers. It cannot replace the incumbent 1,958 heads profitably by itself.
- The pinned 31D source together with the saved 3,960 missing-axis candidates
  has no improving subset. Saturated sixteen-dimensional fibers account for
  all 3,928 axis candidates; the remaining thirty-two have distinct outside
  blockers. This is a finite-union proof, not a bound on other candidates.

The full-coordinate results demonstrate why these hypotheses matter. A
restricted negative result should guide a change of family, not be described
as a global barrier.

## 5. Continuous coordinates remain useful when they change geometry

PackingStar's general framework permits irrational coordinates and continuous
admissible products, but actual searches still use finite alphabets and selected
structural families. Seek new Gram values or relative placements, not a global
irrational rotation of an unchanged configuration.

One concrete exact deformation was derived for a residual 31D fiber. Given
orthogonal length-four axes p, q and a vacant axis r, the three vectors
(8p+q−4r)/9, (p+8q+4r)/9 and (4p−4q+7r)/9 are another orthogonal frame.
Their antipodes replace four old points by six. Projection perpendicular to r
lies in the convex hull of the old parents, preserving compatibility with every
unchanged point having zero r-coordinate. The current pilot still incurs at
least ninety exterior blockers. Future work must deform the coupled blocker
fibers as well, rather than repeat the same isolated four-for-six trade.

For numerical exploration, track the active contact graph and the rank of
its tangent constraints. Use flexible components to propose larger moves,
then reconstruct rational or quadratic coordinates with rigorous signs. A
numerically small angle violation is never rounded into a certificate.

## Verification and stopping rules

Every positive candidate needs exact norm and distinctness checks, exact
changed-pair checks, a complete final-configuration replay, pinned source and
delta hashes, and a deterministic construction receipt. Mixed quadratic files
must retain their denominators and radicands. Integer intermediates must be
widened before multiplication. Floating matrix multiplication is acceptable
only where every intermediate is a bounded exactly representable integer.

Keep discovery separate from proof: optimizers may propose point selections,
but their statuses do not establish kissing geometry. Publish only fully
replayed witnesses. A completed failed heuristic, an interrupted solver and
an exact finite-family exclusion are three different outcomes and must be
recorded as such.

Give each new model a bounded initial budget and a useful checkpoint format.
If repeated runs pay the same exterior costs, change the coupled region,
design, norm split or admitted alphabet. Reserve one logical CPU and cap the
combined thread budget across workers. The present campaign has stopped its
search jobs; none of these proposed continuations is running in the background.
