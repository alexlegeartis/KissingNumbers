# Exact orbits, inherited sections and stabilizer-aware screening

These additions are dimension-independent instruments. They do not establish
universal competitiveness, new records, abstract group orders or irreducibility.
There are38 registered operators, all with planner input/output signatures.

## Exterior-square construction

`groups.exterior_square` takes any existing `unitary_representation`. In the
basis indexed by pairs `i<j`, it forms the matrix of two-by-two minors:

`wedge(R)[(i,j),(k,l)] = R[i,k]*R[j,l] - R[i,l]*R[j,k]`.

This representation feeds the ordinary group candidate search. It remains
numerical and cannot promote floating-point generators or candidates to exact
geometry. Output dimension and matrix allocations are bounded before source
array conversion. This permits general representations beyond the existing
Fourier/Weil parameter catalogue.

## Exact orbit construction and proof

`groups.integer_orbit` accepts a nonzero integer seed and exact rational
orthogonal generators, each encoded as an integer numerator matrix and a
positive denominator. It performs actual orbit generation, removes duplicates,
and rejects a full-orbit proposal as soon as a distinct generated point violates
the seed inequality. That rejection does not rule out partial-orbit selections.
A nonintegral image means the supplied coordinate scale is unsupported, not
that the mathematical orbit is impossible. A point/work cutoff is incomplete.

The independently replayed certificate has backend `exact_integer_orbit`.
It checks integer norms and distinctness, every generator's exact orthogonality,
every point's image under every generator, full reachability from the seed,
seed separation, and exact span rank. It does not trust saved closure flags.
Closure alone is insufficient: a disjoint union of closed orbits is rejected.

For any two orbit points, move the first to the seed with an inverse word.
Exact orthogonality preserves their product, while closure ensures the image
of the other point is still stored. Thus checking `N-1` seed comparisons covers
all `N*(N-1)/2` pairs. Reports distinguish comparisons computed from pairs
covered by this argument. Stabilizers and duplicate group images are allowed.
Full span additionally proves finiteness of the generated matrix group because
its action on a finite spanning set is faithful. Lower-dimensional orbits do
not establish finiteness on their ambient orthogonal complement.

Point, orbit, graph, transform-work and integer-size budgets are enforced.
Block transforms use int64 only with a proven partial-sum bound; larger inputs
use Python integers. A small spanning prefix suffices to prove full ambient
rank; deficient prefixes fall back to exact rank of the complete set.

## Central sections without a new pair sweep

`geometry.section` accepts a batch `normals` parameter as an alternative to
its existing `normal`. It also recognizes exact orbit/section sources.
The output preserves the metric and all source coordinates, selecting all and
only points satisfying the exact rational central-plane predicates.

Backend `certified_coordinate_section` independently replays the source proof
and checks complete subset membership, copied metric/norm, and exact rank.
Separation and distinctness follow from the verified source subset. Rank uses
an exact lower witness and the exact source/normal upper bound, with a complete
elimination fallback where needed. Nested sources have a depth guard. This
allows large sections without silently spending a quadratic pair budget.

The Leech control generates196560 points from four exact Golay/sextet actions.
It checks786240 generator images and196559 seed comparisons, covering
19317818520 unordered pairs. A uniform signed-pair-normal search then obtains
93150/rank23 and49896/rank22. Its next outputs27720/rank21 and17400/rank20 are
below the captured targets31560 and19448. These are known-control and bounded
section-search outcomes, not new SOTAs.

## Exact class-average screening

The classical class-sum argument uses Schur's lemma; see
[Bump's representation-theory notes](https://sporadic.stanford.edu/bump/group/gr2_4.html).
For a real irreducible orthogonal representation, use the symmetric part of
the class sum. It is self-adjoint and commutes with the representation, so its
eigenspaces are invariant and it is scalar. Its scalar average is `a=chi(g)/d`.
This argument also covers real representations of complex or quaternionic type;
the nonsymmetric class sum itself need not be a real scalar.

If a class has `m` elements and exactly `s` fix a proposed seed, full-orbit
kissing requires `m*a <= s + (m-s)/2`, hence

`s >= ceil(m*(2*a-1))`.

The strict-character-ratio obstruction therefore needs the stabilizer count.
At `a=1/2`, all class constraints must be tight only when `s=0`.

`group_constraints.analyze_orthogonal_actions` verifies an explicit rational
orthogonal action list, its symmetric average, exact fixed images and exact
violating distinct images. It uses no assumed character table, conjugacy-class
membership, irreducibility or group finiteness. A scalar-average identity is
checked directly. Its obstruction is explicitly for the full generated orbit;
partial point pools are not pruned.

In the E8 exterior-square example the120 root reflections have exact average
one half. Their root outer products span all36 symmetric8-by-8 matrices.
In the sector with no fixing root reflections, equality forces `A^T*A` scalar,
so a nonzero skew seed is a scaled orthogonal complex structure.

Stabilizers cannot be omitted. Independently constructed exact whole-orbit
controls include rank-two skew seed7560 points with30 fixing reflections,
rank-four37800 with fixing reflections, and full complex-structure15120 with
none. The rank-two seed has nonscalar `A^T*A` and is a concrete exception to an
unqualified complex-structure conclusion. The three-block seed already has a
distinct-image reflection violation. No partial-orbit ceiling is inferred here.

For sign changes in `B_n` acting on the exterior square, the average is `1-4/n`.
For `n>8`, a full kissing orbit needs class stabilizers. That condition does not
exclude every orbit or any partial orbit. All these28D controls are below221988;
the experiment is bounded and subsequent work concentrates on other dimensions.

## Weighted quadratic coordinates

`geometry.quadratic_code` accepts explicit coefficients `A+B*sqrt(3)` and an
optional positive integral diagonal metric. The native checker verifies rational
squared norms and every normalized pair inequality exactly. This is a canonical
input route, not an autonomous candidate generator. The known18D8358 adapter
uses weights `(1,...,1,2,2)` and squared norm288.

The Q(sqrt3) backend vectorizes sign/squaring cases only when proven bounds cover
the second squaring as well. Larger intermediates retain the scalar Python
integer path. Temporary Gram blocks adapt to the declared matrix-entry budget.
Rank is unchanged by the positive metric. Random weighted pairs are tested
against the independent scalar sign oracle, including overflow fallback cases.
