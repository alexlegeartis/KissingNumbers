# Discrete exchanges in 27 dimensions

Takhanov's **210,108-point** configuration is the benchmark. This campaign
replayed its full geometry and did not find a larger configuration. The
coordinates have squared norm 128 and the largest inner product of distinct
points is 64, over all 22,072,580,778 unordered pairs. The two input archives
are pinned to the same published revision. [Takhanov's source and certificates](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1/certificates)

The continuation gives a 210,108-point witness with **424,418,622 touching
pairs** and a family with at least **283 isometry classes**. The exact
exchanges and their scope are recorded below.
[Coordinates](data/d27_midpoint_switch.npz) · [Midpoint switching and fibre exchanges](#switching-midpoint-tails-and-exchanging-native-fibres)

| Finite search | Completed result |
|---|---|
| Signed native midpoint exchanges, with common residuals in each class of 32 fibres | All 133 one-class and 8,778 two-class windows excluded |
| Full integer grid in the 11D kernel | 5,377,994 antipodal lines clear the fixed exterior. Every one has a conflict with an incumbent kernel line |
| E8 midpoint and head-frame exchanges | The degree bound excludes improvement for all 5,702,400 frame triples in the catalogue |
| E8 midpoint refill after omitting incumbent head lines | All 190,051 choices omitting zero to six lines excluded |
| Two-witness fusion under Golay sign changes and affine translations fixing the octad | Perfect matchings certify zero gain in all 512 unions |
| Pure Barnes–Wall midpoint refill | All 604,800 contacts enumerated. The extra-fibre construction leaves 6,144 admissible tail pairs on 512 endpoints and needs a refill gain of 417. The finite optimization did not reach that target |
| Barnes–Wall auxiliary frame changes | 312 configurations of 210,108 points. Twenty distinct contact totals certify at least 20 isometry classes |
| Every Leech contact with arbitrary real tail | All 452,088,000 contact midpoints are blocked by the source midpoint layer alone |
| Contact faces on deleted roots | All 2,966,566 whole contact chords, 137,159,728 triangle centroids and 2,647,762,920 tetrahedron centroids are blocked by the same layer for arbitrary real tails |

Each exclusion applies to the stated finite family. The Barnes–Wall
optimization ended without an exclusion certificate. Its best recovered
refill had gain 256 and gave 209,948 points with the extra class. On the
original fibre arrangement the same enumeration admits 28,544 tail pairs on
960 endpoints, and the best refill found ties the source gain of 512.

## Native fibres and exact exchanges

Write each coordinate vector as `(w, y)`, with `w` on the complement of
the reference octad and `y` on the octad and three tails. Each class contains
32 Barnes–Wall projections and a common E8 frame. Two classes with
nonuniform source states are frozen. All 16 retained Leech options and the
signed midpoint options of a free class are generated from its frame.
Distinct sign choices in the tail are independent decisions.

For two different classes the pair test is
`y · y' + max(w · w') <= 64`, with the maximum over the two projection
classes. For a pair of residual choices in the same class, the maximum
projection term is 64. Repeating a single residual at different projections
of its own class is compatible. The branch-and-bound solver first filters
outside conflicts, then exhausts the finite graph in each window. One
additional selected residual would add 32 points.

The kernel search fixes all points with nonzero complement projection. It
uses the 2,895 distinct exterior residuals as linear inequalities and
enumerates integer vectors of squared norm 128. A sum-of-squares table
prunes prefixes that cannot reach the sphere. Constraints supported on
already assigned coordinates are applied before extension.

## E8 degree bound

The frame catalogue on each tail axis contains all compatible native E8
frames and the compatible four-axis Hadamard switches of that axis's
incumbent frame. The three counts are 144, 275 and 144. With the exterior
fixed we generate every contact midpoint in the E8 roots, with the six
tail lines and independent signs of body and tail.

Take a signed Leech endpoint `r` and two other endpoints `s` and `u`
at inner product 64 with `r`. The midpoint body product is
`64 + s·u/4`, and `s·u >= -64` by orthogonal projection off `r`.
Using the same tail line twice is impossible. Distinct nonorthogonal
tail lines cost 16 in the worst sign combination. Compatibility in that
case requires `s·u = -64` and hence `s + u = r`. Each tail line has
only one perpendicular partner, so among three lines one would require
the same opposite partner endpoint twice. Thus each endpoint has degree
at most two.

With antipodal bodies, an orbit adds four points and each used root line
costs two points. The net midpoint contribution is twice the number of
root lines of degree two. Improving the source while keeping three
full head frames and the poles therefore requires 65 such lines. For
each frame triple, count root lines incident to a pair of mutually
compatible eligible midpoint orbits. This upper bound is at most 64
throughout the catalogue.

Omitting `k` incumbent head lines loses `4k` points, so the midpoint
refill needs at least `65 + 2k` degree-two root lines. The incidence bound
leaves 76 cases when `k` is five or six. In each, a matching of disjoint
pairs of incompatible root lines forces enough omissions to close the
gap. Two root lines are incompatible when every allowable pair of
midpoint choices for one conflicts with every allowable pair for the
other.

## Two-witness fusion and the refill

For two valid codes, conflicts after keeping their common points occur
only between the exclusive parts. Each matching edge forces at least
one deletion, and the corresponding minimum vertex cover attains that
bound. The fusion checker verifies every matching edge by integer
arithmetic.

The saved refill anchor activates one previously unused native class
with 256 midpoint points and deletes 128 Leech endpoints in it. Some
existing row and pure Barnes–Wall midpoints are displaced. The refill
generator enumerates every pure Barnes–Wall contact and all six tail
lines against the remaining fixed additions. Each surviving midpoint
pair blocks exactly its two Leech endpoints. The objective pays an
endpoint only once and selects pairs subject to exact geometric
conflicts. Any candidate that reaches the target is checked over every
unordered pair before it is saved.

## Arbitrary real tails in the Barnes–Wall layers

For a fixed pure Barnes–Wall midpoint body `h`, group the exterior
points by their signed tails. The complement projections at `t` and
`-t` agree, so a candidate tail `y` must satisfy
`|t · y| <= 64 - max(h · b)`. A negative right hand side rules out
the body. Otherwise the coordinate-axis bounds make the feasible
polytope bounded. Every vertex comes from three independent active
inequalities and is enumerated with integer numerators and a positive
determinant.

In both frozen exteriors, every vertex has squared norm at most 32.
Strict convexity of squared norm forces every feasible length-`sqrt(32)`
tail to be one of the vertices of squared norm 32. All those vertices
are the existing discrete palette. This verifies completeness of the
28,544 and 6,144 tail-pair catalogues even for arbitrary real tails.

In the extra-class refill, the surviving tail lines are perpendicular.
For any feasible signed refill, split the points into two sets by
choosing one sign on each line for the first set and the opposite sign
for the second. Each set can be doubled by tail reflection. Reflected
points from different lines keep their inner product, and those from
one line have a smaller inner product across the reflection. If the
sets use endpoint sets `U_1` and `U_2`, the original refill pays
`|U_1 union U_2|`, which is at least the average cost of the two
doubled refills. At least one doubled refill therefore matches or
improves the signed refill's net gain. This closes the proposed
sign-asymmetry advantage for this anchor.

The reflection argument also applies to the original refill. Orient its
four lines `(0,4,±4)` and `(4,0,±4)` toward the positive third coordinate.
All pairwise products of the representatives are nonnegative. Reflecting
one half then only decreases cross-reflection products, and the same
endpoint-cost inequality applies.

The same vertex test was applied with squared body and tail norms 64
after discarding the pole constraints. It excludes any new head in
three tested body families: half-sums of orthogonal pure Barnes–Wall
roots, half-sums of orthogonal pure E8 roots, and bodies with four entries
`±4` meeting both octad blocks. The fixed midpoint and auxiliary
layers alone exclude these candidates.

Exploratory runs also tested contracted midpoint bodies with equilateral
or tetrahedral tails, native fibre states built from oriented triangles,
and axial midpoint tails with a square-root coefficient. Those searches
gave no improvement. The replayable exclusions above concern the stated
frozen exteriors and body families.

## Equal-count constructions and the enlarged refill

Each of three auxiliary axes carries both signs of a tail of length eight.
At each signed tail the complement bodies are the 32 signed vectors in one
orthogonal frame. With everything else fixed, the 135 native frames give
13 choices on the first axis, 24 on the second, and one on the third.
Between different axes the tail product is zero, so all 312 triples work.

The number of unordered pairs with inner product 64 distinguishes at least
20 isometry classes among these triples. The part of that count involving
the 192 auxiliary points ranges from 126,240 to 127,648. Takhanov's triple
has 126,624. All pairs inside the common part contribute equally.
The [minimum-contact example](data/d27_head_frame_tie.npz) ties the bound
and has been verified over all 22,072,580,778 unordered pairs. An independent
count over the full coordinate matrices gives 424,419,726 touching pairs in
this example and 424,420,110 in Takhanov's source.

On the separate extra-class construction, the frame triple (126, 1, 134)
retains the 209,948-point starting arrangement and admits 9,216 midpoint
tail pairs on 640 endpoints. The original extra-class choice admitted
6,144 pairs on 512 endpoints. The new catalogue still requires a net
refill of 417 to improve the 27D bound. The search found no such refill.

## Exhaustion of contact midpoints

Freeze only the 26,496 signed source midpoint points. They occur in pairs
`(m, ±4n)`, where the six possible `n` have two nonzero entries of magnitude
one. For a proposed body `(r+s)/2` and tail `2u`, let
`a = r·m/16` and `b = s·m/16`. Compatibility requires
`|u·n| <= 8-a-b` and `u·u = 8`. All right sides are integers.
The sphere gives `|u·n| <= 4`, so six integers from zero to four suffice.

For each vector of six bounds, the intersection of the twelve halfspaces
is a bounded polytope containing the origin. It reaches the required sphere
exactly when some vertex has squared norm at least eight. The input builder
checks triples of signed planes by integer adjugates and finds 85
componentwise minimal feasible profiles. The C++ verifier checks profile
feasibility for all 452,088,000 unordered pairs with `r·s = 64`.
None admits a tail. This includes every contact in the entire Leech shell
and uses no constraints from heads, poles, or retained equatorial points.

## Other bodies on deleted contact faces

For a clique of `k` deleted roots with mutual product 64, its centroid
has squared norm `64+64/k`. Write the tail as `4z/k`. Then
`z·z = 4k(k-1)` and `|z·n| <= 4k - sum(r·m/16)`, summing over the
clique. All vertices of the tail polytope have half-integral coordinates
because every nonsingular three-plane determinant has magnitude one or
two. The verifier tabulates feasible bounds and intersects sorted neighbor
lists to visit every triangle and tetrahedron exactly once.

The chord test allows a body `lambda*r + (1-lambda)*s` anywhere between
two deleted contacting roots. Each endpoint has an incident source
midpoint whose body product there is 96 and whose product at the other
endpoint is nonnegative. The paired tails force `1/3 <= lambda <= 2/3`.
The nonnegative product follows from `s·a >= -64` for neighbors `a,s`
of `r`, by orthogonal projection off `r`.

With tail `4z`, the required squared norm is `8*lambda*(1-lambda)`.
For each tail line the checker selects concrete source constraints tightest
at the chord midpoint among the tested rows. Ties supply one affine
inequality on each half of the chord. These necessary constraints give
1,413 patterns over the 2,966,566 edges. Three active planes define an
affine candidate vertex `z(lambda)`. Its feasibility interval has
rational endpoints. Squared vertex norm minus the required squared norm
is convex, so checking both endpoints gives an exact exclusion.

Keeping the source's 12 poles along with its midpoint layer gives a
classification over **all 87,748,128 segments** joining two deleted
roots, including non-contacting pairs. For endpoints `r,s`, put
`g = r·s/32` and use the same tail `4z`. The required squared norm
of `z` is `(16-4g)*lambda*(1-lambda)`. The poles add `|z·n| <= 2`
in all six directions. The two incumbent midpoints incident to each
endpoint bound the interval of possible weights. A subset of the source
inequalities reduces the scan to 50,452 segments. For these, the verifier
builds the full upper envelope of the affine body products for each tail
line. Between intersections of those lines, the vertex test above is exact.

Only 9,684 segments reach the sphere. In every case the endpoints are
orthogonal, the weight is exactly one half, and the tail is a signed
coordinate axis of length eight. The verifier checks all feasible vertices
at both ends of every interval and checks for any strict excess above the
required radius. Strict convexity excludes intermediate weights and
nonvertex points once all equality cases are classified. For each of the
six possible signed tails, bodies must have nonpositive mutual products
in 24 dimensions, so at most 48 fit by the orthoplex bound. The source
already has 48 on each signed tail. Replacing its axial heads by any
collection from these segments therefore cannot improve the count while
the poles and midpoint layer stay fixed.

## Replay

Use Python 3 with NumPy and SciPy, and a C++17 compiler with OpenMP. The
optional refill solver uses OR-Tools. Run from this directory. The
CPU sets below reserve one core on an eight-core machine. Adjust them
to match the available affinity.

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 python source.py
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python fibres.py --width 1
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python fibres.py --width 2
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python fusion.py --workers 7
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python bw_refill.py --extra-class --out /tmp/kn27-refill.npz
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python bw_tail_rigidity.py --workers 7
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python bw_tail_rigidity.py --extra-class --workers 7
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python head_tail_rigidity.py bw --workers 7
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python head_tail_rigidity.py e8 --workers 7
taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python head_tail_rigidity.py mixed_sparse --workers 7

g++ -O3 -std=c++17 -fopenmp kernel_grid.cpp -o /tmp/kn27-grid
g++ -O3 -std=c++17 -fopenmp e8_degree.cpp -o /tmp/kn27-e8
taskset -c 0 env OPENBLAS_NUM_THREADS=1 python kernel_input.py /tmp/kn27-kernel.txt
taskset -c 0-6 env OMP_NUM_THREADS=7 /tmp/kn27-grid 2 /tmp/kn27-kernel.txt
taskset -c 0 env OPENBLAS_NUM_THREADS=1 python e8_input.py /tmp/kn27-e8.bin
taskset -c 0-6 env OMP_NUM_THREADS=7 /tmp/kn27-e8 /tmp/kn27-e8.bin
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python e8_omissions.py

taskset -c 0-6 env OPENBLAS_NUM_THREADS=1 python bw_select.py /tmp/kn27-refill.npz \
  --workers 7 --seconds 240 --out /tmp/kn27-candidate.npz

taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python \
  ../26_27_native_projection_obstructions_20261118/verify_exact.py data/d27_coordinates.npz
```

The November 22 runs can be reproduced from the same directory. The
full contact scan uses approximately 4.5 GB of scratch space and 4.5 GB
of RAM. The face scan uses a separate scratch directory.

```sh
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python head_frames.py --out /tmp/kn27-tie.npz
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python \
  ../26_27_native_projection_obstructions_20261118/verify_exact.py /tmp/kn27-tie.npz
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python bw_refill.py --extra-class --head-frames 126 1 134
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python contact_tail_input.py /tmp/kn27-contacts
g++ -O3 -std=c++17 -fopenmp contact_tails.cpp -o /tmp/kn27-check-contacts
taskset -c 0-6 env OMP_NUM_THREADS=7 /tmp/kn27-check-contacts /tmp/kn27-contacts
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python face_input.py /tmp/kn27-faces
g++ -O3 -std=c++17 -fopenmp deleted_faces.cpp -o /tmp/kn27-check-faces
taskset -c 0-6 env OMP_NUM_THREADS=7 /tmp/kn27-check-faces /tmp/kn27-faces
taskset -c 0 env OPENBLAS_NUM_THREADS=1 python chord_bounds.py /tmp/kn27-faces
g++ -O3 -std=c++17 -fopenmp all_chords.cpp -o /tmp/kn27-all-chords
taskset -c 0-6 env OMP_NUM_THREADS=7 /tmp/kn27-all-chords /tmp/kn27-faces
```

`results.json` records the completed replays. Integral matrix products
used in the searches are exact in binary32 because the coordinates,
products and absolute partial sums fit in its integer precision.


## Switching midpoint tails and exchanging native fibres

The best verified size in the continuation is 210,108. The coordinate
witness in `data/d27_midpoint_switch.npz` has 424,418,622 touching
pairs. Direct enumeration gives 424,420,110 for Takhanov's witness.
The combined family realizes 283 distinct contact totals, each of which
distinguishes an isometry class.

Take the 13,248 different midpoint bodies with their paired tails.
Within the six D3 tail lines, two different lines in the same coordinate
plane are perpendicular. Lines from different planes contribute 16
to the largest signed tail product. A body product above 48 consequently
forces opposite lines in the same plane. The 14,144 such edges give 3,142
connected bipartite components. A body product in (32, 48] forbids
assigning both bodies the same line.

The heads initially pin the plane of 20 components. Within each pinned
plane the weak edges impose parity equations on the two possible line
assignments. A trial plane for an unpinned component can contradict one
of these equations with two different neighbors. Propagating those
contradictions pins all 3,142 planes. Within them, 32 connected groups
can switch their two lines independently. This gives exactly 2^32
assignments of the paired D3 tails with incumbent bodies and heads.
All 312 previously verified Barnes–Wall head frame choices also work
with each assignment. The contact-count changes come from body products
equal to 32 between groups in the same plane. Adding the attainable
changes independently across the three planes and the head choices
gives the 283 totals.

For a new real tail, each group assigns two upper bounds to the two
absolute line projections in its plane. Switching the group allows
either order. Sort both bounds and take the componentwise minimum over
groups in a plane. Comparing the two sorted tail projections with
this pair is necessary and sufficient. This removes the need to
enumerate the 2^32 configurations. The shell and face programs find
zero feasible tails on all 452,088,000 Leech contact midpoints,
137,159,728 deleted-root contact triangle centroids and 2,647,762,920
deleted-root contact tetrahedron centroids. They use only the switched
midpoint rows as obstacles.

Among the 2,953,318 fresh contacts with both parents already deleted,
a fresh midpoint with any real tail conflicts with at least three
midpoint pairs after the best switch. Exactly 512 contacts attain
three when only midpoint obstacles are counted. For this pricing
calculation, sort each of the 85 minimal integer tail profiles within
each plane and discard dominated profiles. This leaves 28. For each
profile, count violations on each group in both assignments and charge
the cheaper assignment. Newly recomputing possible switches after
removing points is outside this count.

A separate finite exchange search treats every native projection
individually. In one fibre it offers all 16 original roots and both
tail signs on each of the six lines at all 112 contact midpoints.
It can switch the 32 exterior groups and delete individual exterior
midpoint signs or incumbent heads. Points gained must come from the
selected fibre. All 4,320 searches finish without a gain.

Classes 24 and 134 give a larger simultaneous exchange. We released
all 64 of their fibres and discarded all their mutual constraints,
as well as all auxiliary head and pole constraints, when bounding
them. The remaining exterior midpoint constraints permit 128 maximal
phase domains. Their local capacities sum to exactly 1,216 points,
already present in the incumbent. The 64 capacities are 17 on 12 fibres,
18 on eight, 19 on 12 and 20 on 32.

The 96 regular classes were also rebuilt simultaneously, keeping the two
exceptional classes, the other original layers and the original heads.
Each regular class or previously empty class can use a complete square
from any two axes of its native frame, any coordinate tail plane and
either assignment of its two perpendicular lines. The square repeats at
all 32 projections of that class. The exterior permits 3,168 states in
98 classes, including classes 28 and 108. A gain in this family requires
97 occupied classes. For each of the 98 choices of a class to omit, the
solver enforces pairwise geometric compatibility and consistency of the
exterior phase switches. Propagation and exhaustive branching exclude
all 98 choices. The source attains 96 regular square classes.

The local capacity test uses the 11-dimensional residual vectors after
subtracting the common 16-dimensional projection. Distinct retained
residuals have nonpositive inner products. Split their Gram matrix on
the connected components of its negative entries. An irreducible
positive semidefinite block with nonpositive off-diagonal entries has
nullity at most one, by Perron–Frobenius applied to a sufficiently large
scalar multiple of the identity minus the block. A singular block of
size two is an antipodal pair, and every other singular block has at
least three points. A selection of m points in dimension d therefore
contains at least 2m−3d antipodal pairs and has m ≤ 2d. The verifier
enumerates the necessary orthogonal pairs and tests the remaining
compatibility graph by exact clique search.

### Polar replacements

There is also an obstruction to replacing the 12 pure-tail poles with
points whose 24-dimensional body has length s ≤ 4. Each signed
coordinate-axis head layer contains a full orthoplex of body vectors
of length 8 in dimension 24. Thus each absolute tail coordinate is
bounded by a = 8 − s/√24. Its three squared coordinates sum to
128 − s² = 2a² + m², where m² = 32s/√24 − 13s²/12. For s ≤ 4 this
gives m ≥ 2s/3. Concentrating the squared coordinates up to their
cap a² shows that their absolute values sum to at least 2a + m.

For two points in a common closed octant, write the absolute tail
coordinates as a minus three nonnegative deficits with sum at most
a − m. Their tail dot product is at least aa′ + am′ + a′m.
Subtract ss′ for the worst body product. Substituting m ≥ 2s/3 and
ss′ ≤ 2(s+s′) gives a result strictly greater than 64 whenever
s+s′ > 0. At s=s′=0 the inequality is also strict if either point
is interior to the octant, since its coordinate sum then exceeds 16.
A new point lies in one open octant and conflicts with each of its
three incident poles. At most one new point fits per octant and each
removed pole serves at most two octants.

For bodies on any of the 196,560 Leech root rays this obstruction
extends throughout the polar region of squared tail norm at least
64. The six midpoint and three axial-head support values take just
35 profiles. On the interval 1/3 ≤ a ≤ 3/4 for the multiplier of a
Leech root, exact four-dimensional polytope vertex enumeration bounds
the squared total norm by 108. Every admissible polar point on these
rays therefore has body length below √128/3 < 4.

### Two-for-three rotations

In each native fibre fix the common projection w and choose an incumbent
opposite residual pair (w, v), (w, −v). Remove these two points. Allow
new residuals of the form sv + u, where u lies in the three tail
coordinates, u·v = 0 and 64s² + ‖u‖² = 64. These replacements can
move the body away from every Leech contact midpoint. Of three such
residuals, two have coefficients s of the same weak sign. If both have
|s| > 1/√2, their residual dot product is positive by Cauchy–Schwarz.
Consequently any three mutually compatible replacements need at least
one candidate in |s| ≤ 1/√2.

For a fixed sign of s, each retained point with tail t gives a bound
t·u ≤ a − b|s|, where a and b are integer dot products and offsets.
On each subinterval, use the looser of the two endpoint bounds for
each constraint. Enumerate the vertices of the resulting bounded
three-dimensional polytope by integer Cramer's rule, including
lower-dimensional intersections. If every vertex falls short of
the smallest required squared radius in the subinterval, no point
in the subinterval works. Subdividing [0, 182/256] covers the full
required belt. This excludes both signs for all 40,870 opposite pairs
in Takhanov's 4,320 native fibres. The replay checks 97,612 intervals.

### Importing points from 28 dimensions

The repository's 223,528-point 28D construction has four tail coordinates.
Intersecting their space with a hyperplane gives a 27D construction. We
enumerated the 120 hyperplanes spanned by three independent source tails.
Keeping every point of each section and restoring every compatible Leech
root gave a maximum of 205,244 points, attained in four sections.
[28D coordinates and verification](../../verifications/improved/dim28-takhanov-root-exchange/)

For the four coordinate and eight diagonal hyperplanes, an isometry
identifies the section with integer 27D coordinates. For each section and
each of the 48 signed permutations of the three resulting tail axes, we
allowed every point of Takhanov's 27D code, every point of the section,
and every one of the 196,560 equatorial Leech roots. Individual points
from all three sets could be selected freely. Each of the 576 unions
has maximum independent set size exactly 210,108.

To check the upper bound, take the 27D source on one side of a bipartite
graph. On the other side put the imported nonequatorial points and all
13,248 missing roots. Join pairs with dot product exceeding 64. A
matching covers the entire candidate side in 534 alignments. For each
remaining alignment, start alternating paths from unmatched candidates.
Every candidate outside the reached set has a distinct matched blocker
outside that set. Within the reached set
the bipartite graph splits into components, each with twelve candidate
points and eight incumbent blockers. In all 608 components, exhaustive
subset enumeration shows that every compatible subset removes at least
as many incumbent points as it adds.

We also dropped the common section and alignment assumptions. Take any
of the 25,600 imported midpoint bodies or 188 imported head bodies with
any signed native 27D tail of the required norm, along with the pure-tail
poles. Retain every fresh candidate that conflicts with at most twelve
incumbent points. This yields 3,644 imports. Adjoin all 13,248 deleted
Leech roots. Arbitrary simultaneous exchanges from this 16,892-point
pool have maximum gain zero.

The certificate uses 15,850 cliques of pairwise incompatible candidates.
Each clique has a positive integer weight and one source point that
conflicts with every member. The weights of the cliques containing
each candidate sum to at least 24. Those charged to any single source
point sum to at most 24. Selecting k candidates therefore charges at
least 24k to the removed source points, which forces at least k removals.
The verifier checks each clique by integer dot products and can rebuild
the entire candidate catalogue from the two source configurations.
[Clique cover](data/d27_section_clique_cover.npz)

### Reproduction

Use the scratch inputs from `contact_tail_input.py` and `face_input.py`
described above. The switch group labels use the same midpoint ordering.
The fibre solvers require OR-Tools, NumPy, SciPy and Numba.

```sh
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python midpoint_switches.py --groups /tmp/kn27-groups --out /tmp/kn27-switched.npz
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python ../26_27_native_projection_obstructions_20261118/verify_exact.py /tmp/kn27-switched.npz
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python midpoint_switches.py --count-contacts /tmp/kn27-switched.npz
g++ -O3 -std=c++17 -fopenmp contact_tails.cpp -o /tmp/kn27-contacts-check
g++ -O3 -std=c++17 -fopenmp deleted_faces.cpp -o /tmp/kn27-faces-check
g++ -O3 -march=native -std=c++17 -fopenmp contact_prices.cpp -o /tmp/kn27-contact-prices
taskset -c 0-7 env OMP_NUM_THREADS=8 /tmp/kn27-contacts-check /tmp/kn27-contacts /tmp/kn27-groups
taskset -c 0-7 env OMP_NUM_THREADS=8 /tmp/kn27-faces-check /tmp/kn27-faces /tmp/kn27-groups
taskset -c 0-7 env OMP_NUM_THREADS=8 /tmp/kn27-contact-prices /tmp/kn27-faces /tmp/kn27-contacts /tmp/kn27-groups
taskset -c 0-7 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python fibre_trades.py --workers 7 --out /tmp/kn27-fibre-wins
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python special_fibres.py
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python regular_squares.py
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python polar_caps.py /tmp/kn27-contacts
taskset -c 0-7 env OPENBLAS_NUM_THREADS=1 python residual_rotations.py --workers 7
taskset -c 0-7 env OPENBLAS_NUM_THREADS=1 python sections_28d.py --workers 7
taskset -c 0-3 env OPENBLAS_NUM_THREADS=4 python section_cover.py --rebuild-catalogue
```
