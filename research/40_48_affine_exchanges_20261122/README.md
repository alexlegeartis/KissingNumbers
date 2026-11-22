# P48 affine searches in dimensions 40–48

The [40D certificate](../../verifications/improved/dim40-affine-centre-refill/)
gives 1,126,327 points. The search enumerated the complete 1,092,000-point
rank-40 minimum shell and considered all 2,240 nonzero isotropic centre
classes. For each class it retained the complete affine core and both
short layers before checking minimum-shell insertions.

| Compatible minimum-shell points | Number of centres |
|---:|---:|
| 0 | 376 |
| 1 | 624 |
| 2 | 336 |
| 3 | 464 |
| 4 | 168 |
| 5 | 144 |
| 6 | 32 |
| 8 | 48 |
| 10 | 48 |

At the preceding centre, the conflict graph between the 1,126,314-point
affine construction and the complete minimum shell had 13,157,568
edges. A matching of size 1,091,992 bounds every exchange in that union
by 1,126,322 points, attained by the eight free additions. At centre 267,
a matching of size 1,091,990 gives the corresponding bound of 1,126,324
for its union. The ten-point refill attains it.

The individual-point search then priced coefficient norm-36 candidates
against the centre-267 construction. Three additional bodies from the
same discriminant class have no blockers. The exact check uses
norm-30 differences to screen their contacts with the affine core.
With the ten minimum-shell points fixed, the integer trade over
11,787 norm-36 candidates costing at most five deletions attained
the same three-point gain.

The whole-family trades in 41D, 42D, and 44D–46D fixed all body families
of radius index at least 14 and proposed radius-12 families on rational
grids in the retained head spans. In each finite model with at most
six displaced radius-12 families per proposal, the best gain was zero.
Short-difference tests also priced individual points when a universal
angle bound allowed only the first possible difference shell to
violate an angle.

For 41D, 42D, and 44D–46D, a further 1,000 Weyl relabellings per
dimension tested minimum-shell insertions at head zero. Complete
radius-12 and radius-14 body shells decided the candidate contacts.
The 43D test used the explicit tiny layer and medium split insertions.
It kept the full parent medium families in the compatibility filter.
These searches gave zero additional points. In 45D, a discrete
trade along the sixteen pure-head directions priced 252 distinct
candidate bodies. Each required a different existing point to be
removed, in addition to a pure head.

The 47D searches grouped the P48 minimum shell by its inner products
with short normal vectors. Whenever a profile difference was absent
from the minimal shell, the cross-profile product bound dropped from
three to two. Rational replacement heads recovered the 23,766,960-point
minimal-vector section.

## Capacity with the 48D zero-class affine core

Consider the rank-40 body construction with the zero-class norm-six
core and an eight-dimensional head space. Use the P48 classwise
gluing bounds and complete body families of radius index
`r ∈ {10, 12, 14, 16, 18}`. The family equations and populations are
from the [P48n head certificate](../../verifications/improved/dim41-46-p48-heads/).

In the scaled coordinates of that certificate, let `z` be a head
and let `μ(c)` be the shortest squared normal representative of
body class `c`. The diagonal gluing test requires
`r + ||z||² ≤ 18`. Its test against the core then requires
`r + μ(c) ≤ 18`.

For a fixed body point in a family repeated at `m` heads, each unit
vector has projection at least `sqrt(r/18)` onto the body direction.
The squared norm of their sum is at least `m²r/18`. Their kissing
angles bound the same squared norm by `m(m+1)/2`. Thus
`m(2r−18) ≤ 18`.

| Radius index | Eligible body classes | Copies per class | Total capacity |
|---:|---:|---:|---:|
| 18 | zero class | 1 | 1,092,000 |
| 16 | 240 root classes | 1 | 27,896,400 |
| 14 | 2,160 classes | 1 | 19,828,800 |
| 12 | 2,240 classes | 3 | 3,326,400 |
| 10 | 1,920 classes with 30,240 bodies in total | 9 | 272,160 |
| Pure heads | eight dimensions | | 240 |
| **Total** | | | **52,416,000** |

The P48 minimum shell attains the capacity. The argument applies to
complete families certified by these gluing inequalities with this
fixed core and radius catalogue.

## Complete families with a free core in 40D

The complete-family model assigns a weight to every body class and radius
index. It permits two families when the classwise upper bound on their
doubled inner product satisfies the kissing inequality. In 40D the heads
vanish, so the compatibility test is $r+s-18+\mu(c-d)\leq\sqrt{rs}$.
The exact maximum in this model is **1,126,314**, with or without a
requirement to start from an incumbent family.

| Families permitted | Ceiling |
|---|---:|
| A radius-18 family, with any centre | 1,126,314 |
| At least one radius-16 family and no radius-18 family | 1,126,313 |
| Radius indices 10, 12, and 14 | 963,270 |

A radius-18 family can coexist only with radius-10 and radius-12
families whose labels differ from its centre by zero or a root class.
Every nonzero isotropic centre attains the first ceiling. For the second
case, an exact weighted clique search covers 2,160 branches, one per
possible first radius-16 family. Each branch uses proper colourings to
bound the remaining weight. For the third case, a four-dimensional
totally isotropic subspace partitions the discriminant group into 81
cosets. A radius-14 family occupies its coset exclusively, displacing 27
radius-12 families. Separately, 1,117 disjoint incompatible pairs bound
how many radius-12 families can survive. Both searches may charge the
entire radius-10 population of 30,240 to the upper bound.

Run `python research/40_48_affine_exchanges_20261122/whole_family_bound.py`
from the repository root. It uses the existing discriminant data and
recomputes the bounds with integer arithmetic.

## Joint point exchanges in 40D

For each of the 2,240 nonzero isotropic centres, keep the complete core
and its two complete short layers. Enumerate all norm-6 lattice vectors
and all other coefficient-norm-36 bodies that fit those fixed layers.
Select the two refill types jointly, with exact integer cross products.
The best refill has 13 points and is attained at 48 centres.

The centre-267 exchange model allows proposals to displace points from
any incumbent layer. Short difference shells enumerate the blockers.
The test covered the available radius-10 and radius-12 catalogue, 1,431
radius-14 classes, 42 radius-16 classes, and 25 radius-18 classes,
including the zero-class minimum shell. The pool of 211,940 distinct
proposals with at most eight blockers admits a matching that assigns
every proposal a different incumbent blocker. Consequently any trade
in that pool loses at least as many points as it inserts, even before
enforcing compatibility between the proposals.

For 24 alternate nonzero centres, the two configurations remain
independent within their own sides. Their union therefore poses a
bipartite independent-set problem, including exchanges with no bound on
the number of removed points. A maximum matching covers every
nonshared point of the alternate configuration in each of the 24 unions.

Outside these residue families, an exact search tested 19,112,706
reflected core proposals. A pole is a shortest dual vector, and each
core vector with raw product at least 21 with the pole is reflected
through its perpendicular hyperplane. Every proposal had at least two
incumbent blockers. A separate test of 2,391,775 rational rays around
one- and two-blocker proposals found no insertion compatible with the
entire incumbent. These ray counts include repetitions across seeds.

## Binary-code route

Ten cosets of a rank-16 binary length-40 kernel, with mutual Hamming
distance at least ten, would give 655,360 dense sign points. Adding
128 even signings on each of the 3,755 octads and the 3,120 signed
two-coordinate roots would give **1,139,120**. All dense–sparse
interfaces follow from the weight-eight support bound.
[Constant-weight code data](https://aeb.win.tue.nl/codes/Andw.html).

An exact syndrome distance table and a clique search gave a maximum
of nine cosets in each of 72 tested shortenings of the extended
length-48 quadratic-residue code. Splitting the kernels by two or
three marked coordinates gave 100 finite tests. Splitting by one
additional parity equation gave 72 more. In these split tests each
proposal with at most four blockers conflicted with exactly one
incumbent fibre, and the proposals sharing that blocker were a
clique in the conflict graph. The result rules out a gain in these
exchange pools.

## Changing the rank-40 lattice

The published \(P_{48p}\), \(P_{48q}\), and \(P_{48n}\) complements
all have 1,092,000 minimal vectors. Their discriminant classes have
different distributions of shortest dual vectors. Exact enumeration
in each complement gives 30,240 dual vectors of squared norm \(10/3\).
Those vectors occupy 1,920 classes in \(P_{48p}\) and \(P_{48n}\), and
1,904 in \(P_{48q}\). The distribution is invariant under lattice
isometries, so the \(P_{48q}\) rank-40 complement is a different lattice
from the other two complements. [Published section data](https://github.com/wcgbg/sphere-packing/tree/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim40-delta41).

The 16 additional empty classes have shortest normal representative of
squared norm \(8/3\). They permit norm-36 candidates to fit the core
automatically when their difference class is one of these 16. The
rank-40 minimum-shell sweep over 2,240 centres in each complement
found maxima of eight available points for \(P_{48p}\), nine for
\(P_{48q}\) and ten for \(P_{48n}\). For \(P_{48q}\), joint selection
over the strongest centres initially combined 30 norm-36 points with
seven norm-six points. Checking the classes with three short
differences then added three more points to that configuration.

Using the stronger difference bound for the 16 empty classes leaves
the complete-short-family maximum at 41,724 points beside the core
in each of the 2,240 \(P_{48q}\) centres. A pointwise search over the
empty difference classes in all 2,240 centres attains a best net
refill of 37 when each proposal can displace up to five short-layer
points.

Including difference shells with three or six vectors makes it
profitable to displace short-layer points jointly. Sixteen centres
admit a net gain of 51 over their complete affine construction. The
published witness adds 93 norm-36 points and removes 30 norm-36
and 12 norm-30 points. Exact searches over all difference classes
whose shortest normal representative has squared norm \(8/3\)
give a maximum net gain of 51 across the 96 centres tested. Adding the
difference classes of normal squared norm two in 18 tested centres,
or allowing eight blockers in six centres, gives the same maximum.
[1,126,365-point certificate](../../verifications/improved/dim40-p48q-affine-refill/).

## Distance-ten binary cosets in 40D

Shortening the extended binary quadratic-residue [48,24,12] code on any eight coordinates gives a dimension-16 kernel in length 40. Nine translates with union minimum distance ten always exist. Use the syndromes of zero and of the eight deleted unit columns. The remaining columns span the syndrome space. Any body representing a single deleted syndrome extends to a nonzero codeword with one added bit, and a difference of two such syndromes extends with two added bits. The respective body distances are at least eleven and ten.

An exhaustive computation gives **exactly nine** for every such kernel, allowing arbitrary translates. The translation and inversion generators on the projective line over the field of 47 elements preserve the code. Orbit–stabilizer counts cover all 377,348,994 eight-subsets with 7,497 representatives. For each representative, start with Hamming weights on 24 independent syndrome columns and update all 16,777,216 distances once for each of the other 16 columns. An independent breadth-first search checks the entire distance table for the first representative.

Translate one prospective coset to zero. Join two syndromes of leader weight at least ten when their difference also has leader weight at least ten. No representative admits a nine-clique in this graph. The largest vertex count is 463. This excludes the proposed 40D construction using ten full translates of these kernels. Other length-40 kernels and partial translates remain open.

The same upper bound of **589,824 dense points** holds if any three of the nine standard translates are kept whole and all other binary words can be chosen individually. This covers 629,748 choices of shortening representative and retained triple. Any eligible syndrome must differ from all three fixed labels by syndromes of distance at least ten. The zero translate and the forty one-bit translates partition into radius-one balls about the kernel words, contributing at most 65,536 points in total. For every other eligible nonstandard syndrome, the verifier finds a distinct nonzero standard target at distance less than ten. Adding a fixed shortest error gives a perfect matching between that candidate translate and its target. Their combined contribution is at most 65,536 points per target. The remaining unpaired standard translates have the same individual bound. Thus even arbitrary partial replacements require changing at least seven of the nine standard translates to improve the dense layer.

Two intact translates of opposite Hamming parity suffice for the same bound, for arbitrary translate labels. Translate one to zero and denote the odd label of the other by `y`. The other eligible odd labels form a set `T`, and the eligible even labels are exactly `y + T`. Match pairs of opposite-parity candidate cosets at distance below ten. As above each matched pair contributes at most one coset's worth of points. Across all 60,121 choices of odd difference and shortening representative, the two intact translates, matched pairs and unmatched translates contribute at most nine coset populations in total.

Replay on an x86-64 host with AVX2, BMI2, GNU libstdc++ and OpenMP, with assertions enabled. The optional final argument to the executable limits the number of representatives for a shorter run.

```sh
g++ -O3 -std=c++17 -mavx2 -mbmi2 -mpopcnt -fopenmp research/40_48_affine_exchanges_20261122/qr40_all_shortenings.cpp -o /tmp/qr40-all
taskset -c 0-126 env OMP_NUM_THREADS=16 /tmp/qr40-all
```

## Unequal coordinate weights and lifts from 39D

Shortening on nine coordinates gives a length-39, dimension-15 kernel. Put a real height on the fortieth coordinate and normalize each resulting fibre. Taking the original two levels at heights −1 and +1 gives eighteen fibres. Twenty would suffice to beat the current 40D bound with the complete sparse layer. The 21 quarter-integer heights from −2.5 to +2.5 all preserve the sparse layer. Even for an octad through the weighted coordinate, the largest possible unnormalized dot product is at most seven plus the absolute height, and its normalized value is at most one half.

All 32,345 shortening orbits and all nine choices of the marked coordinate were checked. Among the resulting 291,105 marked cases, none allows twenty fibres after removing up to three of the original eighteen and choosing replacements at any of the 21 heights, anywhere in the full syndrome space. For each height one of the two old levels requires Hamming distance at least ten. Its surviving fibres bound the candidate list by translates of the deep syndromes. The program exhausts the compatible extensions for every removal set.

A separate calculation lifts two copies of the 39D construction with 36 sign coordinates and a four-dimensional head. Each copy has eight distinct central heads of length at most two, ten body cosets per central head, and four exceptional fibres. Use the same 4,096-word kernel and allow any relative translation in the full quotient. Of the 16,777,216 translations, all but 64 produce a body difference of Hamming weight at most five that forces at least sixteen disjoint conflicts. At distance six, compatibility requires antipodal heads of length exactly two. For the remaining 64 translations, counting these conflicts by head pair forces at least 34 deletions. This bounds the two-copy model by **152 fibres, or 622,592 dense points**. Its total with the 3,755-octad layer and all 3,120 roots is at most 1,106,352. The argument permits arbitrary orientations of the two head systems. [39D source construction](../../verifications/improved/dim39-qr-three-coordinate-heads/).

```sh
g++ -O3 -std=c++17 -mavx2 -mbmi2 -mpopcnt -fopenmp research/40_48_affine_exchanges_20261122/qr40_one_coordinate.cpp -o /tmp/qr40-one-coordinate
taskset -c 0-126 env OMP_NUM_THREADS=16 /tmp/qr40-one-coordinate
python research/40_48_affine_exchanges_20261122/qr40_double_lift.py
```

## Norm-eight sections in 47D and 48D

Take a norm-eight vector in an extremal even unimodular lattice of dimension 48. Index projected classes by the inner product with the normal modulo eight. Choosing its representative between −3 and 4 shows that projected squared radius at most six requires only the norm-six and norm-eight shells (apart from the origin). The degree-ten spherical moments give the class populations. The norm-eight parent shell contributes a class of 20,489,368 projected vectors of squared radius six. [Spherical designs of extremal lattices](https://www.mathnet.ru/php/archive.phtml?jrnid=mmj&option_lang=eng&paperid=245&wshow=paper).

Use the minimum squared distance of each difference class to certify pairs of whole families. In the 47D model with all families at zero head, the maximum accepted by this test is 20,876,648. In 48D, allow a one-dimensional head, repeated copies of families at different heights and two possible pure head directions. Constrain the total squared norm before normalization to at most six. The same difference-class test reaches exactly **52,416,000**, attained by the original minimum shell.

To see the 48D limit, a fixed nonzero body of squared radius greater than 4.5 lies at latitude angle less than 30 degrees, so that family can occur only once. The radius-four family admits at most two copies, since two copies on the same side of the equator would be separated by less than 60 degrees. The two radius-six families have zero height and their difference class has squared distance four, so the test excludes their joint selection. Summing the remaining populations, even with both pure head directions, falls short unless the original radius-six core and both radius-four copies are present. With that core, the class-distance test forces each radius-four copy to have total squared norm six. Their heights are then opposite and each excludes the pure head in its direction. These restrictions apply only to complete classes selected by this distance test.

```sh
python research/40_48_affine_exchanges_20261122/norm_eight_section.py
```

The complementary search in ranks 41–43 produced three explicit improvements. [Configurations and exact replays](../../verifications/improved/dim41-43-affine-sections/).

## The determinant-128 P48m section in 45D

The published P48m section with normal Gram matrix `[[6,2,2],[2,6,-2],[2,-2,6]]` has 128 discriminant classes. Its three affine families of squared radius 23/4 contain 7,379,838 points at the origin. An additional family can be normalized separately and contributes 12,746 points at squared radius 9/2. The resulting **7,392,584** is the exact maximum for complete projected families of squared radius at most six, with zero heads and compatibility decided by the minimum squared distance of the difference class. The repository's 45D configuration has 7,655,056 points. [Source certificate](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim45-kissing7379838/P48m-d914e7d8.cert.txt).

To recover the populations, enumerate the integer triples of inner products with the normal basis. Within each discriminant class the least squared normal length is at most 5/2. An ambient representative for any eligible projected family can therefore be chosen on shell six or eight. The spherical 11-design equations relate the populations in the 68 classes up to sign. The source certificate gives the extra equation needed on shell six. For a harmonic polynomial of degree twelve, the weighted theta series has weight 36 and begins in degree three, so it is a multiple of the modular discriminant cubed. The degree-four coefficient is −72 times the degree-three coefficient. One such equation in normal direction (1,2,4) completes the system on shell eight. The replay solves both systems over the rationals, checks every moment, and searches the resulting graph of 149 families exhaustively. An independent enumeration in the published 45D body Gram matrix returned the predicted populations of 49 and 941,720 at normal inner products (−3,−3,−1).

```sh
python research/40_48_affine_exchanges_20261122/p48m_45_section.py
```

## Choosing a scaled D4 in 44D

For seventeen projected families in the complement of a scaled D4, the population is 4,940,961 minus seventy times the sum of three small shell counts. In the previous 4,939,281-point model the sum was 24. Enumerating the 36,848 minimal vectors at inner product three with a fixed minimal vector allowed a search over mutually orthogonal triples within that set. P48m, P48p and P48q all supplied triples for which the three counts vanish. The P48m coordinates give **4,940,961** and require just 1,656 listed vectors to certify the calibrating slice. [Exact construction and population proof](../../verifications/improved/dim44-d4-section/).

Nine additional projected families of squared radius 13/3 pass the class-distance test against sixteen of the seventeen selected families. Each contains 4,488 points and conflicts with one retained family of squared radius five. A matching assigns every one of the **40,392** proposed points a distinct incumbent point at an incompatible angle. Consequently the entire union, even with unrestricted point deletions, has maximum 4,940,961. [Exact matching replay](d4_refill_matching.py) and [integer coordinates](d4_refill_matching.npz).

One-, two- and three-dimensional replacement heads on the new section were also searched on finite grids. The searches found no improvement in 45D–47D within their branch budgets.
