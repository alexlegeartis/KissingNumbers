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
