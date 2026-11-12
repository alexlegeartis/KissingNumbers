# Gamma72 continuation: stronger counts from extra normal coordinates

This research directory strengthens existing geometric constructions by counting
their source slices with a larger, explicitly embedded normal root lattice.
The extra normal coordinates constrain the population; they are not removed from
the resulting configuration.

## Exact candidates ready for independent review

| Dimension | Bound at campaign start | New exact candidate | Increase |
| --- | ---: | ---: | ---: |
| 68 | 471,884,288 | **471,886,208** | 1,920 |
| 69 | 743,522,304 | **746,325,051** | 2,802,747 |

The 69-dimensional value includes two additional pure head poles. The source
population bound is the ceiling of `6716925440/9`, namely 746,325,049.
The maintained registries and older proof packages are not edited here.

The published inputs are the integrality and minimum squared norm 8 of Nebe's
Gamma72 lattice, its 6,218,175,600 minimal vectors, and the spherical 11-design
property of that shell. The checker verifies their finite-dimensional application;
it does not enumerate the ambient shell or reprove those published inputs.

## Reproduce the exact checks

From the repository root:

```sh
python research/dim32plus_continue_20261113/gamma/verify_extension.py \
  research/dim32plus_continue_20261113/gamma/d7_compact68.json --enumerate
python research/dim32plus_continue_20261113/gamma/verify_extension.py \
  research/dim32plus_continue_20261113/gamma/d7_compact69.json --enumerate
```

These checks require only Python's standard library. The optional `--enumerate`
flag directly sums every one of the 162,969 possible profiles and independently
cross-checks the compact combinatorial orbit sums and all selected moment rows.
It does not enumerate the 6.2 billion ambient minimal vectors.

The numerical search can be reproduced with NumPy, SciPy and SymPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  research/dim32plus_continue_20261113/gamma/extension_model.py 5 6 7
```

That command derives the candidate duals from LP active columns and then checks
their slacks in exact rational arithmetic. The exact checker does not run an LP.

## Actual normal embedding

The first five normal vectors and the catalogue Gram matrix come from the pinned
input at:

`verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json`

Its SHA256 is
`5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842`.
In the catalogue's zero-based integer basis the first seven normals are:

```text
u0 = e0
u1 = e3
u2 = e4
u3 = e16
u4 = -e0 - e16 + e17
u5 = e1 + e2 - e5 - e18
u6 = -e2 + e4 - e6 - e7 + e9 + e16
```

`d7_embedding.json` records the complete integer vectors. Direct multiplication by
the 72-dimensional Gram matrix verifies that their Gram matrix is that of `2D7`:
diagonal 8, entries -4 at edges `(0,1)`, `(0,2)`, `(0,3)`, `(3,4)`, `(4,5)`,
`(5,6)`, and zero at all other off-diagonal positions.

In orthonormal coordinates the corresponding root-basis rows start with

```text
(1,-1, 0, 0)
(0, 1,-1, 0)
(0, 1, 1, 0)
(-1,0, 0, 1)
```

and append successive coordinate differences. Multiplying these rows by 2 gives
the actual normal vectors up to an isometry. The basis determinant has absolute
value 2, so it generates the complete D7 root lattice.

## Why the larger model counts the old geometry

Write the orthonormal normal projection as `h/4`. Normal duality makes every
coordinate of `h` integral, with all coordinates of the same parity. Every minimal
vector outside the 84 known normal roots satisfies

\[
 |h_i|+|h_j|\le8\quad(i\ne j),\qquad |h|^2\le128.
\]

These restrictions give 162,969 possible profiles and 60 signed-permutation
orbits. The known roots have two entries equal to `+/-8` and all others zero.

For a fixed D4 profile `a`, let `O(a)` be its orbit under signed permutations and
the half-Hadamard triality of D4. If `S` is the selected set of old profiles, put

\[
 f(a)=\frac{|O(a)\cap S|}{|O(a)|}.
\]

For a D7 orbit `G`, the objective coefficient is exactly

\[
 c_G=\sum_{h\in G}f(h_1,h_2,h_3,h_4).
\]

First average the selection over internal D4 orientations, and then average the
D7 orientation over signed permutations. Every orientation is an actual integer
change of normal root basis. No symmetry of the ambient lattice is assumed.
The internal D4 transformations need not extend to D7: they already give valid
four-normal tuples inside the first four orthonormal coordinates.

The 84 known D7 roots contribute zero to both selected objectives, checked
explicitly. Their first four coordinates are zero, a single `+/-8`, or a D4 root;
none of these is in a selected orbit.

Only the first four normal coordinates enter the geometric projection. The
remaining three D7 coordinates stay inside the 68-dimensional body. Thus
marginalizing those three coordinates counts exactly the source bodies used by
the existing 68-dimensional projection and 69-dimensional head compression.

## Exact moments and duality

If the even powers are `2*a_i` and `s=sum(a_i)<=5`, the design identity is

\[
 \sum_{v\in\operatorname{Min}(\Gamma_{72})}\prod_i h_i(v)^{2a_i}
 =N\,128^s\,
 \frac{\prod_i(2a_i-1)!!}{72\cdot74\cdots(72+2s-2)}.
\]

The known normal-root contributions are subtracted. Signed-permutation averaging
leaves 18 independent moment rows. Each row is reduced to primitive integers;
the exact rows, right sides, objective, rational dual, and all slacks are recorded
in the compact certificate JSON files.

For nonnegative averaged per-profile populations `x`, the equations are `A*x=b`.
The verifier checks `A.T*y<=c` in exact rational arithmetic. Therefore the averaged
selected population is at least `b.T*y`. Some actual orientation has an integer
population at least its ceiling.

| Normal system | Exact 68D population bound | Exact 69D population mean |
| --- | ---: | ---: |
| D4, previous package | 471,884,288 | 743,522,304 |
| D5 | 471,884,928 | 6,706,829,312 / 9 |
| D6 | 471,885,568 | 2,238,440,960 / 3 |
| D7 | 471,886,208 | 6,716,925,440 / 9 |

The earlier explicit-cell D6 LP had poor numerical scaling and reported a value
below its exact dual. The compact orbit-population scaling fixes this: numerical
and rational values agree. None of the earlier inconsistent floating-point values
is used as a proof.

## Geometry and the two additional poles

For a selected D4 profile `h`, put `q=|h|^2/16` and let `b` be the source point's
orthogonal component in the 68-dimensional complement. The 69D map is

\[
 v\longmapsto(\alpha_h b,t_h),\qquad
 \alpha_h^2=\frac{8-t_h^2}{8-q}.
\]

The ten rational heights are unchanged from the published head-compression
package. The checker recomputes all 55 diagonal/cross profile cases exactly using

\[
 (4-t_it_j)^2-
 \alpha_i^2\alpha_j^2(4-h_i\cdot h_j/16)^2\ge0,
\]

with the required signs checked before squaring. The 68D construction uses its
six inherited profiles and zero heights. Every output has squared norm 8 and all
distinct pair products are at most 4.

For 69D, each height has absolute value at most 1. The two pure head points
`(0,+sqrt(8))` and `(0,-sqrt(8))` therefore have product at most 4 with every
existing point. Every old image has a nonzero body, so the poles are distinct.

## Bounded exploration, not certificates

- An actual A4 normal subsystem was tested with a finite scalar-head grid; its
  best numerical population was about 703.2 million, below the 69D incumbent.
- D5, D6 and D7 all reproduce exactly 1,360,986,112 for the incumbent three-slice
  A2 construction in dimension 70. More counting information alone gives no
  improvement for that selection.
- Four-term and then six-term sparse normal searches found D6 and D7. The direct
  six-term D8 search and differences of the finite D7 candidate pool found no D8
  extension. These are bounded search outcomes, not nonexistence theorems.
- `search_a2_head71.py` tested 4,009 rational-scaled head vertices on 61 A2
  profiles, with three bounded population-weighted rounds. All three exactly
  reproduced the incumbent 2,603,658,750 in dimension 71.
- `alternate_pair_search.py` tested actual two-normal Grams with scalar products
  0, 1, 2, and 3. It checked all 81, 85, 67, and 75 maximal compatible affine
  slice cliques, respectively. The best certified candidate counts were
  `3571023872/3`, `2518815561/2`, `1322796635`, and `1291273230`, all below the
  70D incumbent. Finite scalar-head grids of 6,657, 6,383, 5,803, and 4,989
  vertices recovered the 71D incumbent; the orthogonal case requires its two
  normal-span poles. These are finite search results, not optimality proofs for
  all possible deformations.
- An actual extra orthogonal normal `e1-e10-e19` gives D7 x A1. Mixed moments
  yield 296 unknown orbits and 43 rows. The full projection bound is imposed.
  The only norm-8 boundary orbit has absolute D7 profile `(0,4,4,4,4,4,4)` and
  last scalar product `+/-4`; exact catalogue-coordinate membership rejects all
  896 candidate projections. The exact dual bounds nevertheless reproduce D7 alone
  for dimensions 68, 69, and 70.
- A further finite D8 search tried sums of 3,594 verified minimal vectors, each
  supported on at most four signed catalogue basis coordinates. None of the 60
  profile-compatible root pairs produced the requested D8 extension.

Generated vertex pools, heuristic outputs, temporary Gram exports, compiled
binaries, and numerical LP diagnostics are research data. The compact exact
certificates, explicit normal vectors, and the standard-library verifier are the
reviewable proof route. The maintained standalone package is
`verifications/improved/dim68-69-gamma72-extended-moments/`; it includes a second
checker written independently and package rejection controls.
