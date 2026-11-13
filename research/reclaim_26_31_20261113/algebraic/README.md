# A larger 31D module and two geometric alternatives

13 November 2026. The successful construction in this directory gives
**251,395 points in dimension 31**, **141 above** the pinned PackingStar
251,254-point configuration. It is a transplantation of a known
2,564-point 15-dimensional construction, reassembled so that its
240-point E8 equator and 84 D7 tail roots agree with the incumbent.
No improvement to the 15-dimensional kissing number is claimed.

The final configuration passes a complete exact replay of
**31,599,597,315 unordered pairs**. The portable verification package is
[`dim31-leech-sloane-module`](../../../verifications/improved/dim31-leech-sloane-module/).
The results below distinguish the successful finite construction from
the two unsuccessful geometric pilots.

## 1. What changes in the successful construction

Use squared norm 32 and product limit 16. The eight core coordinates are
`I = {0,2,4,5,6,10,11,23}`, and the seven tail coordinates are 24 through 30.
Their 15-dimensional coordinate subspace contains 2,282 incumbent points:
240 equatorial E8 points, 1,958 mixed heads, and 84 pure D7 tail roots.

The preceding Barnes–Wall replacement had only 2,340 points in this
subspace. Here the known Leech–Sloane 2,564-point recipe is reassembled
to give exactly the same E8 equator and D7 tails, together with **2,240 heads**.
Every head has core squared norm 16 and tail squared norm 16.

The first placement kept the whole new module and deleted all 178 external
blockers. It gave

```math
251254 - 2282 - 178 + 2564 = 251358.
```

The subsequent bounded search tested 3,530 E8/D7 orientations and solved
the exact bipartite union problem for each. The best union keeps 2,559
of the 2,564 module points, deleting only 136 external points:

```math
251254 - 2282 - 136 + 2559 = 251395.
```

After cancelling common points, the actual final delta is **2,219 additions
and 2,078 deletions**. The deletions consist of 1,942 old module points and
136 points outside the module. Exactly 340 old points already belong to
the oriented module and remain in the final configuration. Five otherwise
new module points are omitted.

### The explicit 2,564-point recipe

Start with 15 coordinates in the quadratic form `diag(1^14,2)`, using
squared norm 16 and product limit 8. Pair the first fourteen coordinates
into seven pairs. The four point classes are:

| Class | Count |
|---|---:|
| Single entries ±4 in the first fourteen positions | 28 |
| All sixteen sign patterns ±2 on each of 91 four-supports | 1,456 |
| Two affine binary length-14, size-512 glue codes, at final coordinate −1 and +1 | 1,024 |
| On each coordinate pair, all choices ±2, ±2, ±2 including the final coordinate | 56 |
| **Total** | **2,564** |

Twenty-one of the four-supports are unions of two coordinate pairs.
The remaining seventy choose one position from each of four different
pairs. They form an exact cover of the 280 triples using three different
pairs. The 91 supports therefore cover all 364 coordinate triples exactly
once, and distinct supports intersect in at most two positions.

The two glue codes are specified without an external data dependency in
[`ls15_recipe.json`](ls15_recipe.json): their origins are 16,383 and 12,149,
with the common nine binary generators
`680, 2210, 8330, 5160, 1290, 360, 240, 60, 15`.
The recipe records the source data hash and the coordinate pairing used
to obtain this presentation. Generation checks these affine code sets
against the original source rather than assuming their identity.

For `i = 0,...,6`, make the Euclidean coordinate change

```math
z_i=y_{2i}+y_{2i+1},\qquad
z_{8+i}=y_{2i}-y_{2i+1},\qquad z_7=2y_{14}.
```

This doubles the quadratic form, giving squared norm 32 in Euclidean R15.
The first eight Euclidean coordinates now contain the E8 equator.
The resulting populations by tail squared norm are exactly
`{0: 240, 16: 2240, 32: 84}`.

Every distinct pair of the 2,564-point module is checked directly with
integer arithmetic. The construction is not justified by an exact-cover
solver status alone: the selected seventy blocks, both actual glue codes,
the coordinate change, norms, uniqueness, and all geometric inequalities
are rebuilt and checked.

### Orientation and union proof

The final E8/D7 isometry has an integer matrix R and denominator 2, with
`R R.T = 4 I`. The convention is
`oriented = reconstructed @ R.T / 2`.
The full oriented module, R, the admitted module indices, and the
minimal delta are preserved in `ls15_orientations/`.

Each source and each oriented module is already a valid configuration.
Consequently their union has conflicts only between old-only and new-only
points. Maximum bipartite matching and the alternating-path vertex cover
give an exact optimal subset of each particular two-copy union. This does
not prove global optimality over all orientations or all configurations.
Every accepted positive changed-pair set is also checked directly against
the complete source. The final complete verifier removes any dependency
on the original source's claimed geometric validity.

The initial 251,358 witness is retained in `ls15_transplant/` with its
direct int64 changed-pair proof and complete geometry report. The final
251,395 witness and complete geometry report are in `ls15_orientations/`.
Intermediate full coordinate files are omitted; their orientation and
search records remain in the report.

## 2. A genuine deformation of the large body fibers

Write a source point as `(a,b,t)`, with `a` on I, `b` on the other sixteen
core coordinates, and `t` on the seven tails. When `||b||^2 = 16`, the
residual vector `(a,t)` has squared norm 16. Two points in the same fiber
must have nonpositive residual product.

There are 4,004 source fibers with 28 points. Direct checking shows that
every one consists of fourteen antipodal orthogonal pairs in residual
R15. Of these, 3,988 have a vacant coordinate tail axis. The other sixteen
have a vacant four-support tail direction; they are outside this pilot.

Choose any two old orthogonal residual axes p and q and a vacant pure-tail
axis r, all with squared norm 16. Replace the four parents `b ± p, b ± q`
by the six points `b ± v1, b ± v2, b ± v3`, where

```math
v_1=(8p+q-4r)/9,\quad
v_2=(p+8q+4r)/9,\quad
v_3=(4p-4q+7r)/9.
```

The coefficient matrix divided by nine is orthogonal. Hence all six
points have squared norm 32, and their pair products are at most 16.
The unchanged points of the same fiber are orthogonal to p, q, and r.

More generally, for positive integers m,n with `4mn <= 2n^2+m^2`, use

```math
\frac{1}{2n^2+m^2}
\begin{pmatrix}
2n^2&m^2&2mn\\
m^2&2n^2&-2mn\\
-2mn&2mn&2n^2-m^2
\end{pmatrix}.
```

After projection perpendicular to r, each new point lies in the convex
hull of the four old parents. Therefore every unchanged old point with
zero r-coordinate remains compatible. In particular the entire remaining
Leech equator remains compatible. This argument allows p and q to have
both core and tail coordinates; they need not be pure-core axes.

The finite pilot tested three types of axis pairs per fiber, ten rational
parameters, and four sign choices: **478,560 exact states** across all
3,988 coordinate-axis fibers. The best six-point bundle had **90 external
blockers**, in addition to its four parents, giving net gain −88.
No positive bundle was found. The figures do not exclude shared deletions
between bundles, joint changes of other fibers, different parameters, or
an unrestricted change of an entire residual frame.

The precise larger problem is to assign full orthogonal 15-frames to the
norm-16 body fibers while controlling products between different fibers
and with the retained exterior. For a fixed body b and an outside point
`y=(b_y,r_y)`, each selected residual axis v must obey
`|v·r_y| <= 16-b·b_y`. A full 15-frame necessarily satisfies
`16||r_y||^2 <= 15(16-b·b_y)^2` for every retained outside point, by Parseval.
This is a feasibility condition, not an achieved construction. No large
continuous frame optimization was launched in this campaign.

## 3. An irrational cube family and its exact capacity bound

Another tested grammar uses all 2,160 E8 squared-norm-16 core vectors h,
paired with tails `(4/sqrt(7))*s`, where `s` ranges over the 128 signed
seven-dimensional cube vertices. These 276,480 candidates have squared
norm 32 and clear all but 2,136 of the incumbent points by exact finite
core/tail inequalities.

The complete exact pricing found **at least 21 incumbent blockers per
candidate**. Comparisons use integer squares and integer square roots;
no floating tolerance decides the sign of an expression involving sqrt(7).

There is also a stronger structural reason not to pursue a complete
replacement by this one cube alphabet. Each fixed tail permits at most
sixteen residual core points, since their pair products are nonpositive
in R8. Equality forces an orthogonal cross-polytope. If a neighboring cube
vertex, differing in one sign, carried another core vector h', its products
with all sixteen points of that cross-polytope would have magnitude at
most `32/7`. For its eight orthogonal axes q_i of length four, Parseval gives

```math
\sum_{i=1}^8(h'\cdot q_i)^2=16\|h'\|^2=256,
```

whereas the proposed inequalities would bound this sum by
`8(32/7)^2 = 8192/49 < 256`. Thus a tail occupied by sixteen points forces
all seven neighboring cube vertices to be empty.

Every cube edge consequently has total endpoint occupancy at most thirty.
Summing over the 448 edges of the seven-cube yields
`7N <= 448*30`, or **N <= 1,920**. This holds for arbitrary real norm-16
core vectors, not merely the 2,160-point E8 catalogue. It excludes replacing
all 1,958 incumbent module heads by a larger set from this cube grammar.
It does not exclude a union of some cube points with retained old heads.

## Reproduction and source pinning

The source is PackingStar commit
`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`, file
`25D-31D/new_bounds_20261111/D31_251254_coordinates.npz`, SHA-256
`714c9a0189cd1ff2f9363f6672e84506606ca54cb806a18d155d817481019943`.

From the repository root, with a neighboring PackingStar checkout:

```bash
taskset -c 2-3 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python \
  research/reclaim_26_31_20261113/algebraic/leech_sloane_module.py \
  --source verifications/improved/dim73-95-gamma72-caps/data/pole_src_15.npy \
  --out /tmp/ls15-module --seconds 25 --exact

taskset -c 2-3 env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python \
  research/reclaim_26_31_20261113/algebraic/transplant_ls15.py \
  ../PackingStar/25D-31D/new_bounds_20261111/D31_251254_coordinates.npz \
  research/reclaim_26_31_20261113/algebraic/ls15_exact/module.npz \
  --out /tmp/ls15-transplant
```

Adapt affinity to the host and leave a logical CPU free. The exact-cover
discovery may return a different valid seventy-block design. For exact
reproduction of the registered module, use the saved recipe and the
portable verifier, which does not run an optimizer. The orientation search
uses a wall-clock budget; its saved best witness is independently verifiable
without reproducing its timing.
