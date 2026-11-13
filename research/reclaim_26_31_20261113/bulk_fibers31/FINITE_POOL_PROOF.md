# Exact no-gain certificate for the 3,960-point residual-axis pool

This certificate closes only the supplied candidate pool against PackingStar's
pinned 251,254-point source. It gives no upper bound for the unrestricted
31-dimensional kissing problem or a larger candidate catalogue.

All coordinates have squared norm128 and maximum allowed pair product64.
Every candidate has zero coordinates on the normal octad
`I={0,2,4,5,6,10,11,23}`. Its complementary sixteen-coordinate core has
squared norm64, as does its seven-coordinate tail.

There are **3,928 candidates on twelve axis tails** and **32 further
candidates on two signed-tetrad tails**. Thus the pool is not entirely an
axis-tail pool, although a direct full product scan found that distinct-tail
candidate groups are mutually compatible. The proof below does not need
that cross-group compatibility observation.

For each of the twelve axis tails, the source already contains32 points with
the same tail and zero normal coordinates. Their cores are a complete16D
orthoplex. Any feasible final selection at that fixed tail can contain at
most32 such cores: the fixed tail product is64, so different core products
must be nonpositive. Consequently, across these twelve disjoint groups,

```math
\text{selected axis candidates}\leq
\text{deleted points among the384 old axis-fiber points}.
```

For completeness, the dimension-d orthoplex bound is elementary. Normalize
the core Gram matrix G to diagonal1. Its off-diagonal entries lie in[-1,0],
so `sum_ij(G_ij^2+G_ij)<=2m`. On the other hand,
`sum_ij G_ij^2 >= m^2/d` by the frame-matrix trace inequality, and
`sum_ij G_ij >= 0` is the squared sum of the core vectors. Hence `m<=2d=32`.

The stored certificate assigns each of the32 remaining candidates a
different incompatible source point **outside** those384 old axis-fiber
points. The checker verifies the injection and each forbidden product
directly with integers. Therefore

```math
\text{selected nonaxis candidates}\leq
\text{deleted points outside the384 old axis-fiber points}.
```

Adding the two disjoint charges proves that inserted points never outnumber
deleted points. Every feasible subset of this source union this candidate
pool therefore has at most251,254 points. This is an exact finite-family
bound, not an optimizer status or a timeout conclusion.

The supplied source and pool hashes are checked before verification. The
matching need not be trusted: `--verify-only` independently checks the saved
assignment, complete exception coverage, distinct blocker identities,
separation of the victim sets, all relevant source norms and fiber geometry,
and every matched forbidden angle.

```bash
taskset -c 6-7 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python research/reclaim_26_31_20261113/bulk_fibers31/finite_pool_select.py \
  ../PackingStar/25D-31D/new_bounds_20261111/D31_251254_coordinates.npz \
  --verify-only
```

Adapt CPU affinity to the host and leave one logical CPU free. Omit
`--verify-only` to reconstruct the matching certificate. A further MILP
attempt is unnecessary after this exact zero-gain upper bound.
