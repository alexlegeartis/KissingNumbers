# Orbit limits and phase-block admission

`groups.search` is numerical discovery, even when the generators describe a
finite group mathematically. Rounded orbit keys, eigenspaces, phase grids and
the finite conflict graph do not prove an exact spherical lower bound.

An orbit reaching `linecap` or `vectorcap` has status `limit_reached`. The
symmetry-screening backend does not use that partial orbit. Its representative
check needs the complete vector orbit: moving a pair back to a representative
can move its other endpoint outside a truncated list. A partial list can be
searched as an ordinary point pool only with independent full within/cross
checks; this operator does not silently substitute that different backend.
Reaching a cap does not prove the full orbit infinite or geometrically invalid.

The 28D coordinate-seed reproduction uses the even `p=3, rank=3` Weil action.
At projective caps 1600 and 4000, coordinate seed 0 closes at 1120 lines;
the other 13 coordinate seeds reach their cap. The closed vector orbit has
2240 vectors; an admitted three-phase lift has 6720 points. This describes
that restricted orbit pool, not a dimension-28 ceiling.

A larger follow-up closes coordinate seed 1 at **14560 projective lines** with
cap 20000. Its first 500 actual vector images already contain a product about
`0.6666666666666669` against the seed, above the allowed half-angle threshold.
This is a numerical rejection witness for retaining that entire vector orbit,
even before closure. A witnessed violating pair can rule out a whole-orbit
candidate; an absence of violations in a partial list cannot certify one.
Thus merely completing this larger orbit does not supply a valid large block.
No exact algebraic rejection certificate is claimed by this floating-point
diagnostic, and the other capped seeds are not silently covered by it.

The original phase admission could charge identical blocks repeatedly. An
antipodal vector orbit makes orders 1/2 identical, often also 3/6; invariant
global phases introduce further repeated blocks. The operator now fingerprints
complete numerical point sets independently of row order at seven decimal
digits and admits each such block once. It tries each admissible phase order
at a phase step before moving to the next step. Bounded phase construction
deduplicates incrementally, so a raw `order*orbit_size` upper estimate does
not reject a smaller actual union. These are numerical pool heuristics, not
exact identities or certificate promotion.

`data.block_coverage` reports admitted unique blocks, their charged point rows,
duplicate blocks skipped, numerical identity precision, and scheduling order.
Each inventory item records admitted, duplicate, and capacity-skipped blocks.
Orbit discovery order is unchanged. This is fairness among phase orders within
one closed orbit; it does not implement global fair seed scheduling, orbit
checkpoint/resume, or a persistent search planner.

With the 28D coordinate-only pool, 100000 materialized rows and a 30000000-pair
verification budget, the before/after reproduction keeps the same 6720-point
numerical result. Admission changes from 30 partly duplicated blocks and
`materialization_limit_reached=true` to 16 unique blocks, 53760 charged rows,
32 skipped duplicates, and no materialization cutoff. Thus the cap-waste bug is
fixed, while competitive cardinality remains an unresolved generation problem.

The repository's 221988-point 28D configuration uses a large Leech core,
auxiliary points, coordinated point deletions/insertions, and compatible
nonorthogonal isometric copies. Selecting whole closed Weil phase orbits cannot
express those point-level coupled replacements. Raising an orbit cap alone does
not add that construction grammar. Performance-directed experiments must retain
a competitive source and verify the changed and inherited geometry separately.
The bounded full-sextet continuation in
`research/discovery_campaign_20261104/competitive28` explores such a different
pool; it is research evidence, not a registered large-configuration backend.
