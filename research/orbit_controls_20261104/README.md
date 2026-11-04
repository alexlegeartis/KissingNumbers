# Transferable orbit proofs and missing1–24D controls

**No new SOTA was found.** This continuation adds general construction/proof
instruments, fills five exact target-control gaps, and searches central sections
across four other dimensions. It follows commit118193b without replacing its
archive, receipts or negative research results.

New shared capabilities are numerical exterior-square representation lifting,
exact rational orthogonal orbit construction and compressed certification,
inherited exact central sections, explicit stabilizer-aware action screening,
and weighted quadratic-coordinate input. There are35 registered operators.
See `discovery/docs/EXACT_ORBITS.md` for exact scopes and proof obligations.

## Known targets now reconstructed through the shared engine

| Dimension | Count | Constructor and independent proof |
|---|---:|---|
|17|5730|Boolean quadratic/RM formulas generate the odd BW16 equator, two512-point tiers, disjoint192-point tiers and poles. All16413585 weighted integer pairs and rank17 check exactly.|
|18|8358|Existing bent-hexagon component codes are reconstructed as weighted Q(sqrt3) coordinates. The native checker verifies all34923903 pairs and full rank18.|
|22|49896|Searched central section of the freshly certified Leech source; complete subset membership and exact rank22.|
|23|93150|Searched central section of the same source; complete subset membership and exact rank23.|
|24|196560|Actual orbit generation from four small algebraic actions, then independent closure, transitivity, seed separation and exact rank24.|

These controls use the same registered DAG. The17D/18D programs supply explicit
coordinate coefficients to native input operators; their formula/code generation
is documented by `build17.py` and `build18.py`. They are not held-out discoveries.
The Leech program supplies a seed and four exact actions, not a record point
array. It uses the cyclic Golay permutation, a bounded deterministically seeded
Golay-automorphism search, a Golay sign change and an exact sextet map.

The Leech construction/proof took about20 seconds. Its compressed proof checks
786240 exact generator images and196559 seed comparisons, covering all
19317818520 unordered source pairs. No floating-point distance tolerance is
used. Full point-array materialization remains explicit and budgeted.

## Actual section construction search

`section_search.py` tests all signed coordinate-pair normals at every stage,
requiring a new independent central constraint and choosing the largest section.
This catalogue does not contain target counts. Exact source/subset replay then
certifies the selected point sets.

| Dimension | Normal candidates tested | Best output | Captured target |
|---|---:|---:|---:|
|23|552|93150|93150|
|22|551|49896|49896|
|21|549|27720|31560|
|20|546|17400|19448|

The last two are clear failures. No whole-family optimality or record claim is
made. The canonical24-row matrix now has exact target recipes in
1–10,13,14,16–18,22–24; gaps remain11,12,15,19–21. The earlier200-program held-out
benchmark remains a different measurement, matching ten targets with formula
controls and four through actual support search. New control replay is not
retroactively credited as blind discovery.

## Character/stabilizer audit and bounded28D work

The exact class-average principle was checked independently, including the real
irreducible/complex-type caveat. A full orbit with class average `a` and class
size `m` needs at least `ceil(m*(2*a-1))` fixing elements. Equality at one half
forces tight equations only with no fixing elements. Partial orbits are outside
this prune.

The E8 reflection average and root-outer-product spanning were verified by exact
matrices. Independently generated whole-orbit controls reproduce7560 rank-two,
37800 rank-four and15120 complex-structure points in28D. The rank-two example
has30 fixing root reflections and nonscalar `A^T*A`; it confirms that the
stabilizer exception is essential. The three-block seed has a directly checked
reflection violation. No reported60-clique partial-orbit bound was used without
an available independent witness. These controls are all below221988, and no
further28D campaign was run after the instruction to broaden dimension coverage.

## Resource profiles and reproduction

The coverage matrix explicitly assigns20 million pair checks to17D,40 million
to18D, and a250000-point/orbit cap plus one billion transform/rank-work units to
the Leech source/sections. These are declared replay profiles, not general
dimension-based construction rules. Large self-contained runner receipts are
stored as `.json.gz`; the reader bounds their decompressed size.

```sh
OPENBLAS_NUM_THREADS=1 python research/orbit_controls_20261104/build17.py
OPENBLAS_NUM_THREADS=1 python research/orbit_controls_20261104/build18.py
OPENBLAS_NUM_THREADS=1 python research/orbit_controls_20261104/build_leech.py
OPENBLAS_NUM_THREADS=1 python research/orbit_controls_20261104/section_search.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --check
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s discovery/tests -t .
OPENBLAS_NUM_THREADS=1 python -O -m unittest discovery.tests.test_group_constraints \
  discovery.tests.test_exact_orbits discovery.tests.test_certified_sections \
  discovery.tests.test_weighted_quadratic
```

The broader objective remains a genuinely stronger certified construction.
No installation, subagents, Actions modification or PR is part of this work.
