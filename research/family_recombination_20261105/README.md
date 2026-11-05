# Whole-family GF(32) recombination in 32–35D

This searches substantially different complete octad-support families, not
individual point additions. Three nonlinear coordinate permutations of GF(32)
are applied to the entire existing weight-eight support code; extra coordinates
in dimensions 33–35 are fixed. The engine combines all generated families
and solves the set-packing problem using one inequality per five-subset.
Every new support is required to intersect each protected exceptional support
in at most four coordinates, so the existing exceptional lines may be retained.

    OPENBLAS_NUM_THREADS=1 python research/family_recombination_20261105/run.py --dimension 32 --seconds 120

Each additional valid support adds 128 kissing points to the supplied signed
construction. The program independently checks the source and, if a numerical
proposal increases the support count, constructs a complete certificate and
runs the separate signed-octad verifier before reporting a point count. A
solver timeout or infeasibility status is not an exact exclusion or a record.
Comparison against current published and repository records remains necessary.
Only three specified maps and a finite union of families are searched.

For each map considered alone, `verify_two_family_matching.py` checks a
matching from every retained transformed support to a different conflicting
source support. Thus none of those three two-family unions can exceed the
1,704-support source, even after arbitrary switching. This exact obstruction
is specific to those pairs. It does not rule out recombination across three or
more families, other maps, partial fibers, or changes to protected layers.
