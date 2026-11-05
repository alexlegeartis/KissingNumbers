# 13-dimensional mixed-norm exchange experiment

This is a finite, reproducible test of pointwise exchanges, not a claimed record.
The supplied 1,122-point source is a hyperplane section of an E7-frame-based
14-dimensional construction. The 6,224 candidate directions come from projections
of 24 alternate E7-frame assemblies. They have different squared norms, so the
search uses exact squared comparisons, not rounded normalized dot products.
The normal and source/candidates are included in `pool.npz`; the construction of
that input pool is not reproduced here. The script checks the hyperplane condition
and the search operator verifies the source. An exact replay checks any proposed
exchange. The numerical MILP solver's status is not an exact optimality proof.

Run from this repository root:

    OPENBLAS_NUM_THREADS=1 python research/e7_projective13_20261105/run.py --seconds 120 --gain 33

A gain of 33 would reach 1,155, one above the 1,154 baseline used for this
experiment. It would still require comparison with the then-current literature
and independent verification before any record claim. A timeout or no proposal
is not evidence that no such exchange exists. The data contain a fixed finite
pool, not all possible 13-dimensional constructions.
