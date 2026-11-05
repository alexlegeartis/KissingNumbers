# One projected E7-frame construction in 13D

This is an exact restricted result, **not a bound on the 13D kissing number**.
Take the supplied 1,932-point E7-frame construction in 14D, orthogonally
project every point to the hyperplane perpendicular to (1,-1,-1,0,...,0),
then normalize the nonzero projected vectors. The original 1,122 points in
the hyperplane remain mutually compatible. Conversely, the certificate gives
810 disjoint pairs of incompatible projected points, so no compatible subset
of these 1,932 projected directions has more than 1,122 points. Thus this
particular direct projection cannot reach the 1,154-point 13D baseline, let
alone improve it.

`verify.py` checks the 14D source, the 1,122-point equator, and every conflict
in the matching using exact integer squared comparisons. It applies only to
this source and projection: other directions, new vectors, moving the source,
and different constructions remain open.

    OPENBLAS_NUM_THREADS=1 python research/e7_projection13_20261105/verify.py
