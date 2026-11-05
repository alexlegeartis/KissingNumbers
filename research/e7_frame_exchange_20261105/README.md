# 14D E7 frame reassignment and a mixed-shell obstruction

This is an exact result in a **restricted search neighborhood**, not a new lower bound or a kissing-number upper bound.

The seed selects one of 135 orthogonal E7 root bases separately for each of 63 first-coordinate root lines. The selections satisfy the cross-basis constraints. Their product vectors, together with two D7 root shells, give 1,932 points in dimension 14. This independently checks the kissing inequalities but does not establish whether this seed is non-isometric to an earlier configuration.

`verify.py` rebuilds that seed and enumerates every integer norm-8 line: 182 lines of type (2,2), 160,160 of type (2,1,1,1,1), and 384,384 of type (1,1,1,1,1,1,1,1). There are 672 non-seed lines conflicting with exactly four seed lines and 7,056 conflicting with exactly eight. There are no other non-seed lines with fewer than 14 seed-line conflicts.

Even if antipodal symmetry is broken, allow arbitrary deletion of seed points and arbitrary addition of points from those two candidate classes. Every four-conflict candidate point conflicts with four seed points, and an affected seed point can block at most four mutually compatible such candidates. The analogous numbers for the eight-conflict class are eight and eight. The sets of seed points affected by the two classes are disjoint. Counting candidate–removed-point incidences separately for the two classes shows that the number of added points cannot exceed the number removed. The verifier checks both local capacity claims directly by exact bitset enumeration, including every affected seed line and both signs by symmetry. Thus this whole mixed 4/8-conflict neighborhood is exhausted, with pointwise as well as multi-layer changes permitted.

This leaves candidates with 14 or more seed conflicts, non-integer shells, changes outside this seed, and entirely different constructions open. Solver results are not used as proof.

Run `python research/e7_frame_exchange_20261105/verify.py` with NumPy and SciPy. The two small binary files are the seed basis and assignment, both regenerated into full coordinates and independently checked at runtime.

## Pointwise search through 14-conflict candidates

The next shell contains 1,664 lines with 14 seed-line conflicts. The following exports all 9,392 non-seed lines with at most 14 conflicts, then searches with both signs and arbitrary point-level deletions. Exact capacity cuts for the exhausted 4/8 classes reduce the model, but the optimization remains a numerical heuristic and may time out.

    python research/e7_frame_exchange_20261105/verify.py --export
    E7_SEARCH_SECONDS=300 python research/e7_frame_exchange_20261105/search.py

A positive incumbent is independently checked in exact integer arithmetic and saved under `cache/`. Search failure or a solver timeout is not a theorem. The capacity proof above is complete only for the stated 4/8 classes, independently of this search.

The verifier also checks every 14-conflict candidate with either sign. Each has 11 forced deletions at which compatible eight-conflict candidates cannot attain their previous local capacity; any two compatible 14-conflict points share at most two of these witnesses. For h = 1, 2, 3, 4, 5 such points, inclusion–exclusion gives at least 11, 20, 27, 32, 35 witness incidences, respectively. Each missing cheap point accounts for at most eight such incidences, so the cheap class loses at least 2, 3, 4, 4, 5 points. The graph connecting compatible pairs with two shared witnesses is also triangle-free. Mantel’s bound gives at most nine such pairs among any six candidates, so their witness union has at least 66 − 15 − 9 = 42 points and forces a cheap deficit of at least six. Hence an improvement in the 4/8/14 shell requires at least seven 14-conflict candidate points. For seven or more, the cheap deficit is at least six. The pointwise search incorporates both exact restrictions.
