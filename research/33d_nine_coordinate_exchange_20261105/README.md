# A nine-coordinate exclusion and 33D neighborhoods

This is a restricted search result, **not a kissing-number upper bound**.
Consider antipodal lines represented by vectors with exactly one zero and
all remaining eight entries ±1 in a fixed set of nine coordinates. Suppose
all distinct signed points have pairwise inner product at most 4 in this
unscaled representation. There are at most **64 antipodal lines**.

To prove it, complete the zero in both possible ways. The two completions
of a point form an edge of the 9-cube. A line and its antipode produce two
edges. These edges, over all lines, form an induced matching: for two
points with different zero coordinates the seven common signs must disagree
in at least two positions, and for the same zero coordinate the eight common
signs must disagree in at least two. The antipodal absolute-product condition
covers all choices of signs.

An induced matching in an n-cube, n≥2, has at most 2^(n−2) edges. For m
matching edges, its 2m endpoints have 2m(n−1) edges to outside vertices.
An outside vertex cannot have all n neighbors among these endpoints: if it
did, the partner of any one neighbor would also be adjacent to another,
contradicting the induced-matching property. It therefore has at most n−1
endpoint neighbors. There must be at least 2m outside vertices, so 2^n≥4m.
For n=9, m≤128; each antipodal line contributes two edges, proving the bound.

Starting with the repository's 33D certificate, enumerate every support
obtained from a full octad by replacing one coordinate, and count its
intersections of size at least five with the full fibers. There are 369,600
distinct non-incumbent supports. Of the 131 supports with exactly one such
incumbent fiber, their 94 owner groups have the following sizes and union
sizes, where the union includes the owner:

- 69 groups: one candidate, union size 9.
- 15 groups: two candidates, union size 9.
- 3 groups: two candidates, union size 10.
- 2 groups: three candidates, union size 9.
- 5 groups: four candidates, union size 9.

The nine-coordinate bound excludes a positive exchange using only the owner
fiber and its candidates in the 91 groups with union size 9, even if both
parities and arbitrary subsets of their individual sign lines are allowed.
Deleting exceptional points cannot improve this count. The other three groups share six coordinates across the owner and both
candidates. They are also excluded by a second bound: any family of signed
octads whose supports share a fixed six-set has at most 64 antipodal lines.
Orient every line so one fixed common coordinate is positive. There are 32
sign patterns on the remaining five common coordinates. For each pattern,
lines on different octad supports conflict, because their signed overlap is
at least 5. On any one support, at most two lines with that common pattern
are compatible, because their two remaining signs must both differ. This
gives at most 32·2=64 lines.

Thus no positive exchange exists in any of the 94 groups when only its one
owner fiber and these single-coordinate replacement supports are changed.
Candidates with multiple owner fibers, other support replacements, and
completely different constructions remain open.

`replay_neighborhoods.py` independently replays the base certificate and the
complete enumeration using exact integer intersections. Both inequalities
are also included as valid optimization cuts in `search.signed_fiber_exchanges`.

## Matching bound for joint-fiber neighborhoods

For distinct octads meeting in 5, 6 or 7 coordinates, the 64 incumbent
antipodal lines of either parity and all 128 possible lines on the other
support can be partitioned into 64 mutually conflicting triples. Match the
signs on the common coordinates, match all but one of the outside signs,
and use the last outside sign to distinguish the two new parities. Therefore
the number of selected lines on that candidate support is at most the number
of deleted lines on the incumbent fiber. This bound permits arbitrary partial
fibers and either incumbent parity.

More generally, if every candidate support can be matched injectively to a
distinct conflicting flexible incumbent fiber, summing these inequalities
proves that the exchange cannot have positive net gain, even before charging
for deletions from exceptional lines. The engine uses this exact matching
rule only for the explicitly supplied finite pool. If no such matching exists,
the search remains open and proceeds to the pointwise optimization.
