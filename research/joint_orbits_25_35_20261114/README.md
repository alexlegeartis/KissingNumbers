# Joint exchanges across 25–31D and a 35D extension

14 November 2026. Combining older certificates with the latest native design
produces **248,212 points in 30D**, 32 above 248,180 and 220 above pinned
PackingStar. The [portable complete proof](../../verifications/improved/dim30-certificate-unions/)
checks all 30,804,474,366 unordered final pairs. All record coordinates are
rational; the net exchange is 52 additions and 20 removals.

`union_certificates.py` optimizes a union of two verified codes by exact
bipartite matching. Source deletions can be shared by multiple additions.
The two successful 30D steps restore 24 Barnes–Wall points without deletion,
then add 28 joint-refill points while removing 20. A second complete donor
pass finds no further gain. This is optimality within each pairwise union,
not a global kissing-number upper bound. Cached blocker lists are reused
only while the source configuration and algebraic field remain compatible.

`joint_midpoint_tails27.py` implements the requested joint model: midpoint
orbits, restorable roots and rational spherical-tail heads are selected
simultaneously. It prices each blocked incumbent profile once with its true
multiplicity. Four-sign orbit products use exact rational dot comparisons.
The largest run has 7,878 candidates and 4,159,960 conflicts, requires a
fractional tail, permits shared exterior deletions, and minimizes the deleted
population subject to positive gain. It reaches its 300-second solver limit
without a feasible improvement. The smaller cutoff-16 model is infeasible;
the cutoff-32 model also times out. Those finite results do not exclude larger
catalogues or different symmetry choices.

`axis_frames.py` jointly changes midpoint/root/radical orbits while replacing
one complete E8 head frame. It reconstructs the current certificate exactly
before searching. `model26.npz` retains the 26D orbit catalogue; the 27D input
is the published `model_orbits27.npz` in the earlier proof package. Both 26D
axes return the existing 388-point sector in their 275-frame families. The
27D axis-1 run finds a different 444-point sector but no larger one within
480 seconds. Other axes and simultaneous frame changes remain open.

`union_rational25.py` checks older rational general points and restorable
Leech points against the current 25D certificate using rigorous coordinate
intervals and exact signs in Q(sqrt(2),sqrt(3)). Its 111-point donor pool
has a complete matching, so it cannot improve the current source by itself.
`sections.py` generates the coordinate and pair-coordinate hyperplane
sections retained in `sections.json`, using the campaign input certificates
in 26–31D (the 30D input was the preceding 248,180-point code); none improve the neighboring dimension.
The certificate-union scans also find no gains in 26–29D or 31D.

`partial_stars35.py` extends the search beyond 31D: it selects individual
signed star points and individual root deletions, retaining the current
middle layer and six head lines. Its 992-candidate finite model is infeasible
above the existing 415,334-point count. This does not exclude changed heads,
changed middle supports or other point shapes.

`results.json` retains reports and `logs/` retains solver transcripts.
`manifest.json` pins the campaign files; these research receipts do not
replace the independent 30D geometric proof. `test_searches.py` checks a real
shared-deletion gain, a zero-gain complete matching, cached graph consistency,
algebraic signs and the orbit-product identity.

Run from the repository root, with NumPy, SciPy and OR-Tools installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 python -m unittest discover -s research/joint_orbits_25_35_20261114 -p test_searches.py -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 python research/joint_orbits_25_35_20261114/joint_midpoint_tails27.py --out /tmp/joint27 --cost 64 --body-cost 32 --require-fractional --local --seconds 300
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 python research/joint_orbits_25_35_20261114/union_certificates.py --dimension 30 --out /tmp/union30
```

Adjust affinity for the machine and reserve at least one logical CPU. All
campaign jobs used one worker on CPUs 0–3, leaving CPU 4 free on this host.
The 35D catalogue additionally requires a C++17 compiler. Numerical solver
bounds and timeouts are not independent exact exclusion certificates.
