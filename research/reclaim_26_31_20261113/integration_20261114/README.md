# PR #3 integration with the verified November 14 codes

PR #3 contributes the component-wise 26D exchange optimizer, the full E8/D3
27D search and rational stereographic sphere tails. Main contributes the
verified **201,188-point 26D** and **209,884-point 27D** configurations.
The searches now consume those actual certificates; no further record is
claimed by this integration.

The import-breaking literal newline in `solve_joint26.py` was repaired.
`search_sources.py` checks metadata and exact norms and converts the current
certificates to integral norm-32 coordinates. The 26D `--rebuild-blockers`
option recomputes the old pool's interactions with the new source. Solver
budgets start after the 26D replay, and HiGHS uses one worker. The 27D engines
remove unnecessary deleted profiles and retain the empty exchange whenever
a timeout incumbent would lower the count. Nine tests pass, including both
pinned and current source formats and rejection of fractional coordinates.

`results.json` retains exact source hashes, run parameters and the solver
reports. All commands used `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and
CPU affinity restricted to CPUs 0–3, leaving CPU 127 free as required.
The commands in `../SEARCH_26_27.md` reproduce the search setup.

The 26D run reprices all 2,004 candidates from the 759-octad pool against
201,188 points. Its selected four-for-four exchange has zero gain and passed
all changed-pair checks. The 27D runs test D3 and rational sphere tails
against all 209,884 points, with full projected-profile multiplicities.
Limits and timeouts are recorded explicitly. `--seconds` is the HiGHS
optimization limit; preprocessing and some presolve work can make total
elapsed time longer. These runs do not exclude
other pools, longer searches or different constructions.

A useful next model should jointly reselect the successful midpoint/root
orbits and the new tail candidates, pricing shared exterior deletions.
The present 27D engines hold retained midpoint coordinates fixed unless
they are deleted, so they cannot reroute them to accommodate a new head.
The verified 444-point normal sector and its frame-33 seed are preserved in
`../../../verifications/improved/dim26-27-radical-orbits/` for that purpose.

| Run | Pool | Outcome |
| --- | ---: | --- |
| 26D repriced octad pool | 2,004 | Four additions/four removals; gain 0 |
| 27D D3, blocker cost at most 64 | 368 | Solver returned gain 0 |
| 27D rational tails, K=2, initial cap | 4,000 | Timed out; no positive proposal |
| 27D D3, blocker cost at most 256 | 12,000 | Timed out without an incumbent |
| 27D rational tails, K=2, larger cap | 8,000 | Timed out without an incumbent |

The initial rational-tail incumbent deleted unnecessary profiles, giving a
negative objective. The integrated engines trim such deletions and preserve
the empty exchange as a known feasible baseline. The initial raw report is
retained for transparency; it was never published as a coordinate improvement.
