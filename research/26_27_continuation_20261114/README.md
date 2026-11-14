# 26D and 27D continuation checkpoint, 14 November 2026

**No improved kissing configuration has been obtained.** The comparison targets
are PackingStar's **201,180 in 26D** and **209,862 in 27D**. This research
checkpoint does not register a new lower bound.

The checkpoint preserves the following completed searches, their inputs, and
code for continuing them. It follows
the [earlier finite-pool certificates](../26_27_exchange_exclusions_20261114/).
Search receipts describe finite scopes; unsuccessful searches do not imply
global optimality.

| Search | Completed scope | Result |
|---|---|---|
| Six-support Golay midpoint family | 5,440,512 point tests in 26D; 16,321,536 in 27D; at most eight source blockers per point | 14,264 and 18,464 new candidates |
| Flat midpoint family, 26D | 4,521,984 point tests; cutoff eight | 2,944 new candidates |
| One absent and one retained Leech parent, 26D | All 9,775,320 bodies; all root-tail colors | No states with at most two blocked source orbits; 144 with three |
| Same parent enumeration, 27D | All 26,998,768 bodies; all six root-tail colors | 546 states with at most two blocked source orbits |
| Squares between octads meeting in four coordinates | 13,601,280 cycles per dimension | 942,728 cycles admitted in 26D and 1,309,200 in 27D with the fixed head and normal layers |
| Direct cycle replacement and released-parent restoration | The admitted cycle catalogues above; all available orthogonal root-tail pairs | No improvement |
| Connected cycle cavities | 1,524 scopes in 26D, freeing at most four incumbent cycles; 64 scopes in 27D, freeing at most five | No larger compatible clique; no recursion timeouts |
| Mixed-field individual-point case, 27D | 8,592 candidates, 2,774 incident source points | CP-SAT returned **2,774**, with matching objective bound |

The last row concerns the literal saved finite union. Its independently checked
graph contains 214,378 genuine conflict edges. The exact CP replay is supplied;
no standalone dual certificate is claimed for this new case. The earlier
package retains its independent rational dual certificates.

The direct restoration scans were resumed across overlapping index ranges.
The reports and logs preserve those ranges: 26D starts at zero and resumes at
655,360; 27D starts at zero and resumes at 122,880 and 266,240. The 16-point
cycle tests and complete-cycle cavity tests do not exclude arbitrary partial
point exchanges.

## Verify and replay

Use Python 3.12. `requirements.txt` records the discovery environment. Reserve
at least one available logical CPU and limit nested numerical threads. On the
five-CPU execution environment, CPU 4 was reserved.

From this directory:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python verify_checkpoint.py --controls

taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python point_mis_cp.py data/point_case27 \
  --out /tmp/kn27-point-replay --max-cost 4 --seconds 1200 --seed 367

taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python cycle_cavity_clique.py data/cycle26 \
  data/priced_six26.npz data/priced_intersection26.npz \
  --out /tmp/kn26-cycle-replay --max-free 4 --seconds 1200

taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python one_missing_midpoints.py 26 \
  --out /tmp/kn26-parent-replay --budget 3
```

`verify_checkpoint.py` hashes the archived files, checks the Leech seed norms,
reconstructs both complete source point sets exactly from the cycle models,
recalculates the saved one-missing-parent orbit costs, checks every norm and
conflict edge of the reduced 27D graph, and checks all 182,106 pairs in the
published 604-point module. Radical comparisons in this checker use enclosing
rational intervals with integer endpoints, independently of the discovery
code's sign routines. Corruption controls reject a changed norm, an incorrect
orbit cost, wrong parents, a duplicated module point, and a false conflict edge.

This checker does **not** check all pairs of either PackingStar source and
does not establish a new lower bound. Candidate generation coverage is recorded
by the enumeration programs and their receipts, rather than inferred from
checking the saved candidates alone.

The large intersection-cycle catalogues are regenerated, not bundled:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python intersection4_cycles.py --dim 26 --out /tmp/kn26-squares
```

Use `--dim 27` for the other dimension. `octad33_pool.py` regenerates the
six-support and flat families with `--family six` or `--family flat`.
`data/cycle26/` and `data/cycle27/` contain the model inputs needed for the
saved cavity and reconstruction checks; large heuristic-search graph binaries
are omitted.

## Continuing searches

The joint midpoint/head models allow removal of head orbits, partial parent
graphs, and, in 26D, the four irrational orbits. They maximize the actual
population formula and add violated midpoint-pair constraints after each
proposed solution. Current search logs are snapshots of running experiments.
**A relaxed model score above the source score is not a kissing configuration.**
No such proposal has passed the complete pair constraints at this checkpoint.

The bounded 26D continuation can be started with:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python midpoint_head_cp.py 26 data/joint26.npz \
  data/missing_integer26.npz data/six26.npz data/aa26.npz \
  data/module500_pool26.npz data/flat26.npz data/one_missing26.npz \
  --out /tmp/kn26-joint-heads --seconds 1800 --seed 359
```

Its default limits are 32 deleted original midpoint orbits and two deleted
head orbits. The archived running 27D experiment used no such limits; reproduce
that scope with `--max-old-deletions 0 --head-deletions -1`, inputs
`joint27.npz`, `missing_integer27.npz`, `six27.npz`, and `hadamard4_27.npz`,
and seed 347. Time-limited solver trajectories may differ across machines.

The next individual-point continuation should explicitly include restoration
of signed Leech parents when their last blocking midpoint points are deleted.
This permits asymmetric exchanges that the complete-orbit models miss.

## Sources and provenance

The starting repository commit is `8ca0d210ebb391d906e423f4940871b17d91def6`.
PackingStar inputs are pinned to
[`2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6`](https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111).
Their unchanged bytes are bundled, with the same SHA-256 values as the preceding
certificate package. Parent graph indices refer to the bundled Leech ordering.

The 604-point 11D module is due to Together's EinsteinArena team, from
[`c388c6f7408c886311940896713339a1a70c2394`](https://github.com/togethercomputer/EinsteinArena-new-SOTA/tree/c388c6f7408c886311940896713339a1a70c2394/kissing-number),
file `solutions/solution_n=604_d=11.json`, SHA-256
`0bde9ca2c434d7beb5af63a64a291498e8c264d9a5718b738e92dd70b9ba7761`.
The 500-point 10D input is from
[Sloane's coordinate table](https://neilsloane.com/packings/dim10/kiss.10.500.txt).
These are published inputs, not results of this campaign. Their rotated
candidate pools and bounded selection receipts are retained for continuation.

`manifest.json` binds the archived code, data, and receipts. Historical paths
inside original reports identify the discovery workspace; replay scripts use
only this package's data and the repository's Golay implementation.
