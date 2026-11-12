# Bounded 32–39D search and independent 69D review — 13 November 2026

**No new lower bound in dimensions 32–39 resulted from these experiments.**
The existing registered values remain unchanged. This directory also contains
the third independent mathematical/arithmetic review of the successful 69D
construction, described separately below.

## Large support replacements in 34D and 35D

The ordinary Edel–Rains–Sloane construction combines 131,072 dense signs,
128 sign patterns on each weight-eight support, and 2n(n−1) roots. Any family
of octads with pairwise intersection at most four gives such a configuration.
Improving the support packing therefore changes entire 128-point layers.
The experiments begin with the published support witnesses retained in
`verifications/improved/dim32-35-signed-layers/source/public/`.

| Experiment | Exact search pool | Bounded search outcome |
|---|---:|---|
| 35D candidates requiring at most six original support deletions | 187,254 supports, screened from all 23,535,820 octads | 300 seconds; best remained 2,202 supports |
| 35D replacement of every support containing coordinate 0 | 4,687 candidates replacing 251 original supports | 90-second heuristic retained 251; 60-second MILP timed out without a feasible 252-support witness |
| 35D replacement of every support containing coordinate 1 | 6,653 candidates replacing 279 original supports | 120-second heuristic retained 279 after more than 272,000 steps |
| 34D replacement of every support containing coordinate 0 | 3,388 candidates replacing 220 original supports | Bounded heuristic stopped after more than 264,000 steps; best remained 220 |

These outcomes are **not infeasibility or optimality proofs**. The complete
coordinate-star candidate pool can contain supports outside that coordinate:
it includes every octad whose original blockers all belong to the deleted
star. Every retained original support is protected. This allows hundreds of
support layers to change simultaneously.

`star_pool.cpp` enumerates all octads and uses the 5-subset owner map to test
this exact candidate condition. It generated all 35 coordinate pools in 35D
and all 34 coordinate pools in 34D; their population logs are retained.
`support_lns.cpp` uses exact 5-subset incidences for feasible search states.
Its final version charges **distinct deleted supports**, rather than counting
the same conflict once for every shared 5-subset. The first full-pool and
coordinate-0 35D runs used the earlier incidence-priced proposal heuristic;
their feasible-state acceptance was still exact. Heuristic random trajectories
are not part of a proof.

The twelve-blocker pilot generated 6,413,012 candidates but would use about
2.9 GB for incidence data, so it was stopped and the smaller pool was used.
No conclusion follows from that interrupted pilot. Large candidate pools,
temporary binaries, and redundant selected-set files are excluded by the
local publication ignore rules and are regenerated from the sources.

### Reproduction

Compile with a C++17 compiler; OpenMP is used only for full-pool generation:

```sh
g++ -O3 -fopenmp -std=c++17 support_lns.cpp -o support_lns
g++ -O3 -std=c++17 star_pool.cpp -o star_pool
OMP_NUM_THREADS=2 taskset -c 2-3 ./support_lns 35 seed35.txt run35 6 300
taskset -c 2 ./star_pool 35 seed35.txt star35
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 3 python star_milp.py 0 --seconds 60
```

Use CPU affinity appropriate to the host and leave one available CPU free.
The scripts do not update the result registry. A candidate support code needs
an independent exact check of the entire spherical construction before any
new kissing-number bound can be registered.

## Complete sparse-head screen in 38D

`polar38.py` enumerates all 96,096 antipodal sign lines with six nonzero
coordinates, each equal to ±1, in the fourteen-dimensional tail. It checks
the sufficient inequalities `|y·w|≤4` against every original cap head and
every retained rotated flat axis. Exactly **2,016 head lines on 63 supports**
survive, each support carrying all 32 antipodal sign lines.

These are exactly the 49 Fano-grid supports and fourteen six-subsets of the
two seven-coordinate halves already covered by the existing model. Thus
arbitrary signed variants of this specific head shape produce no extra
candidate under those retained layers and sufficient cross-layer inequalities.
The statement does not cover different amplitudes, different norm splits,
modified cap layers, or head/body pairs whose actual products allow a weaker
sufficient head inequality. `polar38_report.json` records every surviving
support and its sign count.

## Independent 69D review

`review69.py` independently traverses the finite normal-coordinate action,
recomputes moments from factorial identities, compares both coordinate
descriptions, verifies rational weak duality and every geometric pair type,
and runs semantic corruption tests. It now reads the final package-local
certificates, not search outputs. The result is **PASS: 743,522,304 points**.
See `review69_report.json` and the final package
`verifications/improved/dim69-gamma72-head-compression/`.

The final package contains separate reproducible controls for both production
checkers. The independent review preserves the published premises: Gamma72's
catalogue identity, minimum 8, shell population 6,218,175,600, and spherical
11-design property. It does not enumerate the ambient minimal shell.
