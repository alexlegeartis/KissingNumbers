# Weighted-shell and mixed-rotation search in 15D

This campaign adds construction-search machinery and a restricted finite scan. It does **not** claim a new kissing-number lower bound. The verified control has 2,564 points, the known 15D lower bound. A solver timeout is not an upper bound on the unrestricted problem.

The source is represented by integer vectors of squared norm 16 in the metric diag(1^14,2), with distinct inner products at most 8. Its original record coordinates were taken from the Cohn table data; see `verifications/improved/dim73-95-gamma72-caps/README.md` in this repository. `source.npz` contains the exact source and metric.

`SATURATION.md` proves that the four pure families already attain their joint ceiling of 2,564. Of the 30,024 pool points, 29,778 lie in those families and 246 lie outside; a strict improvement requires an outside point. `search.py` seeds from those 246, adds proven per-family saturation cuts, and performs individual-point shared-deletion exchanges. The checked pool has 30,024 distinct vectors with denominator 12 in the same metric: it combines exact rational Householder images, including images that change the weighted coordinate, contact reflections and low-blocker weighted-shell samples. Each source deletion is charged once, and candidate incompatibilities are enforced jointly. A positive result is written as `candidate.npz`; it is a proposal until the independent verifier succeeds. A random seed and a per-neighborhood time budget are fixed in the script, so an unchanged result is not a certificate of impossibility.

## Check and run

Requires Python 3, NumPy and SciPy. Keep at least one logical CPU free as required by the repository's `AGENTS.md`.

    python research/weighted_shell15_20261105/verify.py
    python research/weighted_shell15_20261105/search.py --trials 1000 --seconds 3
    python research/weighted_shell15_20261105/verify_candidate.py candidate.npz

The first command checks every source pair, every stored shell point and its exact blocker list, and all 30,024 search candidates and blocker lists. The candidate verifier requires more than 2,564 vectors and independently checks every unordered pair in the specified integer metric.

## Complete integer shell

`scan.cpp` enumerates every integer vector of squared norm 16 in this weighted metric against the fixed source. The dynamic-programming count is 135,725,444; the two partitions contain 67,380,484 and 68,344,960 vectors. Neither partition contained a free extension. At most six blockers occur for exactly 17,572 vectors: the 2,564 source points with one blocker, 672 non-source vectors with four blockers, and 14,336 non-source vectors with six blockers. The checked output `low_points.npz` records these vectors and their exact source blockers. The sample output alone does not prove enumeration completeness: rerun the scan to reproduce it. This fixed-source, integer-alphabet result does not bound the 15D kissing number.

In an empty scratch directory:

    python /path/to/research/weighted_shell15_20261105/prepare_scan.py
    g++ -O3 -std=c++17 /path/to/research/weighted_shell15_20261105/scan.cpp -o scan
    ./scan 0 > part0.log
    ./scan 1 > part1.log

The expected final histograms, in blocker-count order 0 through 6, are `0,1540,0,0,672,0,0` and `0,1024,0,0,0,0,14336`. Each run writes `scan_Q<part>.bin` and `scan_B<part>.bin` in the current directory. The scan is computationally substantial; concurrent partitions must leave a CPU free.
