# Weighted-shell and mixed-rotation search in 15D

This campaign adds construction-search machinery and a restricted finite scan. It does **not** claim a new kissing-number lower bound. The verified control has 2,564 points, the known 15D lower bound. A solver timeout is not an upper bound on the unrestricted problem.

The source is represented by integer vectors of squared norm 16 in the metric diag(1^14,2), with distinct inner products at most 8. Its original record coordinates were taken from the Cohn table data; see `verifications/improved/dim73-95-gamma72-caps/README.md` in this repository. `source.npz` contains the exact source and metric.

`SATURATION.md` proves that the four pure families already attain their joint ceiling of 2,564. The saturated families include arbitrary real heads at tail height ±2. Of the 52,373 pool points, 47,253 lie in these broader families and 5,120 lie outside; a strict improvement requires an outside point. `search.py` seeds from those 5,120, adds proven per-family saturation cuts, and performs individual-point shared-deletion exchanges. The checked pool has 52,373 distinct vectors with denominator 12 in the same metric: it combines exact rational Householder images, including images that change the weighted coordinate, contact reflections and low-blocker weighted-shell samples. Each source deletion is charged once, and candidate incompatibilities are enforced jointly. A positive result is written as `candidate.npz`; it is a proposal until the independent verifier succeeds. A random seed and a per-neighborhood time budget are fixed in the script, so an unchanged result is not a certificate of impossibility.

## Check and run

Requires Python 3, NumPy and SciPy. Keep at least one logical CPU free as required by the repository's `AGENTS.md`.

    python research/weighted_shell15_20261105/verify.py
    python research/weighted_shell15_20261105/search.py --trials 1000 --seconds 3
    python research/weighted_shell15_20261105/verify_candidate.py candidate.npz

The first command checks every source pair, every stored shell point and its exact blocker list, and all 52,373 search candidates and blocker lists. The candidate verifier requires more than 2,564 vectors and independently checks every unordered pair in the specified integer metric.

## Complete integer shell

`scan.cpp` enumerates every integer vector of squared norm 16 in this weighted metric against the fixed source. The dynamic-programming count is 135,725,444; the two partitions contain 67,380,484 and 68,344,960 vectors. Neither partition contained a free extension. At most six blockers occur for exactly 17,572 vectors: the 2,564 source points with one blocker, 672 non-source vectors with four blockers, and 14,336 non-source vectors with six blockers. The checked output `low_points.npz` records these vectors and their exact source blockers. The sample output alone does not prove enumeration completeness: rerun the scan to reproduce it. This fixed-source, integer-alphabet result does not bound the 15D kissing number.

In an empty scratch directory:

    python /path/to/research/weighted_shell15_20261105/prepare_scan.py
    g++ -O3 -std=c++17 /path/to/research/weighted_shell15_20261105/scan.cpp -o scan
    ./scan 0 > part0.log
    ./scan 1 > part1.log

The expected final histograms, in blocker-count order 0 through 6, are `0,1540,0,0,672,0,0` and `0,1024,0,0,0,0,14336`. Each run writes `scan_Q<part>.bin` and `scan_B<part>.bin` in the current directory. The scan is computationally substantial; concurrent partitions must leave a CPU free.

## Complete denominator-three shape family

The pool also includes all 17,472 signed permutations of the numerator shape `(8,2,2,0^11;6)` divided by 3 in the weighted metric. Against the fixed source, exactly 1,344 have four blockers, 2,688 have five and 13,440 have seven. All 17,472 are new shapes but lie in the broader saturated height layers; the pooled verifier checks each point and its source blockers directly. `enumerate_den3.py` reconstructs this whole family and writes `generated_den3_pool.npz`. The search alternates low-blocker outside-family seeds with seeds from the whole outside-family pool.

For a fixed tail sign, at most 28 points from this family can coexist: their 14-dimensional heads have pairwise nonpositive inner products. If their normalized head Gram matrix is G of rank at most 14, then tr(G²)≤2m and m²≤14·tr(G²), so m≤28. This bounds only this one candidate family and does not constrain general 15D constructions.

## Complete outside-family sparse-head scan

A further exact scan examines all 43,051,008 signed vectors with two head entries ±2, six head entries ±1, six zero head entries and last coordinate ±1. Every vector has squared norm 16 and lies outside the four saturated families. Neither partition found a free extension of the fixed source. Exactly 5,120 vectors have at most ten source blockers: 768 have eight, 1,280 have nine and 3,072 have ten. All 5,120 are in `sparse_low.npz` and are merged into the search pool at load time, and the verifier checks every saved blocker list. They are the search's outside-family seeds; individual saturated-family points can still be selected jointly. This is not a bound on vectors in any other shape or on the 15D kissing number.

To reproduce, use the same `source.bin` prepared above, compile `sparse_scan.cpp`, and run parts 0 and 1 as for `scan.cpp`. The partitions enumerate 23,181,312 and 19,869,696 vectors; their blocker-count histograms 0 through 10 are `0,0,0,0,0,0,0,0,256,256,1024` and `0,0,0,0,0,0,0,0,512,1024,2048`. Output files are `sparse_Q<part>.bin` and `sparse_B<part>.bin`.

`search_direct.py` is a complementary global model over all 5,120 sparse-head points and arbitrary source deletions, without adding other replacement shapes. It enforces a strict gain of at least one; infeasibility or timeout applies only to this fixed finite model. Build or run it with:

    python research/weighted_shell15_20261105/search_direct.py --build-only
    python research/weighted_shell15_20261105/search_direct.py --seconds 900
    python research/weighted_shell15_20261105/verify_candidate.py candidate_direct.npz

## A mixed rational and square-root family

`radical50_low.npz` contains 1,024 candidates from the complete signed family
with two head entries ±3, six head entries ±2 and weighted tail ±2. Their
weighted squared norm is 50; after multiplying by 4/sqrt(50), they have the
source norm. A complete 43,051,008-vector scan found no free source extension;
exactly 1,024 have at most 12 source blockers, each with exactly 12. The two
partition counts are 23,181,312 and 19,869,696, retaining 256 and 768 points.
All retained norms and blocker lists are checked by `verify_radical_low.py`.

`search_hybrid_radical.py` searches those candidates together with the existing
rational pool, charging shared source deletions. Cross-family comparisons use
exact integer thresholds for the different norms. It samples at most 600
eligible candidates in a local model and can time out; no result is a global
upper bound. Independently check any output using
`python verify_radical_candidate.py search15_hybrid_lns_potential.npz`.

To reproduce the full scan from this directory, run `python prepare_scan.py`,
compile `scan_radical50.cpp`, and run its partitions 0 and 1. The scan writes
raw candidate and blocker files. A positive candidate is not a record until
independently checked and compared with current published bounds.
