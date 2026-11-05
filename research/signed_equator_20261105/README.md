# Signed equators and coupled 13D search

This is a construction-search campaign and a restricted-family proof. It does **not** claim a new kissing-number lower bound. The verified control is a switched 1,154-point 13D configuration. The generic kissing-number problem remains open in these dimensions.

`BOUND.md` proves that a compatible family of signed four-support vectors in 12 coordinates has at most 816 members, even when individual signs are selected independently. All 24 coordinate axes may be added, giving at most 840 points in this restricted family. The proof uses a leave-parity argument for 3-(12,4,1) packings and averages over global signs; it is not an upper bound on the unrestricted 12D kissing number.

`SATURATION.md` proves that the saved 13D source plus arbitrary points in three specified pure families (signed four-support equators, two 12-bit head codes with tail ±2, and coordinate axes) cannot exceed 1,154. It uses the classical exact bound A(12,4)=144, plus the verified 24 disjoint axis-versus-source-blocker pairs. This conclusion does not apply to other vector shapes.

`search.py` is an actual mixed construction search, not a proof of a global upper bound. It jointly selects individual points from all 8,192 heads, 7,920 equator vectors, 26 axes and 43 additional quadratic-field candidates, while charging each deleted source point once. It encodes head conflicts by Hamming-ball cliques and equator conflicts by signed-triple owners. The source-saturation inequalities are valid even in the presence of the 43 additional candidates. A solver timeout or unchanged objective is not an impossibility certificate. The result `candidate.npz` is only a proposal until the independent verifier passes.

The source is represented exactly as A+B√3, with integer arrays in `source.npz`. Each source vector has squared norm 16; the kissing threshold for distinct vectors is 8. The extra candidate file uses the same format; its first 26 entries are the axes. The control derives from the Zinoviev–Ericson 13D construction, followed by exact neutral exchanges.

## Reproduce

Requires Python 3, NumPy and SciPy. Keep at least one logical CPU free, as required by the repository's `AGENTS.md`.

    python research/signed_equator_20261105/verify.py
    python research/signed_equator_20261105/verify_low.py
    python research/signed_equator_20261105/search.py --build-only
    python research/signed_equator_20261105/search.py --seconds 3600
    python research/signed_equator_20261105/verify_candidate.py candidate.npz

The first verifier checks every source pair, the 816/144/144/2/48 decomposition, all 4,096 sign-pattern incidences, the 24 disjoint axis blocker pairs, and extra-candidate norms. The second checks the saved shell output and every recorded blocker directly. The candidate verifier requires more than 1,154 points and checks norms and every unordered pair using exact integer arithmetic in Q(√3).

## Full field-shell scan

A separate finite scan enumerates every vector in Z[√3]^13 with squared norm 16 against this source. It does not cover rational denominators or other algebraic fields. Its dynamic-programming count is 505,844,378. The two partitions contain 245,829,626 and 260,014,752 vectors. No vector was a free extension. The 2,762 vectors with at most four blockers consist of the 1,154 source points, 264 two-blocker candidates and 1,344 four-blocker candidates; all non-source points are in the pure head/equator/axis families. `low_points.npz` preserves those outputs. Checking that file is not by itself a proof of scan completeness: rerun the enumerator to reproduce the enumeration.

From an empty scratch directory, run:

    python /path/to/research/signed_equator_20261105/prepare_scan.py
    g++ -O3 -std=c++17 /path/to/research/signed_equator_20261105/scan.cpp -o scan
    ./scan 0 > part0.log
    ./scan 1 > part1.log

The expected final histograms by blocker count 0 through 4 are `0,866,24,0,384` and `0,288,240,0,960`. Each run writes `scan_Q<part>.bin` and `scan_B<part>.bin` in the current directory. The search is computationally substantial and can run partitions concurrently if CPU affinity leaves one logical CPU free. The output is a result for this fixed source and alphabet only.

`search_den5.py` expands the same search to a common denominator of 5, including 50 additional rational Householder-image candidates in `extras_den5.npy`. It imposes a strict gain of at least one and requires at least one candidate outside the saturated pure families. Its source, candidate norms and all conflicts are still exact; verify `candidate_den5.npz` with the same candidate verifier. For example:

    python research/signed_equator_20261105/search_den5.py --build-only
    python research/signed_equator_20261105/search_den5.py --seconds 3600
    python research/signed_equator_20261105/verify_candidate.py candidate_den5.npz

`search_local.py` uses the same denominator-5 candidate families in adaptive individual-point neighborhoods. It seeds only from the non-pure candidates and applies the exact source-category saturation cuts, while allowing coordinated replacements in every pure family. Its output `candidate_local_den5.npz` must pass `verify_candidate.py`; as with the global searches, a timeout or failure to improve proves no global bound.

    python research/signed_equator_20261105/search_local.py --trials 1000 --seconds 3

## Complete six-blocker field scan and radical extensions

The same full 505,844,378-vector field shell has also been scanned with a six-blocker cutoff. Its partition histograms, in order 0 through 6, are `0,866,24,0,384,132,1464` and `0,288,240,0,960,1044,1720`. In addition to the ≤4-blocker points described above, the complete output contains 360 non-source points outside the saturated head/equator/axis families: 216 have five blockers and 144 have six, and all 360 have nonzero √3 components. `low6_points.npz` stores all 7,122 output points and `verify_low6.py` checks all norms, every stored source blocker, the classification and inclusion of the 360 radical candidates in `extras_complete_den5.npy`.

The full expanded denominator-five search has 476 extra points, including the axes, the 360 radical points, previously sampled field points and 50 rational Householder-image points. `search_complete6.py` is the global mixed model with a strict gain target; `search_local_complete6.py` is its individual-point neighborhood counterpart. Both retain the pure-family saturation cuts. Verify any output with `verify_candidate.py`.

    python research/signed_equator_20261105/verify_low6.py
    python research/signed_equator_20261105/search_complete6.py --build-only
    python research/signed_equator_20261105/search_local_complete6.py --trials 1000 --seconds 3

To rerun the complete ≤6 scan, prepare `source.bin` as above and compile `scan6.cpp`; run partitions 0 and 1. As before, saved output verification is separate from the time-consuming completeness rerun. No new lower bound is claimed.
