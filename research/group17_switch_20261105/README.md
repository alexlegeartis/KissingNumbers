# Three inequivalent 5,730-point configurations in 17D

This is a pair of **neutral 128-for-128 exchanges**, not a new lower bound.
The source, once-switched and twice-switched configurations each have 5,730
points, matching the known 17D lower bound. Their overall pairwise inner-product
histograms agree, but their distributions of per-point inner-product profiles
differ pairwise, proving that these three configurations are not isometric.
The two switches can also be applied independently; the verifier checks all
four combinations. The certificate does not claim that the new configurations
are inequivalent to every configuration in the literature.

`coordinates.npz` contains integer source, once-switched and twice-switched rows for the metric
`diag(1,...,1,2)`; every row has squared norm 93,312 and every distinct
inner product is at most 46,656. Divide transformed coordinates by the square
root of 93,312 after identifying the weighted metric with Euclidean space.

The switched points were selected individually from the source and seven
rational isometric images. In the weighted metric, use reflections with
normals e_0+e_16, e_1+e_16, and e_2+e_3+e_4+e_5. The finite words used are
[], [0], [1], [2], [0,1], [0,2], [1,2], [0,1,2], applied left to right. The
normals mix the weighted coordinate into the other coordinates. Whole orbits
were not required: the search charged shared deletions and chose individual
points from the image union.

Run `python verify.py` in this directory. It independently checks every norm
and every pair for all four switch combinations in exact integer arithmetic, checks the replacement points belong to the stated group-image union, and verifies the
per-point-profile inequivalence certificate. `generate_pool.py` rebuilds the
finite image graph, and `python search.py 17w 170705` runs an individual-point
shared-deletion search. The search is heuristic and a failed run gives no upper
bound. The related `search.rational_image_exchanges` operator is registered in
the main discovery engine.

## Complete six-support tier scan

For the twice-switched source, `scan_tier.cpp` tests the entire 1,025,024-vector
family with six head coordinates ±6 and tail ±6, in the norm-288 weighted
integer representation. No free extension occurred. Counts through 12 source
blockers are `0,1024,0,0,0,0,0,0,9216,4096,11264,4096,3072`; all 1,024
one-blocker vectors are source duplicates. The remaining 992,256 candidates
were stopped once more than 12 blockers were found. `tier_low.npz` contains all
32,768 retained vectors and their blocker lists; `verify_tier.py` independently
checks every saved norm and blocker list, and distinguishes duplicates.

To reproduce the completeness scan, run in this directory:

    python prepare_tier.py
    c++ -O3 -std=c++17 scan_tier.cpp -o /tmp/scan17-tier
    /tmp/scan17-tier
    python verify_tier.py

This excludes a free addition from that stated family only. Coupled replacement
and other cap shapes remain open.

## 16D equator eight-support scan

The twice-switched configuration has a 4,320-point 16D equator. In its
norm-eight integer representation, `scan16_eight.cpp` enumerates all
3,294,720 signed eight-support vectors and finds no free extension. There are
3,712 source duplicates and 120,704 other vectors with exactly eight blockers;
the remaining 3,170,304 were stopped after more than eight blockers. All
retained norms and duplicates and 1,000 sampled blocker lists were checked
independently. The complete scan is an exact fixed-source, fixed-family result,
not an upper bound on the 16D kissing number.

From this directory run `python prepare16.py`, compile `scan16_eight.cpp`,
run the binary, and then run `python verify16_low.py` to replay those checks.
