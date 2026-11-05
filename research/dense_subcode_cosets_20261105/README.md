# Larger binary-sign layers from subcode cosets

The 32–35D source constructions share a 131,072-word binary length-32,
minimum-distance-eight sign code C. Any other binary code of the same length
and minimum distance can replace that layer without changing the ordinary
octad or root layers. In 32D and 33D a replacement larger by 16,384 words
would raise the point count by 16,384. In 34D and 35D one can always remove the
four extra heavy heads and restore the two omitted roots: the same replacement
would then gain 16,382 points. The heavy heads could alternatively be kept if
separate compatibility is verified.

This search samples index-eight subcodes H of C and seeks nine compatible
H-cosets. H has 2^14 words, so the union would have 147,456 words, 16,384 more
than C. It searches quotient syndromes using exact shortest-path distances,
then looks for a clique. It never assumes the resulting union is linear.

    OPENBLAS_NUM_THREADS=1 python research/dense_subcode_cosets_20261105/run.py --trials 10 --seconds 60

A found proposal must pass the independently recomputed syndrome and kernel
minimum-distance checks and a full construction/record comparison. No new
bound is claimed. Failure or skipped large candidate graphs excludes only the
specific attempted subcodes, not this construction paradigm.

With `--codimension 4`, the search instead tries to replace the 16 original
cosets of a 2^13-word subcode by 17 or more; a success would gain at least
8,192 dense points. This is still a finite, resource-limited search.

A campaign using seed 323333 and 80 row operations tested 1,000 sampled
index-eight subcodes. For each, all quotient syndromes and the corresponding
compatibility graph were checked; none contained the required eight nonzero
cosets. This is a finite-sample finding, not an exclusion for other subcodes
or more general nonlinear constructions. To repeat the same campaign, use
`--trials 1000 --seconds 7200 --max-candidates 15000`; a shorter time limit
may stop it before all trials are checked.

If a proposal is found, the script runs `verify_kissing.py`, which rechecks the
binary coset code and each original signed-octad source independently. It uses
the universal dense-to-octad and dense-to-root inner-product bounds, and drops
the optional heavy heads in 34–35D before counting. The resulting counts still
need comparison against current records before a claim of priority.
