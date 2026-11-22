# Weighted-distance subspaces in the 62D omitted fibre

**Lower bound: 126,570,636 (+1,945,600).**

Keep the stock of the preceding 62D omitted-fibre weighted-spike
construction, and use the same amplitude ratio (a=9/2), zero tails,
omitted label 2, and all 60 marks. Replace all sixty rank-15 sign
fibres by rank-16 fibres, with a separate subspace for each mark.

For a fixed mark (p), a binary difference (u) has weighted Hamming
distance (|u|+(a^2-1)u_p). A word of ordinary Hamming weight 16 is
therefore allowed if it changes the marked sign. All nonzero words of
each proposed 16-dimensional difference space have weight at least
16. Every word of weight below 20 changes the marked sign. Thus
the same-mark weighted dot product satisfies the kissing bound.

For different marks the original distance-16 code is sufficient
provided no sign word is used twice. The verifier enumerates the
sixty translated rank-16 spaces and retains each word once, on
the earliest physical axis using it. All remaining pair checks
and old-stock comparisons are inherited: signs come from the same
omitted label and their magnitudes are unchanged, while the same
120 poles are retired.

The sign fibres contain 3,911,680 distinct retained words after
discarding 20,480 repeated occurrences, giving
(122,659,076+	ext{3,911,680}-120=	ext{126,570,636}).

Run `OPENBLAS_NUM_THREADS=1 python verify.py` from this directory.
The verifier uses the pinned omitted-fibre parent verifier and checks
all kernel spans and translations directly. The JSON files specify
explicit generator masks and translations in the ordered basis of the
24-dimensional shortened kernel.
