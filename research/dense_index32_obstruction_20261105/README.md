# A large-replacement requirement for one binary subcode

This is an exact statement about one specified index-32 subcode H of the
32D dense sign code C. It is **not** an upper bound on the kissing number or
on other nonlinear-code constructions.

C is a union of 32 cosets of H, each with 4,096 words. For this H, the exact
syndrome graph shows that every H-coset outside C that could join a code of
minimum distance eight conflicts with at least 20 of C's 32 cosets. A complete
check of five- and six-original-coset anchors then shows that any larger union
of whole H-cosets can retain at most four of the original 32. Thus it must
replace at least 28 original families and add at least 29 new ones. This
explains why small coset exchanges cannot succeed in this particular model.

The verifier independently enumerates the full 20-bit quotient, checks the
parent and subcode, and exhaustively peels the relevant compatibility graphs.
It normalizes by an original selected coset, which preserves the number of
original cosets because those cosets form a subgroup. It proves no restriction
on a larger union containing no original coset, on a different H, or on codes
that break H-cosets.

    OPENBLAS_NUM_THREADS=1 python research/dense_index32_obstruction_20261105/verify.py
