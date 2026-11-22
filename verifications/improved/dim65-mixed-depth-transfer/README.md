# Mixed-amplitude dense-owner transfer to 65D

**Lower bound: 336,258,432 (+740).**

The preceding 65D construction uses the extended BCH parent P, its
doubly-even radical D (the equatorial binary code), the height-four
owner cap E, and the GF(8) ordinary-axis spikes. Permute the first
64 coordinates of these four layers together by identifying their
F64 coordinate labels with the existing physical affine labels.
Their pair comparisons are unchanged. The old sparse stock can be
kept in its previous positions: all of its comparisons with the
permuted layers use only amplitudes and a tail-occupancy bound.

The relabeling permits using owner words in the common subspace of
the radical D and the original cyclic binary code. The verifier
checks that every selected owner is in both and outside E. Use the
same two mixed clone amplitudes as in 63D and 64D. The cyclic
membership and sign flip on the mark certify ordinary weight-sixteen
support parity; all exceptional supports are checked directly.
Against the E cap, the common BCH parent gives binary distance 14
before the one-sign flip, sufficient after the height-four scale.

Only the GF(8) spikes sharing marks with the new clones are removed.
The surviving GF(8) points occupy different heavy axes and their
exterior signs also differ from each owner before its sign flip by
at least 14. The exact scaled dot-product inequalities are checked.

There are 3,242 new clones, 1,574 retired equatorial owners and 928
retired GF(8) points. The count is 336,257,692 + 3,242 − 1,574 − 928
= 336,258,432. All 5,253,661 new–new comparisons are checked exactly.

Run with Python 3 and NumPy. No optimizer or discovery files are needed.
