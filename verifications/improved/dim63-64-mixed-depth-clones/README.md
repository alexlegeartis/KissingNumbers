# Mixed-amplitude dense-owner clones in 63D and 64D

This replaces the preceding dense-owner clone construction in each dimension;
it does not add on top of its points. Use its unchanged cyclic/sparse stock, its
section in 63D, and its exact exceptional weight-sixteen sign classes.

Each chosen cyclic dense owner is removed once and replaced by two or more
spikes at distinct marked coordinates. Reverse its sign on each spike's
mark. A shallow spike has unscaled magnitudes 11 on the mark and 2 on the
other 63 coordinates; a deep spike has 14 on the mark and 2 elsewhere.
Normalize every spike to norm eight. Each mark uses only one amplitude and
each owner has at most one shallow spike.

Shallow marks strictly improve the same-mark compatibility: when the heavy
sign agrees, ordinary binary distance 24 suffices for shallow pairs, while
deep pairs require 28. When the heavy sign disagrees, the old cyclic distance
16 suffices for either amplitude. Two shallow clones of the same owner would
fail, so always use at least one deep mark in an owner's clone set.

All other old layer comparisons use exact rational envelopes. If the mark
is contained in a weight-sixteen sparse vector, one mismatch suffices; the
same sign reversal as in the previous construction certifies the default
parity, and the exceptional classes are checked explicitly. All distinct
owners at distinct marks are compatible by the cyclic distance 16. The two
original poles on each used axis are retired.

The Python verifier pins the earlier clone certificate, reconstructs its
actual old stock, checks every owner and exceptional parity class, and
verifies every same-mark pair using integer arithmetic. The two win*.txt
files contain a dense owner bit mask, the mask of its deep marks, and the
mask of its shallow mark (or zero), all in decimal.
