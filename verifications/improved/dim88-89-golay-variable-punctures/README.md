# Enlarged shortened Golay layer in dimensions 88 and 89

| Dimension | Prior bound before the new 89D dense construction | New lower bound | Gain |
|---:|---:|---:|---:|
| 88 | 12,977,854,464 | **13,029,906,432** | 52,051,968 |
| 89 | 14,051,596,290 | **15,177,390,082** | 1,125,793,792 |

The 89D count includes the separately verified 14-fiber construction in [the two-color parent](../dim89-two-color-latitudes/) (15,125,338,114 including the previously published sparse layer). The new sparse contribution is 52,051,968 in either dimension, with every old sparse point retained. These are lower bounds only.

The [preceding shortened Golay–dodecacode layer](https://github.com/alexlegeartis/KissingNumbers/tree/3eefe082c42bfdb58ca29d9621c0d58dc07e3bfa/verifications/improved/dim88-95-golay-shortening) used 629 weight-11 outer patterns from one fixed puncture pivot. The new layer checks all weight-12 quaternary words from pairs of binary Golay words and allows the position discarded to get weight 11 to vary with the word. A search over the finite candidate graph gives 353 additional outer words, stored directly as two bit masks each in `patterns_extra.json`. The first two 4D body rows are zero throughout.

An active row is assigned one of three coordinate pairings, with an axis and sign as in the original layer. The row has exactly two coordinates of absolute value 2 and squared norm 8. For two patterns of weights (r,s), the number of jointly active rows plus the number of equal-tag active rows is (r+s-d), with quaternary Hamming distance (d). The mixed row dot is at most 4 for unequal tags, 8 for equal tags; it therefore sums to at most (4(r+s-d)). Both the mutual new pattern distance and the distances to the 629 old weight-11 words and the 11 old weight-10 words are at least 11. The normalized dot constraints follow from raw caps 44 and 40, respectively, under squared radii 88 and 80. The inequality for weight 10 is (4cdot40^2le88cdot80).

The 1,024 shortened additive dodecacode axis words have distance six, and the 144 Julín signs have distance three. So within a single pattern, unequal axis words agree at no more than five positions and contribute at most 40; equal axis words with different signs give at most (8(11-2cdot3)=40). Every new pattern contributes 147,456 points. The dense–sparse pair estimate is unchanged: each new active row contributes at most 4 against any of the three old dense row alphabets, and the dense squared radii are at least 88. Sparse points have zero height. The row supports and coordinate magnitudes also prove that the new vectors do not duplicate other layers; row pairings, axes and signs can be recovered from each new point.

`python3 verify.py` (also requires Node.js) checks every new outer pair and every mixed old/new outer pair directly, plus counts and compatible inherited geometry. The proof of the pinned Z8 homogeneous minimum remains an explicit source premise.

`claims.json` records the exact totals and `verify.log` records a successful replay.
