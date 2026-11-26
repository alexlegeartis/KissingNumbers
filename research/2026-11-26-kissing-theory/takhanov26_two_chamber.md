# A two-chamber reduction for a possible +1 in 26D

**Conditional construction, not a record.** The exterior, exact axes and original chamber charts come from [Takhanov v0.2](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2) at 201,212. The [previous theorem](takhanov26_paid_tetrads.md) caps the clipped replacement chamber at 24; its attaining 24 was the open step at first write. It is now [rigorously excluded for the 2513 quartet](takhanov26_2513_23.md); see the [balanced-amplitude alternative](takhanov26_balanced_quartet.md) for a live target.

Pin Takhanov v0.2 77f0860, independently read data d26_blocks.json.

Let Q_A be the known exact 2513 parity quartet with head (0,1,2,3), original tail 9, alternate 8; Q_D the identical (25,±21,±21,±21;23,6)/√2513 parity quartet on the complementary head (4,5,6,7), same tails.

Both coexist: cross numerator identically 23²+6²=565<2513/2.

Within each known exterior chart, 100²<4·2513 gives the 200864 outer points, axes by caps.

Both quartets share mixed tetrad #5=(0,4,8,9), exact deleted sign masks respectively {0,2} and {0,1}, hence only THREE original skeleton rows paid.

No other skeleton losses.

The opposite rebuilt chambers share only tail 9, so any capped rows are cross-safe; each is also cross-safe with the opposite quartet.

The unchanged chambers B/C keep all 26 rows each, but not in arbitrary orientations.

In upstream support order, identity coordinate permutations, sign masks 16 for B and 4 for C are jointly safe against all eight donors, independently checked using 4·max(0,numerator)²≤2513·denominator² for each of the 416 pairs.

Two of those donor/support intersections do have three coordinates.

Each rebuilt chamber has the same clipped 5D shape, constraints 25u+23v+21 max_even s·x ≤√2513/2, and passes all original surviving skeleton and fixed B/C chambers by support intersections ≤2; exterior passes by the exact four-chart exterior lemma (the alternate-tail coordinate h=0 gives Σa_i≤2 automatically under the caps).

Thus a single genuine 24-point clipped 5D code X gives BOTH 24-point chambers by isometry, and exactly 200864+20+224−3+52+8+48 =201213, +1 over Takhanov.

The 24 upper bound is proved, but sharpness is not.

Killing 25 was not the end of this multi-chamber strategy.

Independent replay: run `verify_paid_tetrads.py` on the SHA-256-pinned upstream `d26_blocks.json`. The two middle chambers use zero-based signs in the upstream indexing. Exact arithmetic checks the donor norms, cross products, the three paid old skeleton rows and the two orientations. Original five-coordinate supports intersect pairwise on at most two coordinates; the surviving four-coordinate supports each intersect a rebuilt original-face support on at most two. Thus the coordinate caps suffice for these interactions. The 24-row code still needs to be supplied and certified.
