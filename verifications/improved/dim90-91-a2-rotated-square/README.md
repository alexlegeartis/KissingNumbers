# Rotated-square and root augmentations in dimensions 90 and 91

| Dimension | Preceding bound | Net increase | New lower bound |
|---:|---:|---:|---:|
| 90 | 21,943,640,934 | 459,540 | **21,944,100,474** |
| 91 | 22,359,345,896 | 459,898 | **22,359,805,794** |

Use the [pinned 90D ternary QR/Pless construction](https://github.com/alexlegeartis/KissingNumbers/tree/0470d09c6a62bddb0aa5294d5bc322332148d7b7/verifications/improved/dim90-92-hexagonal-phase-layers) and the [pinned 91D latitude continuation](https://github.com/alexlegeartis/KissingNumbers/tree/0470d09c6a62bddb0aa5294d5bc322332148d7b7/verifications/improved/dim91-92-direct-latitude-layers). Keep all 90D incumbent points and all 91D incumbent points except the two pure-head poles. The remaining incumbent points have squared norm 45, each body plane is empty, a unit triangle direction, or a root of squared norm 3, and the 91D latitude coordinate has absolute value at most sqrt(3). Every retained point occupies at least 14 body planes.

## New square layer

By [Hanani’s construction](https://doi.org/10.4153/CJM-1960-013-3), a Steiner quadruple system on 46 points has exactly one block through each triple. We construct it using the order-16 binary affine design followed by the 3n−2 expansion. Drop the distinguished point and every block containing it, then identify the remaining 45 points with the physical planes. The 3,465 remaining supports meet pairwise in at most two planes. On each support choose the 128 sign strings of length eight with even parity, and put coordinates of absolute value 1/sqrt(2) on both axes of every occupied plane. These vectors have squared norm 4. Different signs on the same support differ in at least two coordinates, giving dot at most 2; different supports also give dot at most 2 by their at-most-two overlapping planes.

In each body plane take orthonormal axes at angles −30° and 60° relative to one of the fixed sparse roots. The unit square directions are then at 15°, 105°, 195°, and 285°. The six signed old root directions lie at multiples of 60°. The minimum angle between these sets is 15°, so the worst square/sparse dot is at most 4 sqrt(3) cos(15°), whose square is 24 + 12 sqrt(3) < 45. For the old dense layer the dot is at most 4. With squared radii 4 and 45 both estimates give normalized dot at most 1/2. The 91D extra coordinate is zero for all squares, so the same check covers the retained latitude layers.

## New root layer

Use the standard D_n roots in this orthonormal frame, including the single added head axis when n=91. These are the 4 binomial(n,2) choices with two nonzero coordinates, each ±1, and squared norm 2. Distinct root dots are at most 1. Each new root occupies at most two incumbent blocks, even if it involves the head coordinate. The squared norm of the retained old point restricted to those blocks is at most 6; Cauchy–Schwarz gives raw dot at most sqrt(12), strictly below half the product of radii. Against a square the raw dot is at most sqrt(2), exactly half the product of radii. Root and square rays are distinct by coordinate supports, and both are distinct from all retained incumbents by body-plane support. The old pure poles are omitted in 91D; they are not counted together with the root layer.

The result is 443,520 square rays and respectively 16,020 or 16,380 root rays. Run `python3 verify.py` to enumerate supports, exhaust the finite parity and angle tables, and check the exact numerical inequalities. The old witnesses and their proofs remain pinned separately. These are lower bounds, not an exhaustive global-priority claim.
