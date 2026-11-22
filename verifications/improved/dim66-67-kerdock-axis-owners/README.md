# Kerdock-axis dense-owner lifts in 66 and 67 dimensions

| Dimension | Prior | Kissing lower bound | Gain |
|---:|---:|---:|---:|
| 66 | 342,373,178 | **342,381,366** | +8,188 |
| 67 | 355,170,821 | **355,183,103** | +12,282 |

Start from the already certified 66D and 67D stocks. Keep their signed weight-sixteen, four-support and nine-support layers. On the 64 ordinary coordinates use the doubly-even extended BCH radical already checked in the 65D code. It has the same size and minimum distance as the old flat dense layer. All old sparse cross inequalities depend only on amplitudes and support, so changing the dense code does not disturb those layers.

## Dense owner groups

The [65D dense-owner proof](../dim65-dg-dense-owners/) checks the actual RM(2,6) embedding in the BCH radical. Here take the 32 alternating-form masks of its subfamily with second field parameter zero, along with all 128 affine Boolean words. The checker verifies minimum Hamming distance 28 and constructs disjoint translates of the quadratic-mask set with offsets 0, 1, 2. Each chosen extra coordinate gets a separate group of 4,096 dense words, all in the radical.

For each owner word interpreted as its ±1 ordinary body u and assigned tail axis j, remove that equatorial point and add the two points (s u, ±r e_j), where s=√37/8 and r=3√3. Retire the two ±8 pure poles on axis j. Each point has squared norm 64. Same signed axis products are at most 27+(37/64)·8 < 32, by owner distance 28. Opposite signs on one axis give at most 37−27. On different axes the owners are distinct and the product is at most (37/64)·32. The remaining dense points also have the ordinary cap 32, which is contracted.

## Compatibility with the actual prior sparse layers

The verifier reconstructs each published shortened weight-sixteen support in its actual coordinate ordering. Their tail occupancy histograms are 32,794 of occupancy zero and 3,200 of occupancy two in 66D, and 32,669 and 9,600 respectively in 67D. The signed refills do not add any weight-sixteen points. Thus when a lifted point sees its own axis in a weight-sixteen support, at most 14 ordinary entries can contribute. The product is bounded by max(32s,28s+2r)<32. For the four-support and nine-support layers the bounds are max(16s,12s+4r) and (8/3)max(9s,8s+r). The pure ordinary axes give at most 8s. Every bound is checked using rational upper enclosures for s and r.

All lifted points replace their own equatorial word; the new points, retained dense points, and retained sparse points are distinct. Each lifted tail axis therefore gives exactly 4096−2=4094 net points. These are ordinary kissing-number lower bounds, not general optimality claims.

## Replay

Run `python verifications/improved/dim66-67-kerdock-axis-owners/verify.py` with NumPy installed. It checks the pinned siblings, the exact code arithmetic, the true occupied tail counts, the rational geometry, and the total. The old calibrated sparse proofs are inherited.

Optional `owner_caps.py` verifies exact normalized binary Krawtchouk polynomial certificates. It shows that in the fixed uniform 65D two-latitude lift, the doubly-even dense owners cannot exceed 131,072; in each of the fixed uniform 66/67D Kerdock-axis lifts they cannot exceed 4,096. For a distance-code bound, binary Krawtchouk kernels are positive semidefinite by the full binary characters of each fixed degree. The polynomial has constant coefficient 1, positive nonconstant coefficients and nonpositive values at every allowable nonzero distance, giving its value at zero as the size bound. These upper bounds apply to these dense owner classes, not to other tail geometries or other kinds of owners.
