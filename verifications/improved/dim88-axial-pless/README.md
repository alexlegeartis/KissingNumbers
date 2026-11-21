# Axial Pless and square layers in dimension 88

| Layer | Points |
|---|---:|
| Antipodal QR48 chord section | 12,397,455,648 |
| Sparse Pless layer | 959,965,380 |
| Induced square design | 333,312 |
| D88 roots | 15,312 |
| **Total** | **13,357,769,652** |

This improves 13,029,906,432 by **327,863,220**. Run `python verifications/improved/dim88-axial-pless/verify.py` with NumPy and assertions enabled. The replay uses the earlier [86–87D axial Pless package](../dim86-87-axial-pless/) and the [rotated-square support generator](../dim90-91-a2-rotated-square/). As in those packages, the two pinned ternary parent codes' minimum distance 15 is inherited.

1. Fix QR48 coordinates 0 and 3 to zero and delete them. Restrict coordinates 1, 2, 24 and 28 to trits 1 and 2 and center their chords. The actual six-column projection has rank six. The 42 unit triangles and four squared-norm-3/4 chord entries give total squared norm 45 and dimension 88. Centering leaves distances unchanged; the possible distances of the parent words are multiples of three between 15 and 45, so one can use both antipodal images, disjoint by the full triangles.
2. Shorten the verified weight-15 Pless orbit union at 0 and 3. Its 194,304 words have exact axial occupancy histogram (37,468, 80,272, 58,152, 16,784, 1,628). Map full-plane nonzero trits to signed, phased A2 roots and axial nonzero trits to ±√3. The 86–87D inequalities apply unchanged: every nonzero slot contributes at most 3/2 against a dense point, every distinct outer trit costs at least 3 squared distance, and phase words for the same outer word differ in at least five of the full-plane roots (cost 9 each). The exact weighted phase census is 959,965,380.
3. Take the induced four-point supports of the order-46 Steiner quadruple system on 42 planes. Each of the 2,604 supports carries 128 even-parity sign patterns, with two coordinates of absolute value 1/√2 in each occupied plane. The basis is as in the rotated-square proof: square unit directions lie at least 15° from the six old sparse root directions. Against sparse vectors the squared dot is at most 24+12√3 <45; against dense vectors the dot is at most 4 <√45. The square/square test uses at most two shared planes for distinct supports, or at least two differing sign coordinates for equal supports.
4. In the same frame, take all D88 roots with two signed unit coordinates. Their mixed dot against either ternary layer is at most √12 because each root meets at most two blocks of squared norm at most 3. Against squares it is at most √2, exactly half the product of the radii. Different layers are disjoint by their different numbers of occupied full planes.

The earlier 88D construction is not counted as a sublayer.
