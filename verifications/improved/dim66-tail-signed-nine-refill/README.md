# Tail-supported signed nine-coordinate refill in 66D

**342,373,178 kissing points (+76,006).** Starting from the verified [342,297,172-point nonary construction](../dim61-66-nonary-affine-caps/), keep all 138 previous nine-coordinate supports, and add 80,829 new points on 5,624 supports involving one or both of the two previously unused nonary tail coordinates. Retire only the 4,823 *actual signed points* of the old four-coordinate layer that conflict with the added points. All other old layers are unchanged.

The new points have amplitude 8/3 on nine coordinates and squared norm 64. The 1,685 full sign classes use the 40-word punctured binary Gray image of the pentacode. The 3,939 additional supports use explicit sign lists of minimum Hamming distance three, verified without assuming that they are subsets of that canonical code. Some were obtained by greedily extending the pentacode restrictions or packing the remaining sign patterns independently. Across supports, if k coordinates coincide and d signs disagree there, the raw dot is (64/9)(k−2d), so the comparison is exactly k−2d≤4. The two checkers replay these separations explicitly.

Each new support meets every one of the 35,994 old weight-sixteen supports in at most six positions. The two tail coordinates are tracked in that check, not dropped from the old masks. Dense and axial dots have ceilings 24 and 64/3. Each selected nine-coordinate support contains at most one of the 11,088 retained Netto four-supports. A signed new point then conflicts with precisely the old four-point having that support and matching all four signs. The checkers compute the exact union of these old points. For a full new support, the four signs are assigned to a projection of size 12 instead of 16; using the same projection on a repeatedly hit old block shares those removals. An explicit partial support only pays for its actually used sign patterns.

Thus 342,297,172 + 80,829 − 4,823 = 342,373,178. This is a lower bound, not a proof of optimality.

Run with Python, NumPy and a C++17 compiler from this directory:

    OPENBLAS_NUM_THREADS=1 python verify.py
    OPENBLAS_NUM_THREADS=1 python audit.py

The Python checker independently compares packed support masks and projected sign sets. The scalar C++ audit compares all 1,543,952 literal signed pairs on relevant overlapping support pairs and rederives the same count. The named parent data are pinned; no numerical optimizer, tolerance, external downloading, or expansion of the entire code is required.
