# Balanced signed nine-support refill in 63D

**184,426,297 kissing points (+8,513).** This construction extends the verified [balanced 63D nonary section](../dim61-66-nonary-affine-caps/). It retains all 2,084 previous nine-coordinate points and adds 14,790 new points on 746 supports, while retiring only 6,277 actual signed directions from the affine four-support layer. All other old points remain.

The new points have amplitude 8/3 on nine coordinates and squared norm 64. The marked eight-position hyperplane is the exact one in the earlier 63D section. On its intersection with a new support the positive and negative entries balance, so all added points lie in that same 63-dimensional hyperplane without any further projection or scale change. The signed lists have minimum Hamming distance three on each support. Across supports, if k positions coincide and d signs disagree there, the raw dot is (64/9)(k−2d), so the exact comparison is k−2d≤4. Both checkers replay the required phase separations explicitly.

Every new support meets each of the 32,801 old weight-sixteen supports in at most six places, so the cross-dot ceiling is exactly 32. Dense and axis ceilings are 24 and 64/3, respectively. The inherited four-support block stock is the affine SQS after the already-pinned 63D permutation, not the ordinary 64D labeling. A conflicting pair requires an old quadruple entirely contained in the new nine-position support, identical signs on those four positions, and a signed old point which survives the section balance condition on the marked octad. Both replays count the exact union of such old points, not whole unsigned quadruple blocks.

Thus 184,417,784 + 14,790 − 6,277 = 184,426,297. This is a lower bound, not a proof of optimality.

Run with Python, NumPy and a C++17 compiler from this directory:

    OPENBLAS_NUM_THREADS=1 python verify.py
    OPENBLAS_NUM_THREADS=1 python audit.py

The Python verifier uses exact support masks, and the scalar C++ audit independently checks all 170,792 signed pairs on relevant overlapping supports and rederives the count. Existing dense and sparse geometries are inherited from the pinned parent; no optimizer, tolerance or expansion of the entire code is needed.
