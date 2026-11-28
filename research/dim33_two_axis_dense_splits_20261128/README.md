# Two-axis dense splits in the 33D kissing construction

The incumbent full 33D code is unchanged. These are exact deformations of its **dense component**, not new lower bounds.

Put the missing dense coordinate on physical axis 0. For selected words of the pinned binary [32,17,8] carrier, choose one body axis u and a head line (1, ±1) in coordinates 0,u. Replace the word by the two opposite signs of that head at amplitude η, and retain the 31 other body signs at amplitude a. Here a=19/25 and η²=8809/1250, giving squared norm 32. The same-parent dot is 62a²−32<16. Any unmodified dense parent has at least seven punctured Hamming differences, so its dot with a replacement is at most 17a+η<16.

The one-pencil certificate assigns two disjoint affine 64-word sets of punctured distance at least 14 to orthogonal heads. This adds 128 dense points. The multipencil certificate assigns 256 selected parents to six pencils. The exact verifier checks all 130,816 new/new dot products in the rational quadratic field, giving 131,328 dense points in total.

An incumbent octad missing either marked axis automatically clears: max dot is ≤2(η+7a)<16. This concentrates conflicts, but does not save the existing sparse stock. For the single-pencil certificate there are 36 old octads containing both axes. On such a support, a fixed head sign faces the body sign projections of one affine 64-word set. Its six-bit projection has rank ≥3, since rank ≤2 would leave a binary [25,≥4,14] shortened code, excluded by the Griesmer sum 14+7+4+2. Every matching head/body center and its six single-bit neighbors conflict. If the projected affine rank is r<6, match the six-cube along one coordinate whose unit vector is outside the projected difference space. The forbidden centers and their partners reserve 2^r distinct deleted matching edges; an allowed code can use at most one endpoint of every remaining edge. Rank 6 rules out the whole 32-slot head cell. Across four head patterns the loss is at least 32 per support, hence at least 1,152 on the 36 fixed supports, even allowing arbitrary signs on each. This rules out a gain of this one-pencil recipe without changing supports, not arbitrary two-axis carrier redesign.

## Replay

Use NumPy and the [pinned incumbent witness](https://github.com/alexlegeartis/KissingNumbers/blob/abda686407447f0f2770d43077bf2b95fb8c29bf/verifications/improved/dim33-root-star-exchange/witness.npz). The checker hashes this witness and proves the dense claims, independently of the conditional sparse cost argument.

    python verify.py witness.npz onepencil.json multipencil.json
