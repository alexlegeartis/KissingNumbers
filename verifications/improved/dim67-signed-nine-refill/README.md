# Signed nine-support refill in 67 dimensions

**355,170,821 kissing points (+75,375).** Starting from the verified [355,095,446-point affine sign construction](../dim61-67-compatible-sign-families/), add 82,485 points on nine-coordinate supports and retire only the 7,110 conflicting points of the old signed four-coordinate layer. The dense, weight-sixteen, axial, and all other signed four-coordinate points are unchanged.

## Geometry and count

Use squared norm 64 and absolute amplitude 8/3 on each of the nine occupied coordinates. The certificate has 1,835 full supports, each carrying 40 sign patterns from a puncture of the binary Gray image of the quaternary pentacode; it also has 2,404 partial supports carrying 9,085 explicitly stored sign patterns. The small full code is constructed from the four cyclic patterns in [Conway and Sloane, 1994](https://neilsloane.com/doc/Me185.pdf) and checked for Hamming distance three, as are each of the explicit partial sign lists.

On one support, the normalized distance condition follows from Hamming distance at least three. On different supports, let k denote the number of shared coordinates and d the number of signed disagreements there. The unnormalized dot is (64/9)(k−2d). For k≤4 this is automatically at most 32; for k=5 or 6 the projected sign sets must be disjoint, and for k=7 or 8 they must also have no cross-pair of Hamming distance one. Both independent replays check the required sign conditions for every overlapping support pair.

Every new support meets every old weight-sixteen support in at most six places, hence dot at most 6·2·(8/3)=32. Dense products are at most 9·(8/3)=24 and axis products at most 8·(8/3)<32; old dense and axial geometry is inherited. For an old four-coordinate point, only a support contained in the new support and with all four signs matching can conflict. Both replays compute the exact union of these forbidden signed old points. For a full nine-coordinate support containing an old Netto block, relabel four positions of the full sign code so it projects onto exactly 12 of the 16 old signings. Identical relabeling on that same Netto block lets several additions share these removals. Partial supports keep their actual admissible signs rather than paying to delete entire phase classes.

These changes give 355,095,446 + 82,485 − 7,110 = **355,170,821**. The construction is a lower bound, not a claim of optimality.

## Replay

Run the following from this directory, with Python, NumPy and a C++17 compiler:

    OPENBLAS_NUM_THREADS=1 python verify.py
    OPENBLAS_NUM_THREADS=1 python audit.py

Both rebuild the 42,269 old weight-sixteen support masks and 11,792 four-support masks from pinned existing sources. The Python checker uses exact packed integer overlaps and projected-sign sets. The C++ audit independently reconstructs literal signed supports and compares all 680,849 signed pairs on support pairs with overlap at least five. Both count the actual signed four-points retired. Numerical optimization, exhaustive enumeration of the old 355-million-point code, and external data downloads are not required.
