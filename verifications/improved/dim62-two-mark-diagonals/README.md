# Two-mark sparse heads in 62 dimensions

**172,134,412 kissing points**, improving the [previous verified 62D bound](../dim62-63-guarded-rephase/) of 170,279,948 by **1,854,464** (about 1.09%). All predecessor dense points, affine-plane points and poles remain. Every predecessor sparse support and its population remain, with some sign cosets rephased.

## Orthogonal two-mark heads

The predecessor uses four marked sites `(0,9,20,41)` in the 64-coordinate Reed–Muller parent. Its geometric body is the 60 ordinary coordinates and a two-dimensional head. On the sparse scale (squared norm 16), the single-site heads are `(0,1),(1,0),(0,-1),(-1,0)`. To a support using two nonopposite marked sites, assign the sum of their heads, a directed diagonal of squared length two, and place signs of magnitude one on its other fourteen sites. This still has squared norm 16. For two source supports of weight 16 and binary distance at least 16, their total overlap is at most eight. The head dot is at most the number of shared source marks: distinct source mark directions have nonpositive dot. Thus many two-mark banks fit immediately.

The dense predecessor uses six heads `(2,0),(1,√3),(−1,√3),(−2,0),(−1,−√3),(1,−√3)` at squared norm 64 including their 60 body signs. A diagonal sparse head has exactly one dangerous dense label: its head dot exceeds two, but never exceeds four. Orthogonality of each support to all parent dense rows is checked. The body signs are then assigned total parity opposite the sum of the two corresponding dense label bits. One ordinary mismatch pays two in raw dot, meeting the mixed limit of 16. The verification treats single-mark banks by the same rule. The unmarked banks meet the dense layer by the raw sixteen-site bound. All new coordinates have magnitude at most one, so old signed affine planes and poles are compatible by their support sizes on the sparse scale.

## Joint sign guard proof

The selected layer contains 3,554 two-mark banks, eight single-mark banks, and thirteen unmarked banks. It includes 277 extra shared-eight-site guards to retain every old bank and repair atypical overlaps. Every new/new and new/old pair is checked exactly. A guard with opposite parity between two banks forces at least one disagreement on an ordinary coordinate, because every deleted dummy bit has sign +1 in both banks. This repairs every raw dot above eight, which the replay verifies is at most ten. Previous guard signs remain in force.

Assign each bank its original 16 positions, with one or two erased mark positions fixed to binary zero as necessary. Include total parity and all its guard rows. Their homogeneous span has only weights zero, eight and sixteen. Its nonconstant rows give balanced columns (all nonzero binary characters have zero sum on the columns); refine these equal-sized column classes by bisection to obtain distinct four-bit labels. Adding these binary checks completes the span to the full extended Hamming dual. The checker tests the assigned affine parities against every dependency, so the completed coset has 2,048 signs, and it loses a factor two per fixed dummy. Indeed the completed dual has no weight-one or weight-two word, so its kernel projects surjectively on any one or two fixed dummy positions. The sign code has distance at least four on the remaining ordinary coordinates. Counts are consequently 512, 1,024 and 2,048 signs per two-, single- and zero-mark bank.

The exact gain is `3554 × 512 + 8 × 1024 + 13 × 2048 = 1,854,464`.

## Replay

Run `python3 verify.py` in the repository (Python 3 and a C++17 compiler). The verifier replays hash-pinned predecessor inputs and checks every new pair and all local parity constraints. The only new witness is `merged62.txt.gz.b64`: gzip then base64 of an ASCII list, beginning with the counts of banks and guards, followed by each 64-bit support mask and each guard's two support masks, common mask, and second-bank phase. No floating-point tolerances or optimizer are needed. The number is a construction lower bound, not an upper or optimality statement.
