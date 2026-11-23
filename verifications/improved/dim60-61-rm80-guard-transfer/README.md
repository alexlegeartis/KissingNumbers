# Transferred signed support banks in dimensions 60 and 61

| Dimension | Published predecessor | Additional | Verified lower bound |
|---:|---:|---:|---:|
| 60 | 68,033,576 | 688,128 | **68,721,704** |
| 61 | 118,639,402 | 210,944 | **118,850,346** |

Keep all dense points, affine four-sets and poles of the [published RM80 puncture constructions](../dim60-61-rm80-puncture). Keep every previous sparse support and its population, jointly rephasing some sign cosets to accept extra supports from existing higher-dimensional guarded witnesses. The proof depends only on the actual support masks, not their discovery or the optimizer.

Both retain the marks (0,20,41,61). In 60D every new support avoids these four marks; in 61D it may also use exactly mark 0 or exactly mark 20, assigning the fixed tail +1 or −1 and shortening the sign code at its removed site. The ordinary sparse banks have odd total parity; the marked banks have even total parity. All new supports are checked orthogonal to the parent dense rows and therefore satisfy the same mixed inequalities as the corresponding old sparse modes. Each sparse point has squared norm 16.

All pairs involving a new support are checked. Their largest unconstrained raw dot is the number of shared ordinary coordinates plus the product of the fixed heights. If it exceeds eight, it is at most ten and there is a shared eight-site guard with opposite sign parity. A removed mark present in the guard is fixed to binary zero in both banks, so the discrepancy occurs on an ordinary coordinate and decreases the dot by two.

In each bank the homogeneous guard space, including total parity, has dimension at most five, and its nonconstant rows all have weight eight. Bisect its equal-size column classes to obtain an extended-Hamming dual. Its affine parities are jointly consistent: every dependency is checked. Thus each bank has 2,048 signings, exactly half retained when it has one fixed dummy coordinate. Earlier guard conditions remain in force.

Run python3 verify.py in the proof directory (Python 3 and a C++17 compiler). It resolves the published 60/61 sibling directory and replays its SHA-pinned files, then verifies every new pair, every sign bank, and both counts in claims.json. Each witness is gzip-compressed ASCII wrapped in base64. No optimizer or floating-point tolerances are needed.
