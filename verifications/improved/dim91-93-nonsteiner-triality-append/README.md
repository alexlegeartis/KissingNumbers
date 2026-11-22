# Four-row triality modes beyond the Steiner supports (91–93D)

This construction adds **14,336** points in each dimension, retaining the
entire preceding complementary-parity Steiner append and every older
layer. The bounds are **25,923,338,000** in 91D, **34,513,272,592** in 92D,
and **47,398,174,480** in 93D.

Each record in `witness.txt` gives three disjoint 22-bit masks, in
the order axial, even flat, odd flat, with four occupied rows in total.
None of the supports is a Steiner quadruple. Each occupied row is a
signed axis from the indicated 24-cell frame. The four axis choices
obey a zero-XOR check, and all 16 sign choices are used, so each mode
gives 1,024 distinct vectors of squared norm 16.

For any two triality-tagged modes, a shared row contributes at most
4 if its two frames agree and at most 2 otherwise. The verifier bounds
the resulting integer sum against the old uniform and new
complementary Steiner modes, the retained seven-row parity modes, and
all other modes in this append. For each short matching mode of weight
w and overlap q, it checks the squared mixed bound q² ≤ 2w. The
new directions have zero head; four occupied rows with per-row product
at most four are compatible with the inherited dense and Golay layers
of squared radius at least 88 and 80 respectively. Ordinary single-row
roots have mixed product at most two.

Run `python3 verify.py` in the published sibling layout, with assertions
enabled. It checks the pinned older inventories, the exact 1,024-word
per-mode Gram bound, and the printed totals. The optimizer is not
a premise of the replay.
