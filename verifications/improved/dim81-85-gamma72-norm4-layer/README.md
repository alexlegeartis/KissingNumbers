# Squared-norm-four tail shells and pure eight-sets in 81D–85D

| Dimension | Constructive lower bound | Gain over complete D_m fibres |
|---:|---:|---:|
| 81 | 6,234,658,944 | 44,192 |
| 82 | 6,238,411,492 | 72,128 |
| 83 | 6,242,536,524 | 72,416 |
| 84 | 6,247,163,912 | 107,520 |
| 85 | 6,251,849,384 | 153,888 |

Keep the entire [complete D_m tail construction](../dim81-85-gamma72-apn-complete-dm/).
Add each tail of squared length four from the list below on a nonzero B-only phase of the same seven-dimensional shift space. Each such phase has exactly 144 shortest horizontal points, squared length four: three local placements of 48 signed norm-eight vectors in the raw metric. Hence these add 144 contacts per tail. The new phase leaves the nonzero colour of every old mixed phase difference unchanged, so every such horizontal distance is squared at least six. Every new tail has its two largest absolute coordinates summing to at most two, hence has vertical squared distance at least two from an old norm-two tail and dot at most four with an old pure tail. The old horizontal norm-eight shell clears by 4+4.

The extra tails are all signed axes of length two and the sixteen signings of each selected four-set, developed from cyclic starters in the verifier. Distinct supports meet in at most two coordinates. The intersection graph colours allow disjoint supports to reuse one phase for each sign parity. All axes share another nonzero phase. Tails on different phases have dot at most two and horizontal separation at least four; on a common phase distinct tails are nonpositive and need no horizontal separation. Points over the same tail have the parent lattice separation.

Add pure tails with eight nonzero signs of amplitude one and even parity, on one eight-set in dimensions 81–83 and three eight-sets in 84–85. Their intersections are at most four, and every pair in one support disagrees at least twice. They have zero horizontal projection and dots at most four with the old pure tails and every new squared-norm-four tail. These add exactly 128 contacts per eight-set.

Run `OPENBLAS_NUM_THREADS=1 python verify.py`, requiring NumPy and SymPy. It checks the pinned parent, the 127 B-only horizontal phases, all new norm-four tail pairs, and the support counts. The construction still uses Nebe's Γ72 lattice and makes no optimality claim.
