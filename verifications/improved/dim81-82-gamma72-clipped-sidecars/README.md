# Orthogonal sidecars for clipped Γ72–E8 fibres in 81D and 82D

| Dimension | Constructive lower bound | Gain over preceding bound | New tilted directions |
|---:|---:|---:|---:|
| 81 | 6,228,705,122 | 325,320 | 1 antipodal pair |
| 82 | 6,229,187,942 | 260,792 | 3 antipodal pairs |

The [80D certificate](../dim80-gamma72-partial-fibres/) gives 6,228,431,712 points using the Γ72 norm-eight shell and an eight-dimensional E8 tail. It uses 39 of the 42 depth-six phases in the seven-chart construction. The three unused phases are all in one chart.

Append an orthogonal copy of A1 or A2 to the tail space and use its root directions. Every new tilted tail has squared norm two. It is orthogonal to all old tilted tails, so the sum of their squared vertical norms is four; every new–old phase is distinct and has horizontal squared depth at least four. Two new tilted A2 directions have absolute vertical dot one. Of the three phase differences, one has horizontal squared depth six and the other two have depth four. Clip the one opposing fibre against both blockers using exactly the earlier integer short-shell rule. The survivor count is 104,704 per sign. The other two are full, with 136,704 per sign. In 81D the single new phase is full. Same-ray antipodes have squared vertical separation eight.

The pure roots in the orthogonal sidecar, scaled to squared norm eight, give two more points in 81D or six in 82D. They are compatible with all old pure roots by orthogonality and with all tilted roots by Cauchy–Schwarz. The total gain therefore includes exactly these extra pure points. No optimality claim.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` from this package. Requires NumPy and SymPy. It checks the pinned 80D certificate, the phase differences, the clipped count and both exact totals.

