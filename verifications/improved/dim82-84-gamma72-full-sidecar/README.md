# Full orthogonal D4 sidecars for clipped Γ72–E8 fibres

| Dimension | Exact lower bound | Gain over preceding bound | New antipodal tilted directions |
|---:|---:|---:|---:|
| 81 | 6,228,705,122 | tie | 1 |
| 82 | 6,229,251,942 | 64,000 | 3 |
| 83 | 6,230,072,172 | 323,352 | 6 |
| 84 | 6,231,712,632 | 319,048 | 12 |

Start from all fibres of the [80D weighted certificate](../dim80-gamma72-partial-fibres/). Take a full D4 orthogonal sidecar, whose twelve antipodal root directions are paired with explicit depth-six phases of Γ72. These phases are the first sidecar block from the [earlier half-root construction](../dim77-84-gamma72-half-root-shells/sides.json), but here we prove they coexist with the much stronger weighted 80D stock. A1, A2 and D3 ≅ A3 are root subsets of D4 and give the three shorter sidecars. Each tilted root uses the entire 136,704-point horizontal contact shell. Each pure sidecar root, scaled to squared norm eight, adds one point, and its compatibility follows by Cauchy–Schwarz. There are no optimality claims.

The new–old vertical tails are orthogonal and each has squared norm two. Their horizontal differences therefore need squared depth at least four. The new phases are in the original integral-block chart. Six of the old phases are also in that chart and are checked by the ordinary label-depth rule. Each of the other 33 old phases has a raw Leech parity vector split into three blocks. There are three distinct patterns, each with all three individual block classes of minimum squared norm eight. Adding a new phase, which has zero raw parity, does not change this pattern.

A squared horizontal depth below four would require a raw vector of total squared norm below 32: the parent metric is the three-block Leech metric divided by two, and the half-phase contributes a further factor of four. A local representative has squared norm divisible by four in these classes, so the only possible norm triples are (8,8,8) and permutations of (12,8,8). In a local class with reference norm-eight vector v, every other norm-eight vector w is v−2r for a Leech root r of norm four, except v and −v: the squared norms of (v±w)/2 sum to eight. For norm twelve that sum is ten, so one of the vectors (v±w)/2 has norm four; testing both signs completes the norm-twelve enumeration. The exact counts are 48 norm-eight and 4,096 norm-twelve representatives per local class. The verifier tests the four possible triples by XOR of local label lists, without storing their triple products. It then checks all within-sidecar depth rules and totals.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` in this directory. Requires NumPy and SymPy. The replay first validates the pinned 80D certificate and then all comparisons and cross-phase separation.

