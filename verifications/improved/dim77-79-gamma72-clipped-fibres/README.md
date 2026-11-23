# Clipped Γ72 fibres in dimensions 77–79

Reuse the exact [80D clipping rule](../dim80-gamma72-partial-fibres/) after restricting the vertical root stock to the D5, E6 and E7 sections of E8. Every clipped edge is directed to the designated endpoint; the rejected contacts are a union of its exact raw norm-eight collision filters. The retained raw norm-twelve shell includes the 4+4+4, 4+8 and 6+6 types. Whole-coset depths handle all other pairs. All fibres are used with their antipodes, the full Γ72 shell and the pure tail roots.

| Dimension | Previous bound | New bound | Full / clipped antipodal fibre pairs | Gain |
|---:|---:|---:|---:|---:|
| 77 | 6,222,276,760 | **6,223,198,744** | 12 / 8 | 921,984 |
| 78 | 6,223,917,240 | **6,225,814,200** | 21 / 9 | 1,896,960 |
| 79 | 6,226,104,558 | **6,226,634,478** | 24 / 9 | 529,920 |

Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` with NumPy and SymPy for all three exact replays. The public [seven-chart parent](../dim77-79-gamma72-seven-charts/) is pinned. These are constructive lower bounds, not optimality claims.
