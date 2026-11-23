# Rephased Γ72 half-root shells with exact nonlattice tail appends

The fixed vertical geometry of the [published half-root shells](https://github.com/alexlegeartis/KissingNumbers/tree/0cabadac8febbdd7b2477a28103f2ccf2792c6df/verifications/improved/dim77-84-gamma72-half-root-shells) admits alternative horizontal phase assignments. The shifts are all in the second polarized Leech half; they preserve the first-label pair and therefore the 136,704-point norm-six count on every valid minimum tail. The new counts come from evaluating all finer shell layers exactly. A separate finite list of nonlattice vertical fibre pairs is then included wherever compatible. Each append has the vertical geometry tested against the resulting rephased shell, not just the original parent.

| Dimension | New lower bound | Gain since pinned half-root shell |
|---:|---:|---:|
| 77 | 6,222,013,324 | 192 |
| 78 | 6,223,119,668 | 224 |
| 79 | 6,224,782,016 | 1,304 |
| 80 | 6,228,106,200 | 1,248 |
| 81 | 6,228,379,802 | 1,216 |
| 82 | 6,228,927,150 | 1,296 |
| 83 | 6,229,748,820 | 1,392 |
| 84 | 6,231,393,584 | 2,000 |

These are lower bounds, not optimality claims. The counts and phase offsets are stored explicitly. Each append row supplies tail squared norm v, its three Leech labels, an eight-integer vertical direction X, and the contact count per sign. The actual directions are ±sqrt(v) times the normalized vector (X_0,…,X_3,sqrt(3) X_4,…,sqrt(3) X_7), with zero sidecar. The minimum inter-fibre separation is checked via integer comparisons by bounding squared horizontal depth plus squared vertical distance, for both signs of every tail. The original half-root shell is unchanged in dimension 79, and rephased in the other dimensions.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` from within the package in the repository. To replay elsewhere, supply the paths to the original half-root and depth-six packages in that order. Requires NumPy and SymPy. Input pins and exact output comparisons are checked.
