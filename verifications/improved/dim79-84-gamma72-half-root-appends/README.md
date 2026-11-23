# Nonlattice Γ72 tail fibres in 79–84 dimensions

The half-root shells at [repository commit 0cabadac](https://github.com/alexlegeartis/KissingNumbers/tree/0cabadac8febbdd7b2477a28103f2ccf2792c6df/verifications/improved/dim77-84-gamma72-half-root-shells) can be extended by mixing full norm-four horizontal fibres with additional pure vertical antipodes. All earlier horizontal and attached shell points are retained.

| Dimension | New exact lower bound | Gain |
|---:|---:|---:|
| 79 | 6,224,782,016 | 1,304 |
| 80 | 6,228,105,432 | 480 |
| 81 | 6,228,379,066 | 480 |
| 82 | 6,228,926,334 | 480 |
| 83 | 6,229,747,908 | 480 |
| 84 | 6,231,392,064 | 480 |

These are constructive lower bounds, not upper bounds or optimality claims.

Each witness row is [v, (a,b,t), X, count]. The tail squared norm is v=4 or 8. For its eight integer coordinates X, write H(X)=(X_0,…,X_3,sqrt(3) X_4,…,sqrt(3) X_7). The two tail positions are ±sqrt(v)/||H(X)|| H(X), with zero sidecar. The horizontal vectors are all the points of squared norm 8-v in the indicated integral-block coset of Γ72. Coset depth in this subgroup is read from the exact Leech norm-four lists and anisotropic mod-2 classes. The counts for horizontal norm four come from precisely one norm-eight Leech block or a pair of norm-four blocks. The zero coset contributes only the zero vector at v=8.

The verifier reconstructs the entire prior active tail stock, including every nonempty fibre in the sidecar. For two tail positions of squared norms v and w, the required squared horizontal depth of the difference of phases is 8-v-w+2 times their tail inner product. Every such condition is checked in both signs by a rational integer square comparison, after checking the nonnegative side. No numerical tolerance or claim of optimality enters replay.

Run: `OPENBLAS_NUM_THREADS=1 python verify.py` from this directory of the repository. To replay elsewhere, append paths to the exact half-root and original depth-six packages, in that order. This requires Python, NumPy and SymPy.
