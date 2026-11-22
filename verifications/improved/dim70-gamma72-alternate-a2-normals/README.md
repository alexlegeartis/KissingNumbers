# Alternate A2 normal embedding in 70D

**Kissing number at least 1,361,142,400**, a gain of **5,184**.

Use the same three A2 pairing profiles as the preceding A2 section,
with a different explicitly recorded pair of pure Leech-block normals.
Their ordinary Leech Gram is [[4,2],[2,4]], so their Gamma72 Gram is
[[8,4],[4,8]]. This leaves every normalized cross-fibre inequality
unchanged, while increasing the sum of the three actual populations.

The verifier checks the two coefficient rows, an exact determinant-one
solve minor, the affine integer parameterization and its unimodular
basis change, and all six rational pairing inequalities. It rebuilds
the pinned norm-eight Leech completion table and safely counts accepted
leaves with exact signed-128-bit norms. Computed populations are
610,923,456, 610,923,456 and 139,295,488. The accepted-leaf argument does
not rely on completeness of floating pruning.

The alternate normals were found by scanning the rare A2 corner
pairing (-4,-4) using exact minimum-shell moments. The verification
does not rely on that discovery argument: it enumerates the three
actual projected fibres directly with the pinned Gamma72 construction.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy, SymPy and g++.
