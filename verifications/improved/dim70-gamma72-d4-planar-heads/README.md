# 70D: D4-profile compression with anisotropic planar heads

**Kissing number at least 1,361,134,720.** This adds 148,608 points (0.0109191%) to the preceding repository construction, 1,360,986,112 at commit e363fc2ce669c24b7b08f057570888db94746b0c. This is a small verified improvement, not an exhaustive worldwide-priority claim. The unrecovered larger historical claim is not used.

## Construction and exact geometry

Use the published minimum-8 Gamma72 lattice and its Leech-block representation identified by the [existing exact source/completion theorem](../dim68-91-lattice-code-continuation/public_gamma72_provenance_20261117/THEOREM.md). The four normals in kernel_0.npz have coordinate matrix R given in replay.py; their span has dimension 4. Each recorded integer profile h has projection h/4 and residual b in dimension 68 with norm squared (128−|h|²)/16. Put B=4b. For the recorded integers u,v, replace each source point by

    (sqrt(3) B, u, sqrt(3) v) / sqrt(T),
    T = 3(128−|h|²)+u²+3v².

These are unit vectors in dimension 70. Distinct source minima have scalar product at most 4, so residual products B·C≤64−h·k. The verifier checks T>0 and T≤384, which gives within-fiber products ≤(T−192)/T≤1/2. Across distinct profiles, put P=3(64−h·k)+u*u'+3v*v'. We check either P≤0 or 4P²≤T*T'. Thus all pairs have scalar product at most 1/2, which also ensures no duplicates. There are 121 same-fiber and 7,260 cross-fiber checks; no pure-head points are added.

## Population audit and reproduction

Run `python verify.py` from any directory (Python, NumPy, SymPy, g++ required). It rebuilds the complete parity-weight table from the existing public integer witnesses, checks its hash and source-matrix binding, then directly re-enumerates all 121 selected profiles. The integral solve minor and complementary kernel minor have determinant ±1, establishing injective saturated coordinates. Exact 128-bit leaf norm checks establish valid nonnegative completion contributions. Pure c=2x contributions are separately counted and disjoint. Floating traversal may omit points: only lower counts are needed, not exhaustive pruning or completeness of the 1,681-profile histogram. No histogram-zero argument is used.

The independent replay reproduced every selected population. The traversal reuses reviewed enumeration code, with freshly rebuilt affine inputs; it is not a second enumeration algorithm. Published Gamma72 minimum 8 and the rootless even-unimodular rank-24 theta-series theorem remain explicit premises, with their exact finite consequences rechecked. This is a spherical kissing construction, not a lattice-only kissing configuration.

Source: G. Nebe, *An even unimodular 72-dimensional lattice of minimum 8*, J. reine angew. Math. 673 (2012), 237–247, https://arxiv.org/abs/1008.2862 ; the exact catalogue Gram and its binding are retained in the linked public source theorem. The prior repository proof remains unchanged at ../dim70-71-gamma72-cross-sections/affine70/.

`check_geometry.py` separately checks the exact inequalities and receipt binding and rejects four deliberately corrupted candidates. This fast check alone is not a replacement for the full source/population replay.

`final.py` emits the registered 70D claim; it is a count driver, not the proof verifier.
