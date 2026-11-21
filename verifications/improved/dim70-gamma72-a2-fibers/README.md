# A2 normal fibers of Gamma72 in 70 dimensions

**Kissing number at least 1,361,137,216**, a gain of 2,496 over 1,361,134,720. The three replayed populations are 610,915,680, 610,915,680 and 139,305,856.

The construction takes three full fixed-pairing classes of minimal Gamma72 vectors, then projects and normalizes. It uses no fitted planar head coordinates. The candidate classes were identified by an exploratory harmonic-theta calculation; the released proof counts their populations directly and does not rely on that calculation.

Let the Leech Gram matrix be `F` in the pinned Gamma72 block model. For column labels (zero based) 5 and 22 take the two pure Gamma72 vectors \(n=(0,0,2e_5)/\sqrt2\) and \(m=(0,0,-2e_{22})/\sqrt2\). Their Gram is \(A=\begin{bmatrix}8&4\\4&8\end{bmatrix}\). Their span has rank two. Keep Gamma72 minimum-8 vectors \(z\) at one of three pairing vectors \(p=(z\cdot n,z\cdot m)\):
\[
(-1,-1),\quad (-1,0),\quad (-2,-1).
\]
Project onto the orthogonal complement and normalize. The projected squared norms are \(8-p^\mathsf{T}A^{-1}p\), namely \(47/6,47/6,15/2\). For distinct source minima with pairings \(p,q\), their original dot product is at most four, hence their projected dot product is at most \(4-p^\mathsf{T}A^{-1}q\). The six exact rational pair-type checks in `verify.py` prove all normalized dots at most \(1/2\), both within and across distinct types. This also proves no two different selected source vectors yield the same normalized point.

The computation uses the [pinned Gamma72 construction and Leech completion theorem](../dim68-91-lattice-code-continuation/public_gamma72_provenance_20261117/THEOREM.md). A source point \((a,b,c)/\sqrt2\) has pairing vector \((cFe_5^\mathsf{T},-cFe_{22}^\mathsf{T})\). For each pairing vector the replay constructs a saturated affine integer representation for \(c\), with two eliminated pivot coordinates in an integral determinant-one minor and all other coordinates free. A checked unimodular basis change accelerates the traversal. It runs the existing signed-128-bit leaf norm checker on this rank-22 kernel and weights each valid leaf by the exact inherited completion table. A pure \(c=2u\) contribution would require both pairing entries to be even, so none occurs here. Distinct affine parameter tuples give distinct \(c\), hence distinct source points across their completions.

Run `python verifications/improved/dim70-gamma72-a2-fibers/verify.py` with NumPy, SymPy and g++. It stages pinned source files in a temporary directory, checks the completion-table hash, then enumerates three fibers. Floating-point interval traversal may omit leaves: only lower counts are needed, and each accepted leaf satisfies the exact integer quadratic bound and the exact affine pairing equations. The published Gamma72 minimum and pinned norm-8 completion theorem are inherited premises. This is not an upper bound or a priority claim.
