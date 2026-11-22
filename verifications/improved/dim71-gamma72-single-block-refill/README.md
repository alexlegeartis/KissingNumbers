# A single-block refill in 71 dimensions

The configuration has **2,623,656,288** points, an increase of **77,352**
over the [previous four-slice construction](../dim68-91-lattice-code-continuation/mixed_shell71_20261117/PROOF.md).
It removes 4,600 points from the small norm-eight slice and adds 81,952
new points. The dominant norm-ten slice is unchanged.

## Construction

Use the [Γ72 Leech-block coordinates](../dim73-82-gamma72-depth6-lattice/).
A vector is a triple of integral Leech blocks, with half the usual direct
sum metric. Modulo two its projected A labels all agree and its B labels
sum to zero. A and B are complementary totally singular halves of the
Leech quadratic space. Only the B-only classes can occur on occupied
blocks when another block is zero. The exact Leech norm-four shell has
no A-zero labels.

Take the normal x=(2λ,0,0), where λ is the explicit norm-four vector in
[pattern.json](pattern.json). Write (s,j) for a lift z of squared norm
s and inner product j with x. Project orthogonally to x-perpendicular
and normalize to unit length. The old slices are (8,3), (8,4), (10,4)
and (10,5). Their old population proof applies to this new normal: the 11-design moments depend only on the norm of x, and the discrete inner-product support is fixed by the minimum-8 arithmetic. The pinned verifier checks those moments and counts, so no previous orientation is assumed. Replacing the last by its translate under z↦z−x gives
(8,−3) with exactly the same projections.

Remove all (8,4) lifts of the form (2η,0,0), where η has Leech squared
norm four and ⟨λ,η⟩=2. Add the projected lifts y=(2μ,0,0), where μ
is a Leech norm-four vector orthogonal to λ and the A label a of λ+μ
is nonzero. Require its B label b not to be the sum of two B labels
of norm-four Leech vectors with A label a. The verifier checks precisely
this predicate across the complete 196,560-vector Leech minimum shell.

## Cross-pair proof

Any raw norm-16 triple in Γ72 has sorted block norms (0,0,16),
(0,8,8), (4,4,8) or (4,6,6). Indeed each Leech norm is even with
nonzero minimum four. If one block vanishes then all A labels vanish;
norm four is forbidden by the shell check, and norms six and ten
are forbidden since B is totally singular. In the single-block case
the nonzero block is even. Two Leech norm-four vectors cannot have
inner product ±3, as their sum or difference would have norm two.

For z in (10,4), a conflict with y requires ⟨y,z⟩=5: the projected
squared radii are both 8, the integer product is unchanged because
y is perpendicular to x, and Γ72 minimum eight bounds it by 5.
So d=z−y is a raw norm-16 triple with first-block products
⟨λ,d₁⟩=4 and ⟨μ,d₁⟩=−3. The possible first-block squared norms
are 4, 6, 8, 16. At norm four, d₁=λ gives μ-dot zero. At norm six,
d₁−λ would have norm two. At norm eight, d₁−λ would have norm four
and μ-dot −3. At norm sixteen the first block is even, inconsistent
with the odd μ-dot. Thus there is no conflict with this entire slice.

For z in (8,±3) or (8,4), a conflict requires ⟨μ,z₁⟩=4.
Indeed a dot of at most 3 is already below half the product of the
projected radii. At first-block norm four this forces z₁=μ and j=0.
At norm six, z₁−μ would have norm two. At norm eight, z₁−μ has norm
four and λ-dot j; ±3 is impossible and j=4 forces z₁=λ+μ.
Our filter excludes this last case: A is nonzero so both other
blocks must be occupied and each norm four, contradicting the
specified B-label condition. At norm sixteen, the only block is
even. Hence j cannot be odd; at j=4 the point was one of the
4,600 deletions.

The new lifts already lie in x-perpendicular. Two distinct ones
have inner product at most 4 by the parent minimum or the Leech
bound, so they form a kissing code. The previous verifier proves
all unchanged pairs and their populations. These cross checks also
exclude coincidences with the new normalized rays.

Run with Python, NumPy and SymPy:
`OPENBLAS_NUM_THREADS=1 python verify.py`.
The script pins and replays both previous inputs, computes the precise
addition and deletion counts, and checks the block-norm possibilities.
The results are lower bounds, not optimality.
