# Three affine slices: tau(70) >= 1,360,986,112

This exact existence bound improves the preceding repository value 1,249,778,250 by **111,207,862 (8.90%)**. The older central-section construction remains unchanged in substance. The paper and workflows are unchanged. No worldwide-priority or global-optimality claim is made.

## Construction and geometry

Let L be Nebe's extremal even unimodular 72-dimensional lattice Gamma72. Its minimum squared norm is 8, its minimal shell has 6,218,175,600 vectors, and this shell is a spherical 11-design. Choose minimal vectors u1,u2 with inner product 4. Their existence follows from the independently checked one-point distribution: for any fixed minimum, there are 127,800 minima with inner product 4.

Their Gram matrix is G=[[8,4],[4,8]], with inverse (1/12)*[[2,-1],[-1,2]]. For each shell vector v, set c=(v.u1,v.u2). Keep all vectors with c in

    (0,-1), (1,0), (1,-1).

Project them orthogonally onto the common 70-dimensional complement of span(u1,u2), then normalize each projected vector. This is one spherical code in the same space.

For offsets a,b, the removed-component product is a^T G^-1 b and the projected squared norm is 8-a^T G^-1 a. Distinct original shell vectors have inner product at most 4. Thus the normalized projected vectors satisfy the kissing condition if

    4(4-a^T G^-1 b)^2 <= (8-a^T G^-1 a)(8-b^T G^-1 b),

or the unsquared numerator is nonpositive. All six offset-pair types, including equal offsets, pass exactly. The pair (0,-1),(1,0) is tight; the other cross squared margins are 5/2. All projected norms are positive. The same angle inequalities exclude coincident normalized vectors from distinct source vectors, so the count is not reduced by projection.

## Exact population proof

The joint cell counts obey all mixed design moments through degree 10, the exact one-point marginals, antipodal symmetry, and the minimal-vector translation identities: if w and v are minima and w.v=4, then w-v is another minimum. Lattice minimum and positive-semidefinite constraints restrict possible cells; a singular cell contains at most one vector.

There are 61 unknown cells and 108 equality constraints. The objective sums the three selected populations. The integer-scaled dual in `certificate.json` proves the objective is at least

    1360986111473559012792944 / 1000000000000000
    = 1,360,986,111.473559012792944.

Integrality gives at least **1,360,986,112** points. This uses exact weak duality, not a rounded optimizer result. For Ax=b, 0<=x<=u and objective c, any rational y with c-A^T y nonnegative on unbounded coordinates proves the lower bound y.b minus the sum of u_j times negative-slack magnitudes on bounded coordinates. The verifier checks this inequality using integers only.

## Self-contained verification

From this directory, using Python's standard library only:

    python verify.py
    python rejection_controls.py

`verify.py` independently regenerates every possible cell and all 108 constraints. Its Gaussian moment calculation uses a two-variable recurrence, independent of the discovery engine's Wick matching enumeration. It checks one-point moments, lattice translations, primitive integer matrix rows, dual feasibility, the integer objective bound, and every geometric pair type. It imports neither an optimizer nor repository common code. `independent_verification.json` records PASS.

`rejection_controls.py` rejects eight corruptions: bound, numerator, Gram matrix, right-hand side, matrix coefficient, offset, cell, and dual multiplier. Its receipt is `rejection_controls.json`. `SHA256.json` seals the package. `source_review.json` records source scope and comparison provenance.

The inherited inputs are explicitly limited to Gamma72's published existence, minimum, shell population and 11-design property. This package does not enumerate billions of shell vectors. It certifies an existence lower bound and a concrete selection/projection recipe, not an explicit list of all resulting coordinates.

## Sources and comparison

- G. Nebe, *An even unimodular 72-dimensional lattice of minimum 8*, J. reine angew. Math. 673 (2012), 237–247: https://arxiv.org/abs/1008.2862 .
- Nebe–Sloane Gamma72 entry: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html .
- The extremal-shell design input and moment machinery are documented in `../../../../common/PROOF-kpoint.md` and `../../../../common/kpoint_lp.py`.
- The preceding repository bound is the central-cell construction in `../derive.py`; its proof is retained. The comparison snapshot was main commit 8bb5028b4437e12e4efc86cd001a55b22e830121, RESULTS.md blob 1912c0d9aef1fafb6f9db2ee8d9a93cc7f6a65ab.

The package checks and targeted record-table checks passed. A complete repository-wide audit was not run in this publication environment because the full checkout and its historical data were unavailable. This limitation does not affect the self-contained proof above.
