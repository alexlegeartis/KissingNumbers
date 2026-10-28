# Phase refill and partial latitude orbits

28 October 2026. The certificates prove **501,687 points in 36D** and
**505,503 points in 37D**, respectively 15 and 15 above the preceding main-branch
certificates. These are improvements over explicit witnesses, without an
optimality or comprehensive literature-priority claim.

The retained [phase proof](../phase/proof.txt) establishes the starting geometry.
Both new checkers freshly enumerate the short vectors with rational branch bounds,
prove the minimum four of K, and verify `M G M^T = 3G`, `|det M| = 729` and the
surjective quotient-label map with kernel L. Consequently L has minimum 12.
The quotient weights and multiplicities are `(0,1)`, `(4,3)`, `(6,18)` and
`(8,81)`, with populations 1, 252, 224 and 252. All norms below are squared norms.

The unchanged fixed base has 403,704 points. The certificate replaces 147 of
the retained original three-point phase orbits by 165 others. Its 32,661 orbits supply
97,983 points, so the new equator has `403704 + 3*32661 = 501687` points.
Base/base and phase/base inequalities follow from exact quotient distances.
Phase pairs with three distinct labels have distance at least `4+4+4=12`.
Every other phase pair is checked by all three relative phase products, or,
independently, by matrix products of expanded coefficient vectors in the Gram basis.

For the latitude body, label **252 denotes a zero block**, rather than the
ordinary quotient label zero. Nonzero indices 0 through 251 index the increasing
list of weight-four quotient labels. An entry `(a,b,c,p,q)` denotes

`(X[a,t], X[b,(p+t) mod 3], X[c,(q+t) mod 3])`, for `t = 0,1,2`.

Exactly one block is zero. Its phase and the first nonzero block's phase are
fixed to zero; this makes each entry canonical. Every orbit has three distinct
norm-eight vectors. There are 636 such orbits over 216 cells, giving **1,908
distinct bodies**. The old 207 complete cells account for 621 of those orbits;
the continuation adds 15 orbits over nine further cells, with no old lift removed.
Allowing part of a cell is essential to this certificate.

Append each height `+2` and `-2` to every body. Norms become `8+4=12`.
For distinct bodies at one height, the exact body dot bound is **2**, so the
37D dot is at most `2+4=6`. At opposite heights Cauchy–Schwarz gives at most
`8-4=4`; this also covers the two heights of the same body. The lift/base
distance bound eight gives body dot at most `(8+12-8)/2=6`.

A weight-four block has dot at most two with a block from a different quotient
label. Thus a lift/equator product can exceed six only if both nonzero lift
labels equal the corresponding equator labels. The verifier explicitly checks
every such pair. In this witness there are **zero** such pairs. It is not
necessary to multiply all other lift/equator pairs: their bound is `4+2=6`.
All lift/lift orbit pairs, including pairs from different zero-block positions,
are checked; the coordinate checker independently checks all 1,819,278 distinct
body pairs. The total is `501687 + 2*1908 = 505503`.

The symbolic checker and coordinate checker share the certificate and the
mathematical lattice-distance argument, but share no enumeration or angle code.
The latter independently expands all phase points and latitude bodies. Its
61,667,370 phase pair products include repeated checks of pairs sharing more
than one label; all pairs not covered by the different-label proof are included.
No floating-point comparison, optimizer output or inherited shell population
is used in either route. The historical coefficient NPZ files describe the
preceding witness, and are not inputs to these new proofs.

## Search logic and relation to PackingStar

[PackingStar's paper](https://arxiv.org/abs/2511.13391v4) searches discrete Gram
data by adding candidates and removing obstructing choices. Its
[26 October repository update](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6)
reports joint point removals, midpoint insertions and auxiliary layers in
26–31D. The transferable idea is to optimize the surviving and inserted pieces
together and allow partial families. The new search here is an independent
finite-graph heuristic; it does not reproduce or train PackingStar's RL system,
and uses no PackingStar implementation or coordinate files.

The old phase configuration occupied 3,979 cells. Checking its fixed base allows
729,000 cells and 6,561,000 three-point equator orbits. Filter candidates by their
current blockers; with at most four blockers the initial pool has 45,923 orbits.
Construct their exact conflict graph, extend greedily, use `(1,2)` replacements,
explore equal-size exchanges, and periodically force outside candidates while
removing their blockers. The resulting equator trade removes 120 and inserts
133 orbits. For 37D, expose individual latitude phase orbits instead of requiring
all nine points of a body cell. Exact compatibility and the same correction/refill
search add 15 orbits. The proof rechecks the complete result independently of
these choices.

`delta.json` and `reproduce.py` recover the shipped certificate exactly from the
retained witness. The optional `search/` tools expose the larger search pools.
Search files and graphs are scratch data; they are not trusted certificates.

## Further refill and broader base choices

Relative to the preceding 501,672 / 505,488 main-branch certificate, the new
phase trade removes 79 orbits and inserts 84. Freezing the latitude labels and
relative phases during preparation preserves all 636 lift orbits. Both exact
proofs and byte-for-byte delta reconstruction pass. `negative_controls.py`
requires both proofs to reject a duplicate phase orbit and false counts.

The optional `search/base_choices.py` exposes a larger base model. For each
of the 756 base cells with component norms (0,4,8), let i hold the norm-four
label a, j the norm-eight label b and k the zero label. Its alternative keeps
b at j and puts -a at k, leaving i zero. Negation is in F_3^6.

Exact coset distances show that each alternative is compatible with every
original base cell except its own owner, and that all distinct alternatives
are mutually compatible. Each replaces 3*81 = 243 points. Thus any subset of
these 756 replacements preserves the 403,704-point base. The script checks
every one of these finite distance facts before writing a proposed seed.
This is a new search freedom; it is not used by the current bound and is not
a claim of a successful refill. Seed retention removes all phase or lift
orbits lacking the corresponding base-distance bound. Both complete proofs
must verify the final result after any search.
