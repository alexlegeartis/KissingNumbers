# Exact triangle-switch construction in dimensions 36 and 37

The configuration has **501,813 points in dimension 36** and **505,629 in dimension 37**, increasing each preceding certificate by **126**. Search results are proposals only. The proof uses integer coordinates in a positive-definite Gram basis, with squared norms 12 and distinct inner products at most 6.

## Retained source

The unchanged source is `../completion/certificate.json`, SHA-256

    b37691e110c73f1100112aea1a4dc39e2a911559a7a49c67c88bb303f6cd8760.

It gives 501,687 points in 36 dimensions and 505,503 in 37 dimensions. Its 636 three-point lift-body orbits are used at both heights, contributing 3,816 latitude points. The older `../source/vectors.npz` contains an earlier equator; neither new verifier uses that old equator or trusts its cached base. Both rebuild the base from the completion certificate and freshly enumerated lattice shells.

Let G be the source's positive-definite integral Gram matrix in dimension 12. The corresponding Coxeter–Todd lattice L has minimum squared norm four, established by complete exact short-vector enumeration. The full Gram matrix in dimension 36 is diag(G,G,G).

The source base has 403,704 points in 1,429 outer cosets:

- The zero outer coset supplies 2,268 axial norm-twelve points.
- 672 outer cosets of squared-norm pattern (6,6,0) supply 217,728 points.
- 756 outer cosets of pattern (8,4,0) supply 183,708 points.

The retained source proof checks every unchanged pair through exact quotient distances and phase products. Both new verifiers freshly replay that proof in an isolated temporary directory. Their short-shell enumerations find exactly 756 norm-four roots, 4,032 norm-six vectors and 20,412 norm-eight vectors. These finite populations are computed, not inherited from a published census.

## Three-point base switches

A norm-(8,4,0) owner consists of a three-point root triangle X_a={x_0,x_1,x_2} and 81 large vectors y in a norm-eight quotient coset. Its points place x_s in a small block i, y in a large block j, and zero in the remaining block k.

For each y independently, replace these three points by the points with -x_s in block k, y in block j, and zero in block i. There are 756×81=61,236 available switches. The frozen certificate uses 90, replacing 270 base points by 270 new ones.

Each point retains squared norm 12. The three triangle-internal products are 8-2=6. Distinct y,y' in the same quotient coset satisfy y dot y'<=2, so any old or switched triangles in that owner have product at most 2+4=6. Different-owner endpoint compatibility follows from the exact quotient-distance checks. The verifiers also directly check every changed base point against all retained base points and every other changed point, so the actual witness does not depend on a heuristic interpretation of this argument.

Triangle IDs are 81 times the owner's position among the source norm-(8,4,0) rows, plus the large vector's index in the owner's lexicographically ordered 81-vector list. `switches.json` independently records all 90 actual large vectors, root indices and coordinate placements. The independent verifier requires its fresh decoding to agree exactly with this record before deleting any point.

## Equatorial replacement

A five-integer phase orbit (a,b,c,u,v) represents

    (x_(a,t), x_(b,u+t), x_(c,v+t)),  t=0,1,2,

with phases reduced modulo three. Each point has squared norm 12. The proposal removes 1,469 source phase orbits and adds 1,511 new ones. Thus the final number is

    403704 + 3×(32661-1469+1511) = 501813.

Every added equatorial point is checked against the entire final base and final equator. Every changed base point is checked against the retained base and final equator. This explicitly includes the new-equator versus unchanged-base category. All unchanged/unchanged pairs inherit the freshly replayed retained source proof. No numerical angle tolerance or sampled pair test enters the proof.

## Extension to dimension 37

All 636 source lift-body orbits survive. They give 1,908 distinct vectors of squared norm eight. Appending either height +2 or -2 gives 3,816 additional points of squared norm 12.

Distinct body products are at most two. Equal-height products are therefore at most 2+4=6. Opposite-height products are at most 8-4=4, including the case of equal bodies. The 36-dimensional configuration lies at height zero, so its product with a latitude point is its product with that lift body.

Every changed base point and added equatorial point is checked against every lift body. All actual latitude pairs are also checked by the independent verifier. The final count is therefore

    501813 + 3816 = 505629.

Both heights are retained explicitly. The lift extension is verified for this witness, not assumed merely because the 36-dimensional count increased.

## Two proof implementations

`verify.py` uses a fresh enumeration from the retained rational shell enumerator, reconstructs the complete source base, and checks the new configuration through exact compressed block-product tables. It also checks norms, distinctness, and exact ranks 36 and 37.

`verify_coordinates.py` imports no producer code, search code, gain verifier, or source enumerator. It implements its own rational LDL shell enumeration, reconstructs the entire base, checks all decoded switches, and directly evaluates every affected point-pair category with bounded int64 products. It checks full norms and uniqueness in both dimensions, including the actual latitude points. Both routes separately replay the unchanged source proof.

`negative_controls.py` restores deleted original base or phase points and detects products 8 and 10, respectively, exceeding the allowed 6. These demonstrate that the verifiers reject specific invalid alterations rather than merely confirming the point count.

Run `python ../verify.py` for the complete package replay in temporary copies. Optional explicit coordinates can be written from a writable copy with `python verify.py --write-coordinates --prefix exported`. The coefficient rows use the accompanying Gram matrix; normalization by sqrt(12) gives unit vectors with distinct products at most 1/2.

The original source certificates and manuscript remain unchanged. The dated source review establishes improvement over the pinned repository construction, not exhaustive worldwide priority.
