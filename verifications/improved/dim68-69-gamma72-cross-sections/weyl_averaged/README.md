# Weyl-averaged affine sections of Gamma72

This package proves the following improvements over the immediately preceding repository records:

| Dimension | New bound | Previous repository bound | Gain |
|---|---:|---:|---:|
| 68 | **471,884,288** | 471,366,656 | 517,632 |
| 69 | **689,452,922** | 689,442,146 | 10,776 |

The construction uses the same compatible six- and ten-offset geometries as the previous affine package, but averages over finitely many orientations of an actual root subsystem. It does **not** assert that the Gamma72 shell is invariant under that group. The claim is existence of a favorable orientation, not enumeration of a particular shell section.

## Reproduce

Python 3 and SymPy are required. No numerical solver or research cache is needed.

    python verify.py 4
    python verify.py 3
    python negative_controls.py

The two small certificates have respectively eight and fourteen nonzero dual rows. The 68D certificate has denominator 24 and exact objective 471,884,288. The 69D certificate has denominator 990,904,320 and exact objective 689,452,922. The verifier uses exact integer arithmetic and checks the ceiling of each rational dual objective; both final objectives happen to be integers. Twenty deliberately corrupted certificates must be rejected.

## Published inputs and exhibited normals

The inherited source is Nebe's extremal even unimodular 72-dimensional lattice Gamma72, with squared minimum eight, 6,218,175,600 minimal vectors, and an 11-design minimal shell. These published facts are not established by shell enumeration here.

- G. Nebe, *An even unimodular 72-dimensional lattice of minimum 8*: https://arxiv.org/abs/1008.2862
- Catalogue: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html
- The exact catalogue matrix is included in `embedding.json`, along with its inherited original repository blob and source information. The verifier pins the complete embedding file's SHA256.

`embedding.json` contains integer lattice vectors in the catalogue basis. The first three and four vectors have the required A3 and D4 Gram matrices. Exact matrix multiplication verifies these embeddings. An additional fifth vector records the D5 research realization and is not needed for either published bound.

## Why averaging is legitimate

Write the chosen normal vectors as rows of U, with Gram G. For each normal i, the root reflection acts on profiles by the integer matrix

    T = I - (G[:, i] / 4) e_i^T.

The verifier checks that every generator is integral and satisfies T G T^T = G. In D4 it also includes the two adjacent permutations of the three mutually orthogonal leaf normals; these preserve G as well. Every transformed tuple T U therefore consists of actual lattice vectors of squared norm eight. The positive-definite integral isometry group of G is finite: each integer row has bounded positive-definite norm, so only finitely many matrices are possible.

For every group element h, let n_h(c) count source vectors having profile c against the actual tuple h U. Each such distribution satisfies the same exact design moments, one-point marginals, antipodal identities, and minimum-vector translation identities. Average these distributions over the finite group. The averaged distribution is constant on profile orbits and still satisfies all those linear constraints. This does not require any group element to preserve the ambient lattice or its shell.

For the selected offset set C, the objective sum over c in C of the averaged count equals the arithmetic mean of the actual union populations over all orientations h. Consequently at least one orientation attains at least that average. The rational dual lower-bounds the average. Since every actual union population is an integer, at least one orientation attains its ceiling, including when the average bound is not integral.

The checker reconstructs the profile domain from all norm-eight roots in the normal span. A coordinate equal to plus or minus eight forces a known normal vector. All other profiles lie in [-4,4]^k and obey the root inequalities. Exact group generators permute this domain and the known roots. Graph components combine these averaging orbits with the exact antipodal and translation identities. A separate standard-library checker, `verify_orientations.py`, enumerates all 1,152 orientations in 68D and 48 in 69D, including simultaneous normal negation. It confirms that the reported groups are exactly the action-only orientation orbits, without relying on the translation identities. The main verifier invokes this checker. There are sixteen orbit variables in 68D and twenty-two in 69D.

## Exact dual and geometry

The independent checker obtains mixed moments from coefficients of powers of the Gram quadratic polynomial, rather than the discovery Gaussian recurrence. It regenerates every certificate row, checks the known root contributions and one-point moments, then verifies the integer dual slack and exact objective. A successful floating-point optimization is not part of the proof.

For offsets a and b, put q(a,b) = a^T G^-1 b. Project each selected source vector onto the normal-space complement and normalize it separately. Its projected squared norm is 8-q(a,a). Distinct source vectors have inner product at most four. The checker verifies positive projected norms and

    4-q(a,b) <= 0, or
    4 (4-q(a,b))^2 <= (8-q(a,a)) (8-q(b,b))

for every same-offset and cross-offset pair. Thus the normalized images have pairwise inner product at most one half. The same inequality excludes image collisions. It is unchanged under every normal orientation because G is unchanged. The projected complement has dimension 72-k.

## Verification scope

Both proofs were independently audited before publication. `verification3.json` and `verification4.json` record exact replays. `negative_controls.json` records twenty rejected corruptions. `SHA256.json` seals the package. `source_review.json` states the comparison and audit scope.

The new package, table/registry integration and retained stronger affine replays were checked. The complete historical repository audit was not run because a complete checkout with historical data was unavailable. The earlier central-section and affine packages are retained. No paper or workflow files are changed. This is an improvement over checked repository records, not a worldwide-priority or optimality claim.
