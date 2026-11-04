# Generated real stabilizer parity sectors

`groups.real_stabilizer_sector(m, parity)` defines a complete family in ambient dimension `2^m`, where `m >= 1` and parity is 0 or 1. Its parameters are mathematical family parameters, not switches between benchmark dimensions. The backend imports no record coordinates or cardinality constants and materializes no states.

For each affine r-flat A in F2^m with r congruent to parity modulo 2, choose its canonical RREF coordinate system. Include every sign function `(-1)^q` for squarefree Boolean polynomials of degree at most two, including the constant term. Set coordinates outside A to zero and divide each sign row by `sqrt(2^r)`. This gives unit vectors, including their antipodes.

There are `[m choose r]_2` linear r-subspaces and `2^(m-r)` affine cosets of each. There are `2^(1+r+r(r-1)/2)` quadratic sign polynomials. Their evaluations are distinct because squarefree Boolean monomials form a basis of the function space. Consequently the exact layer cardinality is

```
[m choose r]_2 * 2^(m-r) * 2^(1+r+r(r-1)/2).
```

The generator counts Gaussian binomials by the product quotient. The verifier independently rebuilds them with the q-Pascal recurrence. It also checks the q-binomial product identity at 1 and -1, giving equal parity-sector counts

```
2^m * product(1 + 2^j for j = 0,...,m-1).
```

The pairwise bound follows from exhaustive cases. Disjoint supports are orthogonal. Distinct equal-dimensional affine supports intersect in at most `2^(r-1)` points, giving absolute normalized inner product at most 1/2. Unequal support ranks in one sector differ by at least two, so intersection size alone gives the same bound.

For equal supports, the sign difference is quadratic. If its polar alternating bilinear form has rank k, its Gauss sum S satisfies `S^2 = 2^r * sum_{h in radical} (-1)^(q(h)+q(0))`. The radical restriction is linear, so `|S|` is zero or `2^(r-k/2)`. An alternating form has even rank; positive rank is therefore at least two and gives bias at most 1/2. Rank zero gives a balanced nonconstant linear function, identical vectors, or antipodes. Identical vectors are excluded from the pair condition and antipodes have inner product -1. Thus every distinct pair has inner product at most 1/2.

Full span is also certified. The even sector contains the signed coordinate axes (r=0). The odd sector contains both sign choices on every two-point affine line (r=1); summing or subtracting those two vectors gives each coordinate axis.

The implicit certificate records the definition, support-rank range, derived layer counts, normalization, all pair-case conditions, and full-span witness layer. Verification reconstructs these data and facts; saved status is never authority. It rejects changed definitions, proof conditions, dimensions, counts, and nonintegral or Boolean parameters. Ambient dimension and count bit bounds are checked before exponentiation. Materialized-point and pair budgets are not consumed because no point list or pair matrix is allocated.

The default recipe `discovery/recipes/stabilizer16.json` sets m=4 and parity=0 and derives 4320 points. The same operator derives 146880 points at m=5 and 9694080 at m=6; these larger controls are below the repository comparisons. This certifies constructions only, with no novelty or unrestricted kissing-number upper-bound claim.

Run `python -m unittest discovery.tests.test_stabilizers -v`. The tests independently enumerate XOR-closed affine flats and quadratic signs at m=1,2,3, checking all normalized pair inequalities with integers. They also check larger derived counts, bounded symbolic execution, altered parameters/counts/proofs, and runner/verification dispatch.
