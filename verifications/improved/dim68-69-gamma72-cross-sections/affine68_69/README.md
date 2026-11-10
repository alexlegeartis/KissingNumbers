# Affine Gamma72 slice unions in dimensions 68 and 69

Two new exact existence lower bounds extend the successful 70D affine-slice method:

- **tau(69) >= 689,442,146**, versus 627,822,180: gain **61,619,966 (9.81%)**.
- **tau(68) >= 471,366,656**, versus 375,222,168: gain **96,144,488 (25.62%)**.

Both have independent standard-library verifiers. These comparisons use the current repository bounds, not older paper tables. No worldwide-priority or optimality claim is made.

## Construction

Use the minimal shell of Nebe's extremal even unimodular lattice Gamma72, with squared norm 8 and 6,218,175,600 vectors. It is a spherical 11-design. The catalogue Gram and the repository's four-vector witness are included in `embedding.json`, with original source blob IDs and hashes. Exact integer multiplication verifies that these vectors span a scaled D4 with Gram

    [[8,-4,-4,-4],[-4,8,0,0],[-4,0,8,0],[-4,0,0,8]].

For dimension 69, use its first three vectors, whose Gram is the upper-left 3-by-3 block. This realizes the required A3 rather than merely postulating its existence. For dimension 68, use all four vectors.

Given normals u_i with Gram G, each shell vector v has integer label c_i=v.u_i. Select the labels below, project every selected vector orthogonally to the common complement of the normals, and normalize. If a,b are two labels, the removed-component product is a^T G^-1 b and the projected squared norm is 8-a^T G^-1 a. Since distinct original shell vectors have inner product at most 4, compatibility follows from

    4(4-a^T G^-1 b)^2 <= (8-a^T G^-1 a)(8-b^T G^-1 b),

or a nonpositive unsquared numerator. Every label has positive projected norm, and every same-label/cross-label pair satisfies the inequality exactly. This also excludes coincidences after normalization. All selected vectors occupy the same ambient complement, of dimension 69 or 68 respectively.

The ten 69D labels are

    (0,-1,0), (-1,0,0), (-1,-1,1), (0,-1,-1), (-1,-1,0),
    (0,-2,-1), (-2,-1,1), (-1,-1,-1), (-1,-2,1), (-1,-2,0).

The six 68D labels are

    (0,1,-1,0), (0,1,0,0), (0,1,0,-1),
    (1,1,-1,-1), (-1,2,0,0), (0,3,-1,-1).

## Population proof

All mixed moments through degree 10, the exact one-point marginals, antipodal symmetry, lattice-range restrictions, and minimum-vector translations constrain the unknown populations of integer-label cells. The translation v -> w-v is valid whenever w and v are minima with w.v=4.

Exact cell identifications are substituted out; this is equality elimination, not a relaxation or rounded rank reduction. The 69D system has 369 cells, 132 remaining variables and 188 moment/marginal constraints. The 68D system has 1,681 cells, 520 remaining variables and 617 constraints.

The rational dual for the 69D label union proves population at least

    689442145992299599289941204 / 1000000000000000000
    = 689,442,145.992299599289941204.

The 68D certificate proves at least

    471366655617735410113015808 / 1000000000000000000
    = 471,366,655.617735410113015808.

Integer populations give the stated bounds. The optimizer only finds witnesses; feasibility and weak duality are checked exactly.

## Independent replay

    python verify.py 3
    python verify.py 4
    python rejection_controls.py

The verifier imports no NumPy, SciPy, source LP engine, or source equality-reduction code. It independently checks the catalogue embedding, all possible cells, every translated/antipodal equivalence class, all moment and marginal constraints, exact dual feasibility and bound, and all 55/21 geometric pair types. Its Gaussian-moment recurrence differs from discovery's Wick matching enumeration.

`certificate3.json` and `certificate4.json` are compact certificates: matrix entries and right-hand sides are regenerated independently, not trusted as inputs. The proof data need only 12 KB and 52 KB. `embedding.json` supplies the exact catalogue Gram and normal vectors. `independent3.json` and `independent4.json` record successful replays. `source_review.json` records the comparison and inherited-source scope. `SHA256.json` seals the package.

`rejection_controls.py` checks eight corruptions per dimension: bound, numerator, Gram matrix, equality group, moment tag, offset, objective, and dual multiplier. All sixteen are rejected; see `rejection_controls.json`.

The discovery solver and cached research systems are not needed for this self-contained proof package. Exact equality elimination and a rational dual certificate are reproduced directly by the independent verifier.

The source shell's minimum, population, 11-design property and catalogue identity remain explicit published inputs. We do not enumerate the billions of shell vectors or claim a new proof of those lattice facts. This is a certified existence construction with an exact selection/projection recipe.

## Sources

- G. Nebe, *An even unimodular 72-dimensional lattice of minimum 8*, J. reine angew. Math. 673 (2012), 237–247: https://arxiv.org/abs/1008.2862 .
- Gamma72 catalogue: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html .
- Catalogue Gram source: `verifications/improved/dim73-95-gamma72-caps/data/gamma72_gram.npy`, blob `5857c4e610c14b6e6e56e734e5102546f356c1dc`.
- D4 witness: `verifications/improved/dim68-69-gamma72-cross-sections/data/gamma72_d4_k4.npy`, blob `51e91c8f489ef600b4e6b8e054654c2ea06650df`.
- Current record table: https://github.com/alexlegeartis/KissingNumbers/blob/main/RESULTS.md . The preceding 70D publication and concurrent 88–94D publication do not change the compared 68D/69D values.

All computation used CPU 0, leaving reserved CPUs untouched. The package checks and targeted registry/table checks passed. The complete repository audit was not run in this environment because the historical full checkout/data were unavailable. The paper and workflows are unchanged.
