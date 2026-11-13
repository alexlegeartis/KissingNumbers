# 26D–27D global fiber search

This branch adds two search engines for the remaining deficits against
PackingStar's 11 November 2026 release. It does **not** claim a new record:
the target coordinate archives are not mounted in this execution environment,
so only the small helper/catalogue tests have been run here. Run the
data-dependent searches on the research machine with the adjacent
`PackingStar` checkout.

## 26D: component-wise optimization of the 759-octad pool

`solve_joint26.py` consumes `continuous26/transplant759.npz`, the 2,004
candidate points collected from all 759 Golay octads. The previous
`joint26_pool.py` replayed blockers exactly but sent one monolithic MILP to
HiGHS with presolve disabled; the process hit its 160-second wall limit before
returning a selection. That failure was inconclusive.

The replacement:
1. Replays every candidate's complete blocker list exactly in Q(√3).
2. Rebuilds all candidate–candidate conflicts using int64 dot products.
3. Decomposes the combined candidate/blocker/conflict incidence graph.
4. Generates 34 greedy feasible selections per component, exploiting shared
   blocker reuse.
5. Solves components separately with presolve enabled and zero requested MIP
   gap, retaining the better feasible result from the solver or greedy search.
6. Replays all selected candidate conflicts and all changed source/candidate
   pairs exactly, then writes a coordinate proposal.

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
python research/reclaim_26_31_20261113/solve_joint26.py \
  --source ../PackingStar/25D-31D/new_bounds_20261111/D26_201180_coordinates.npz \
  --pool research/reclaim_26_31_20261113/continuous26/transplant759.npz \
  --out /tmp/joint26_components --seconds 180
```

The output `proposal_coordinates.npz` is a full coordinate witness over
Q(√3), with changed-pair replay. Unchanged source-source pairs inherit the
frozen source certificate. A positive gain still requires the repository's
independent full final verifier before any leaderboard update. Timeouts and
finite-pool results are not global exclusion proofs.

## 27D: full E8-head alphabet with partial D3 fibers

`search_partial_fibers27.py` abandons the three preselected head frames used
in the earlier 438-point E8/D3 transplantation pilot. It enumerates all 2,160
integral norm-16 heads in the 8D E8 coordinate block and all 12 norm-16 D3
tails `sqrt(2)*(±2,±2,0)` and coordinate permutations: 25,920 candidates
before exterior filtering.

For each candidate, products with the 2,807 distinct projected incumbent
profiles are `p+√2 q`. The predicate `p+√2 q>16` is evaluated exactly with
sign-aware integer squaring. The true multiplicity of each projected profile
is charged as its deletion cost. Candidate–candidate products are
`h·h' + 2 t·t'`, also exact integers. The binary optimization permits
arbitrary partial fibers and shared incumbent deletions.

Run:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
python research/reclaim_26_31_20261113/search_partial_fibers27.py \
  --source ../PackingStar/25D-31D/new_bounds_20261111/D27_209862_coordinates.npz \
  --out /tmp/partial_fibers27 --cutoff 12 --max-candidates 12000 \
  --seconds 180
```

The blocker cutoff and pool cap define an explicit finite subproblem, not an
exhaustive search. Increase them in staged runs. A positive output is encoded
as `rational + sqrt(2)*radical`; it includes exact blocker and candidate-pair
checks, but must still pass an independent full-configuration verifier.

## 27D extension: rational spherical tails

The earlier D3-tail catalogue is still useful as a small controlled baseline,
but its twelve tail directions are unnecessarily restrictive. The new
`search_rational_sphere_tails27.py` replaces that alphabet with rational
stereographic points on the full radius-4 sphere in the three tail coordinates.
For coprime integers `(a,b,c)`, `c>0`, put `d=a^2+b^2+c^2` and

```text
t = (8ac, 8bc, 4(c^2-a^2-b^2))/d.
```

Then `||t||^2=16`. Pair it with any of the 2,160 norm-16 E8 heads. Unlike
the D3 pilot, this varies the tail continuously over a rationally dense
subset of the sphere. Candidate/source products are compared after
multiplication by `d`; candidate/candidate products after multiplication by
both denominators. No floating-point geometry is used.

Run:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
python research/reclaim_26_31_20261113/search_rational_sphere_tails27.py \
  --source ../PackingStar/25D-31D/new_bounds_20261111/D27_209862_coordinates.npz \
  --out /tmp/rational_sphere27 --K 2 --blocker-cost 512 \
  --max-candidates 10000 --seconds 180
```

Increase `K` and the blocker-cost threshold in staged runs. The old D3
engine's default blocker cost was corrected from 12 to 256: projected source
profiles often have multiplicity 32 or more, so a cutoff of 12 would
artificially discard every candidate that blocks even one such profile.
Neither finite rational search excludes algebraic tail points outside its
parameter catalogue.

## Tests

Run:

```bash
python -m unittest -v research/reclaim_26_31_20261113/test_reclaim_searches.py
```

Tests cover Q(√3) and Q(√2) sign predicates, the 2,160/12 point catalogues,
exact norms, connected-component construction, and shared-blocker greedy
selection. They do not replace the data-dependent coordinate runs.
