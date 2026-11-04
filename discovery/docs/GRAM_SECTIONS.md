# Exact rational Gram codes and weighted central sections

These are general checked operators, not a dimension-to-record dispatch table.
Neither operator claims novelty or optimality. Both replay their mathematical
conditions from artifact data and use explicit exceptions under Python `-O`.

## `geometry.gram_code`

Parameters are `gram`, a nonempty square exact rational matrix, and optionally
`norm_squared`. A rational JSON scalar is an integer or `[numerator, denominator]`
with an integer numerator and a strictly positive integer denominator. Floats,
Boolean values, string approximations and zero/negative denominators are rejected.
Raw integer bit lengths are checked before Fraction construction.

The operator returns `gram_spherical_code`, with the matrix as its symbolic
geometry descriptor. Its independent verifier in `discovery/exact_gram.py` checks:

1. Exact symmetry and a common strictly positive diagonal, matching the optional
   supplied norm.
2. Every off-diagonal pair: duplicates are rejected and each inner product must
   be at most half the common squared norm.
3. Positive semidefiniteness and exact rank, using symmetric rational Schur
   elimination. Negative Schur diagonals are rejected. If no positive diagonal
   remains, **every** residual entry must vanish; a zero diagonal alone does not
   establish semidefiniteness.

The exact number of positive pivots is `dimension = rank_exact`. A PSD Gram
matrix has a real realization in that rank, so `ambient_dimension` is set to
this sufficient dimension and `dimension_scope` explicitly says symbolic real
realization. Rational Cartesian coordinates or algebraic coordinate extraction
are not asserted or implemented.

Facts are recomputed. Serialized verification flags, claimed ranks and supplied
cardinalities are never mathematical authority. CLI comparison and scheduler
ranking accept this new exact geometry kind. The default comparison catalogue
still contains only its existing dated dimensions; the separate 1–24 coverage
matrix is broader. A valid 5D receipt therefore needs an appropriate supplied
catalogue for `compare`, or produces the honest missing-baseline error.

`discovery/recipes/gram5_control.json` exercises the general interface on the
known 40-point 5D code recovered in the earlier numerical search. Its exact
result is rank 5, norm squared 2, maximum distinct inner product 1, and 780
checked pairs. This reproduces a known bound and is not autonomous discovery
by the Gram import operator.

## `geometry.section`

This operator takes one explicit exactly verified integer spherical code and a
nonzero rational `normal` vector matching its ambient coordinate width. If rows
are x and the positive diagonal metric has weights w, it selects **all and only**
rows satisfying

    sum_j w[j] * normal[j] * x[j] == 0.

Thus it takes a central hyperplane section in the existing metric. There is no
rounding, coordinate puncturing, weight change, or unproved dimension decrement.
The output retains the integer coordinates and their positive diagonal weights.
Actual span rank is recomputed exactly; if the source already lies in the
hyperplane, its rank can stay the same. Empty sections and zero normals are
rejected.

The ordinary `spherical_code` output embeds the minimal explicit source rows,
metric and normal in its `section` descriptor. Independent replay rechecks the
source geometry, regenerates the complete selected index set, compares the
stored output rows and metric, and rechecks output norms, distinctness, all
pairs and exact rank. Chaining sections stores the immediate source coordinates
rather than an unbounded recursive ancestry proof.

`discovery/recipes/section_e8_chain.json` reconstructs E8 with the general root
operator, then cuts the hyperplanes x6+x7=0 and x5+x7=0. It obtains the known
126-point E7 code of exact rank 7 and the 72-point E6 code of exact rank 6,
while preserving all eight ambient coordinates and the original norm squared 8.
The weights in the tests are deliberately nontrivial so an unweighted dot-product
implementation would select the wrong set.

## Exact span and resource limits

The explicit integer coordinate verifier now uses bounded fraction-free row
elimination in place of an unrestricted symbolic rank call. It scales the
**whole** candidate row during elimination, supports out-of-order nonunit pivots,
and divides by a common integer gcd without changing the rational row span.
Positive diagonal weights form an invertible real transformation, so this
integer row rank equals the rank of the corresponding weighted Gram and the
actual Euclidean span.

Existing `max_materialized_points`, `max_dimension` and `max_pairs` limits remain.
Three additional overridable `RunContext` budgets are:

| Field | Default | Meaning |
|---|---:|---|
| `max_gram_entries` | 2,000,000 | Dense Gram scalar count, checked before rational conversion |
| `max_rank_work` | 50,000,000 | Counted exact elimination scalar updates, checked before each update batch |
| `max_rational_bits` | 4,096 | Raw integer numerator/denominator and stored intermediate bit lengths |

Section source row, width and pair limits precede normal Fraction allocation.
The Gram verifier bounds matrix size and pair count before allocating rational
entries. Rank dimension limits precede the next Schur update. Per-code limits
are replayed for the source and output; this is not a global cumulative
arithmetic allowance across a whole construction DAG. Scheduler subprocess
wall-time limits provide the separate execution deadline.

The new tests include indefinite and singular-invalid Grams, duplicate points,
wrong norms/ranks, forged facts, overbudget preallocation traps, rational metric
sections, altered membership/normal/weights, and Python `-O`. Adversarial rank
fixtures reproduce the partial-row-scaling defect caught during review; random
and deliberately dependent integer matrices are compared to SymPy in tests.
The independent reviewer additionally compared 500 diverse/dependent rank cases
and 120 generated PSD Gram ranks and cleared the corrected implementation.

From the repository root:

    OPENBLAS_NUM_THREADS=1 python -m unittest discovery.tests.test_gram -v
    OPENBLAS_NUM_THREADS=1 python -O -m unittest discovery.tests.test_gram -v
    python -m discovery run discovery/recipes/gram5_control.json --out /tmp/gram5-new.json
    python -m discovery verify /tmp/gram5-new.json
    python -m discovery run discovery/recipes/section_e8_chain.json --out /tmp/sections-new.json
    python -m discovery verify /tmp/sections-new.json --output e6

Receipt paths must be new because the CLI refuses to overwrite previous runs.
