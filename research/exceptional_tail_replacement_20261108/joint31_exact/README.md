# Exact exclusion of the final 31D continuation pool

8 November 2026. **No exchange from this explicit 3,762-point candidate pool can
increase the frozen 251,200-point configuration.** Arbitrary incumbent-point
deletions are allowed. This is a finite-pool exclusion, not an upper bound for
the kissing number and not a claim that all construction families were searched.

This directory is self-contained. Its incumbent, candidates, constraints,
rank certificates and rational dual are in `data/`. The geometric certificate
for the 251,200-point construction remains independent of this selection result.

## Replay

Requirements: Python, NumPy, SciPy and SymPy. The checker uses sparse matrix
arithmetic from SciPy but **does not invoke an optimizer**. It needs no network,
external files, previous PASS receipt, search state, or parent checkout.

```sh
OPENBLAS_NUM_THREADS=2 python verify.py
python negative_controls.py
```

The replay takes a few seconds and prints `"status": "PASS"`. It does not modify
any input. Reserve one available logical CPU as required by the repository's
AGENTS.md when running concurrent computations. On the development host the
replay was pinned to CPUs 4–5, leaving CPU 8 unused by research jobs.

## Candidate pool

The frozen incumbent has 251,200 vectors of squared norm 128. The compatible
inner-product threshold is 64. The pool combines:

- 3,008 tetrad-body candidates from 104 completely priced tail directions out
  of the 574-direction alphabet. This part of the catalogue is partial.
- 750 candidates from the complete signed six-support search with its stated
  two-blocker cutoff, before the four-for-eight improvement.
- The four points removed in that improvement, eligible to return.

All 3,762 points are distinct and absent from the frozen incumbent. Every blocker
was recomputed against the new incumbent; old row indices were not reused.
There are no free insertions. The finite pool is not enlarged by this proof.

## Constraint certificate

Let x_i indicate a selected candidate and d_j a deleted incumbent point. Of the
251,200 incumbent points, 2,650 conflict with at least one candidate. Every other
incumbent point clears the entire pool, as the replay checks. Deleting such a
point is unnecessary: restoring it preserves compatibility and increases size.
Thus allowing just the 2,650 relevant deletion variables loses no improving
exchange, even when arbitrary deletions were originally allowed.

The certificate contains 20,717 valid inequalities:

- 10,665 blocker-clique inequalities. A mutually conflicting candidate clique
  C whose members all conflict with old point j satisfies sum(x_i, i in C)≤d_j.
- 10,000 candidate-clique inequalities sum(x_i, i in C)≤1.
- 52 fixed-tail capacity inequalities. For a tail q of squared norm 64, its
  core vectors have squared norm 64 and must have pairwise nonpositive products.
  In a real subspace of dimension r, Rankin's bound gives at most 2r such points.
  The default r=24 gives 48. For 48 fibers, exact integer null-space certificates
  prove smaller rank bounds on the union of all old and candidate core vectors.

Every clique product and null-space equality is checked directly with integers;
null-space independence is verified by a nonzero exact Gram determinant. The
sparse coefficient matrix is reconstructed from these geometric constraints and
compared with the stored matrix. No optimizer status is treated as evidence.

The replay of non-variable old-point coverage uses binary32 matrix products of
integer coordinates of magnitude at most 8. Every product and partial sum has
absolute value at most 31·8²=1,984<2²⁴, so these products are exact, not tests with
an angular tolerance. Final norms, distinctness, and hashes are checked too.

## Rational dual bound

The minimized integer objective is negative gain:

    c·z = sum(d_j) − sum(x_i).

Write the constraints as A z≤b and 0≤z≤1. The certificate stores nonpositive
integer vectors y and u and denominator D=1,000,000,000. Let

    r = D c − Aᵀ y − u.

For every feasible z, multiplying A z≤b by y≤0 and using the variable bounds gives

    D c·z ≥ b·y + sum(u) + sum(min(r_i,0)) = −428571624.

All quantities on the right are recomputed as integers. Therefore

    negative gain ≥ −428571624/1000000000 > −1.

Negative gain is integral, so it is at least zero. The empty exchange attains
zero. This proves the stated exact finite-pool optimum, even though the floating
LP objective used to discover the dual was about −3/7.

## Files and provenance

- `data/incumbent251200.npz`: complete frozen incumbent, SHA-256
  `b75a8e9d958017fe03926265f5ee9d42f24b16c9ec848414b4a4e7eba0d329ee`.
- `data/repriced_pool.npz`: the explicit candidate pool and complete blocker data,
  SHA-256 `58bcda47bb38581fcd99e6dab809ac1473e68280b1712289918efa89759dd4b5`.
- `data/lp_constraints.json`: geometric meaning of every matrix row.
- `data/rank_certificates.json`: exact integer null spaces for the smaller fibers.
- `data/lp_dual_certificate.npz`: matrix, objective-variable ordering and rational
  dual numerators. No floating dual values are required by the checker.
- `data/report.json`: historical numerical selection report and source-pool
  provenance. Its repository paths are provenance only, not replay dependencies.
- `data/independent_exact_exclusion.json`: the original exact replay receipt,
  retained for comparison only; `verify.py` computes a new result.

The candidate construction and geometric improvement were found in
`research/exceptional_tail_replacement_20261108`. The combined repricing,
clique/rank model, rational dual and optimizer-free replay were developed in
`research/midrange_descent_20261108`. These are research provenance locations;
all required inputs have been copied into this directory.
