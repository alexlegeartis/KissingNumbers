# Clifford obstruction and continuous cyclic orbit search, 2026-11-04

No new kissing-number record was found. All files in this directory are new;
existing discovery code and prior verification archives were only read.
No package was installed, no push or workflow action was performed.

## Exact finite-geometry result

Let S(m,r) contain every normalized real vector supported on an affine r-flat
of F_2^m, with signs (-1)^q for every Boolean polynomial q of degree at most two.
Both global signs are included. The exact size of this layer is

    [m choose r]_2 * 2^(m-r) * 2^(1+r+r(r-1)/2).

The union over even r and the union over odd r are each kissing configurations.
For distinct supports of equal dimension their intersection has at most half
the support size. For support dimensions differing by at least two the normalized
inner product is at most 1/2. On the same support, a nonconstant nonzero quadratic
Boolean polynomial has bias at most 1/2. The constant-one difference gives -1.
These observations prove the within-sector inequalities without floating point.

Between the two sectors, a conflict (inner product >1/2) occurs exactly when one
affine support is a hyperplane of the other and all signs agree on the smaller
support. For a vector of support dimension r the numbers of neighbors below
and above it are

    2*(2^r-1) and 2^(r+1)*(2^(m-r)-1),

respectively. Their sum is 2*(2^m-1), independently of r. The two sectors have
equal size M. Consequently the conflict graph is a regular bipartite graph.
For any independent subset A union B, every edge out of A ends in the complement
of B. Counting edges gives degree*|A| <= degree*(M-|B|), hence |A|+|B| <= M.
One whole sector attains M. This proves the maximum for this entire family,
including arbitrary unions of fragments of its two Clifford orbits.

| Dimension | All real stabilizer states | Exact maximum kissing subset | Verified repository baseline |
|---:|---:|---:|---:|
| 16 | 8,640 | 4,320 | 4,320 |
| 32 | 293,760 | 146,880 | 352,800 |
| 64 | 19,388,160 | 9,694,080 | 335,496,064 |

`clifford_stabilizer.py` independently enumerates all 293,760 states in 32D,
checks distinctness, verifies the six support-dimension representatives against
the full family with integer inequalities, and verifies every degree formula.
Affine coordinate permutations and quadratic sign multipliers are transitive
on each support-dimension layer and preserve the whole family, so these six
checks cover every pair type. Optimized Python mode (-O/-OO) is explicitly
refused before any check or proof output because assertions must remain active. Runtime was 1.58 seconds, peak RSS 60,224 KiB.

`bw32_reconstructed_control.npz` contains 146,880 integer rows from one sector.
Normalize row x by 2^(support_dimension/2). This is a reconstruction of the known
Barnes–Wall control, not a new record. The report calls out that it is not a
separate exhaustive all-pairs certificate. There is no numerical new candidate
being promoted from this control.

Any coordinate section in dimensions 25–31 is also noncompetitive: even the
pool omitting a single coordinate has only 142,290 states. All single coordinates
are equivalent under affine permutations; further removals only reduce the pool.
The additional 32 sampled sections per dimension in `results.json` describe
those particular coordinate sets, not an exhaustive optimization over sections.

This is a bound for real stabilizer states only. It does not bound arbitrary
Clifford orbits, other group representations, deformed vectors, or kissing
numbers. It avoids repeating Leech-shell/backbone extension searches.

## Flexible continuous families

`cyclic_continuous.py` searches independent seeds in real Fourier representations
of cyclic groups, first retaining antipodality. `cyclic_asymmetric.py` removes
antipodality, optimizes every seed's amplitude and phase independently, permits
repeated real irreducibles, and optionally adds a one-sided fixed pole.
`flexible_scan.py` goes to the larger continuous multiplicity spaces of C2 and
C3, including asymmetric trials in dimensions 5,6,7,9,10. The C2 representations
can have multiple trivial coordinates; C3 can act by the same rotation on many
independent two-planes. These are not discrete eigenline searches or a fixed
list of shell coordinates. All starting seeds are randomized; no known record
coordinates are read.

The objective is the sum of squared positive violations above cosine 1/2,
with analytic gradients, row normalization, and L-BFGS-B. Orbit compression
reduces the objective to seed versus full-orbit comparisons. Independent full
pair matrices screen every stored trial. Finite-difference gradient checks
passed. A small positive loss is not a certificate, and failed arrays are
explicitly labeled as failed point sets.

The initial nonantipodal scans tried 41/42 points in 5D, 75/77 in 6D, and 127 in
7D, for 160/160/160/160/120 starts. None was admissible: their best maximum
cosines were 0.523905, 0.543211, 0.580172, 0.553351, and 0.519956. These are
bounded numerical failures, not nonexistence theorems. The 40-point 5D control
was approached within 7.9e-9 in cosine. `refine_control.py` then polished this
independently generated control to maximum cosine 0.5000000000000021
(39 least-squares evaluations, 0.039 seconds). This is a numerical known-bound
reproduction. `certify_control.py` then proposed an integer Gram matrix by
rounding twice the numerical Gram and verified every claim in exact Fraction
arithmetic: 40 points, rank 5, norm squared 2, maximum distinct inner product 1,
and positive semidefiniteness via H=A*B^-1*A^T with B positive definite.
All 780 pair inequalities pass; norm and angle mutations are rejected.
`d5_40_exact_gram.json` and `d5_40_exact_receipt.json` are the exact control
certificate and receipt. Rounding proposes the data and supplies none of the proof.
This exactly reproduces the known 40 bound, not a new SOTA. An initial 242-point 8D diagnostic was infeasible as expected
from the known exact 240 bound; no improvement is claimed in that dimension.

The multiplicity-space pass completed 3,770 starts in 350.35 seconds, with
peak RSS 66,328 KiB. Its per-job starts were 1000/1000/600/600/245/300/17/8.
Final multiplicity-space results and all restart details are in
`flexible_results.json`. Across all three passes the total was 4,640 starts
in 367.33 seconds. The final 10D 513-point trial failed with a near-collision
(maximum cosine rounded to 1); it is explicitly a failed array. The 10D job
used 70.95 seconds against its 55-second target because wall budgets are checked
between restarts; its last individual optimization exceeded that target.
No search process remains running. Independent full-pair replay checked all
20 saved numerical point sets; only the known 40-point 5D control passed its
1e-9 numerical tolerance. `summary.json` records this final census. Saved NPZ
arrays with violations are trial coordinates, not kissing configurations.

## Reproduce

Existing Python 3.12, NumPy 2.3.5 and SciPy suffice. No downloaded matrix catalogue
or additional package is needed. From this directory run:

    OPENBLAS_NUM_THREADS=1 python clifford_stabilizer.py
    OPENBLAS_NUM_THREADS=1 python cyclic_continuous.py
    OPENBLAS_NUM_THREADS=1 python cyclic_asymmetric.py
    OPENBLAS_NUM_THREADS=1 python flexible_scan.py
    OPENBLAS_NUM_THREADS=1 python refine_control.py
    OPENBLAS_NUM_THREADS=1 python certify_control.py
    OPENBLAS_NUM_THREADS=1 python replay_candidates.py

Wall budgets are checked between restarts; each optimization also has a
3,000-iteration ceiling (1,800 in the initial antipodal pass). Thus a single
restart can slightly exceed its wall target. The scripts report wall time, peak RSS, RNG seeds, full-pair numerical cosine,
and exact-versus-numerical status. Floating optimization may vary by platform;
the finite-geometry proof and integer counts do not.

## Primary sources and comparisons

- Nebe, Rains, Sloane, *A Simple Construction for the Barnes-Wall Lattices*,
  https://arxiv.org/abs/math/0207186 (Clifford/Barnes–Wall construction).
- Nebe, Rains, Sloane, *The Invariants of the Clifford Groups*,
  https://arxiv.org/abs/math/0001038 (real Clifford group and its lattice action).
- Kliuchnikov, Schönnenbeck, *Stabilizer operators and Barnes-Wall lattices*,
  https://arxiv.org/html/2404.17677v2 (stabilizer/Clifford connection; read live
  2026-11-04). The specific finite-geometry graph bound above is derived here,
  not attributed to these sources as a new theorem.
- Repository `README.md` verified-baseline table, read 2026-11-04: 5D40, 6D72,
  7D126, 8D240, 9D306, 10D510, 16D4320, 32D352800, 64D335496064. The prior live
  table snapshot is `research/group_search_20261104/baselines.json`.

No external approval or access blocker occurred. The remaining blocker is
mathematical: these tested families have not produced a better admissible code.
