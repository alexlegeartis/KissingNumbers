# Gamma72 profile models and independent 69D proof

This directory contains the root lane's affine-profile experiments, an
independent derivation of the successful 69D population certificate, and a
bounded follow-up with finer head heights. The maintained exact proof is in
[dim69-gamma72-head-compression](../../../verifications/improved/dim69-gamma72-head-compression/).

## Outcomes

| Experiment | Recorded value | Comparison | Status |
|---|---:|---:|---|
| D4 affine projection to 68D | approximately 471,884,288 | Existing 68D bound 471,884,288 | Rediscovery of the existing value |
| D5 affine projection to 67D | approximately 185,406,464 | Existing 67D bound 355,090,822 | Below the record; no exact certificate exported |
| Independent D4 one-head population proof in 69D | **743,522,304 exactly** | Previous 69D bound 689,452,922 | Exact dual exported and independently verified in the final package |
| Finer one-head grid in 69D | approximately 735,182,848 | New exact bound 743,522,304 | Below the new record; numerical research only |

The affine candidate JSON files explicitly have `exact_certificate: false`.
Their displayed numerical values are not promoted to new rigorous bounds.
The existing 68D certificate, reproduced numerically here, remains available
in the earlier Weyl-averaged affine package.

## The profile model

`dk_affine.py` works in orthonormal coordinates on a scaled Dk normal space.
The removed component is h/4, where h has integer coordinates that are either
all even or all odd. Apart from the known normal roots, the two largest absolute
coordinates of h sum to at most eight. The known roots have two coordinates
equal to plus or minus eight and all others zero.

Signed-coordinate permutation averaging identifies profile populations. The
D4 model also uses the half-Hadamard triality map. Minimum-vector translation
identities give further valid identifications. This is averaging over actual
normal orientations, not an assertion that Gamma72's shell has those symmetries.
The final proof separately verifies the normal embedding in the catalogue Gram.

The model reconstructs symmetric even moments through degree ten from the
spherical 11-design formula and subtracts the known normal roots. After
redundant equations are removed, the D4 model has sixteen variables and ten
rows; D5 has thirty-nine variables and eighteen rows. Numerical linear programs
bound selected affine unions, while bounded mixed-integer searches select
compatible offsets. Those research scripts do not themselves certify their
floating-point optima.

## Independent exact 69D derivation

`independent69_model.json` records the ten orthonormal profiles and rational
heights corresponding to the construction discovered in `../theory/`.
`extract69.py` rebuilds the D4 moment model, recovers a rational dual from the
active columns, and checks every dual slack exactly. It writes
`independent69_certificate.json`, whose objective is exactly 743,522,304 and
whose common dual denominator is 1,451,520.

The final package's standard-library `verify_orthonormal.py` reconstructs the
domain, groups, rows, embedding, dual and geometry without importing this
research code, an optimizer, NumPy, or SymPy. Its proof is independent of the
Gram-coordinate verifier and of the numerical discovery trajectory.

## Finer head grid

`fine_heads69.py` reuses the full D4 Gram-profile population model, but searches
only profiles of squared projection norm at most three quarters. Heights use
a grid of spacing 1/128 for norms at most one half and spacing 1/32 for the
remaining retained profiles, with the two exact-radius endpoints also sampled
numerically. A weighted clique heuristic chooses compatible profile/height
pairs, then a separate population minimization evaluates the selected objective.

The search score was approximately 741,459,957; the numerical worst-case
population objective was approximately 735,182,848. Neither exceeds the new
certified 743,522,304 value. In particular, this narrower norm cutoff excludes
the two norm-one profiles used by the successful certificate. It is not an
exhaustion of the head-compression family.

## Reproduction and file roles

For exact verification, run from the repository root:

```sh
python verifications/improved/dim69-gamma72-head-compression/verify_orthonormal.py
```

To reproduce the research computations, install NumPy, SciPy and SymPy, then
run the following from the repository root. Apply host-appropriate CPU affinity
and leave one logical CPU free; the thread limits below bound numerical-library
worker pools.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/gamma/dk_affine.py 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/gamma/dk_affine.py 5
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/gamma/extract69.py
g++ -O3 -std=c++17 research/dim32plus_20261113/theory/head_search.cpp -o research/dim32plus_20261113/theory/head_search
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/gamma/fine_heads69.py
```

The searches use time limits, so reproducing their complete random trajectories
is not required to replay a saved exact certificate. `model_d*.json` records
the profile systems; `candidate_d*.json` and `fine_heads69_result.json` record
numerical selections; the `.log` files preserve search traces. The generated
`fine_heads69.txt` vertex pool is reproducible and excluded from version control.
**Raw pools, logs, numerical multipliers and search scores are research material,
not lower-bound certificates.** No search remains running.
