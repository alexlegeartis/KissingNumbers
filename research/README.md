# Research campaigns

These packages preserve searches, exact certificates for restricted families,
and checks of published constructions. None of these campaigns produced a new
kissing-number lower bound. The registered constructions remain in
[`../verifications/`](../verifications/) and [`../RESULTS.md`](../RESULTS.md).

| Package | Main outcome |
| --- | --- |
| [12D quaternionic campaign](dim12_quaternionic_20261104/) | An exact three-parameter 840-point family and an exhaustive ceiling for a specified whole-cell ansatz; the verified lower bound remains 841. |
| [12D paper audit](dim12_paper_audit_20261104/) | Exact certification of the published 841 decimal rays, an alternative bridge system, and a 48-system outside the paired-floor ansatz. |
| [12D contact campaign](dim12_contact_20261105/) | Exact local exclusion neighborhoods around two invalid 842-point near-misses, with discovery and recovery controls. |
| [17D and 20D campaign](compute17_20_20261104/) | An exact 19448-point ceiling for the saved finite 20D profile family; searches retained 5730 in 17D and 19448 in 20D. |
| [19D and 21D campaign](compute19_21_20261104/) | Independent reconstruction of Li's 12460- and 31560-point configurations and unsuccessful broader searches. |
| [Structured 21D campaign](attack21_structured_20261104/) | Verified recovery, nonlinear examples, exact local trade barriers, and rational bounds conditional on retaining specified holes. |
| [21D upper-bound certificate](upper21_certificate_20261104/) | An exact interval of 960–971 per component in the full line-profile family, plus an exact feasible pair-moment relaxation at 971. |

The package READMEs state each result's scope and attribution. An exclusion in a
fixed family, a solver timeout, or an invalid numerical candidate does not bound
the unrestricted kissing number.

## Replaying compact evidence

From the repository root, the following checks use only committed inputs,
Python and NumPy. They need no optimization solver, GPU, or search checkpoints:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
taskset -c 0-126 python research/upper21_certificate_20261104/verify_upper.py
taskset -c 0-126 python research/upper21_certificate_20261104/verify_relaxation.py
taskset -c 0-126 python research/attack21_structured_20261104/certify_halves.py --verify-only
taskset -c 0-126 python research/dim12_quaternionic_20261104/verify.py \
  research/dim12_contact_20261105/evidence/recovered841.npy
```

The following additional CPU checks use the pinned NumPy, SciPy, Numba and
threadpoolctl environment in
[`dim12_quaternionic_20261104/requirements.txt`](dim12_quaternionic_20261104/requirements.txt):

```sh
taskset -c 0-126 python research/compute17_20_20261104/certify_profiles20.py
taskset -c 0-126 python research/dim12_contact_20261105/local_certificate.py
```

These scripts may refresh their local result JSON files. On a host with a
different CPU layout, adapt the affinity while reserving at least one logical
CPU, as described in [`../AGENTS.md`](../AGENTS.md).

## Full campaigns and provenance

Each package documents its full replay, dependencies, seeds and budgets.
The 12D packages provide downloaders that check public input hashes before
using them. Their ignored inputs and large `runs/` directories must be restored
or regenerated for full campaign aggregation; the compact checks above do not
need them. The other campaigns similarly regenerate ignored graphs and native
executables from the saved inputs and source.

Manifests pin source and evidence bytes. Paths in result reports and numerical
run configurations record the original environment; cross-package paths in
the 21D manifests are relative to their package directories. Wall-clock search
budgets and parallel optimization do not guarantee identical trajectories.
