# Research campaigns

These packages preserve construction searches, exact certificates for restricted
families, and checks of published constructions. The November 4 algebraic
campaign improves 27D and 31D; the 43D packages record a separate discovery and
its reproduction by the engine. The earlier lower-dimensional campaigns below
produced no new kissing-number lower bound. Registered constructions remain in
[`../verifications/`](../verifications/) and [`../RESULTS.md`](../RESULTS.md).

## 5 November lower-bound search

- [Signed equators and coupled 13D search](signed_equator_20261105/) gives a sharp restricted-family equator proof, a source-specific saturation result, a complete fixed-alphabet scan and joint point-level searches. No new lower bound is claimed.
- [Weighted-shell and mixed-rotation 15D search](weighted_shell15_20261105/) records a complete fixed-source integer-shell scan and a checked individual-point exchange pool. No new lower bound is claimed.

## Certified construction improvements

| Package | Outcome and scope |
| --- | --- |
| [Algebraic completion](algebraic_switches_20261104/) | 209,825 in 27D (+18) and 251,196 in 31D (+8). Complete explicit coordinates, every unordered pair and independent direct-integer changed pairs; no lattice-minimum or unenumerated-population input. |
| [43D original discovery](43d_head_exchange_20261104/) | 2,919,785 points (+3) from finite body/head exchanges. Both exact representations retain the documented P48n and population inputs. |
| [High-dimensional engine campaign](highdim_characters_20261104/) | Independently reproduces all three 43D exchanges through a generic operator; also records negative 33–35D, 39D and 40D searches. Seeded discovery replay, distinct from held-out search. |

The source comparisons are dated and limited to the checked sources. Exact
geometry does not certify worldwide priority. Historical certificates and
campaigns are retained unchanged.

The later [partial-orbit continuation](partial_orbits_20261104/) adds individual
sign phases, shared deletion costs and new rational image families. Its finite
17D, 27D and 29D–31D searches have no further improvement; the large 17D
optimizations remain unresolved after their time limits.

## Lower-dimensional campaigns and controls

| Package | Main outcome |
| --- | --- |
| [12D quaternionic campaign](dim12_quaternionic_20261104/) | An exact three-parameter 840-point family and an exhaustive ceiling for a specified whole-cell ansatz; the verified lower bound remains 841. |
| [12D paper audit](dim12_paper_audit_20261104/) | Exact certification of the published 841 decimal rays, an alternative bridge system, and a 48-system outside the paired-floor ansatz. |
| [12D contact campaign](dim12_contact_20261105/) | Exact local exclusion neighborhoods around two invalid 842-point near-misses, with discovery and recovery controls. |
| [17D and 20D campaign](compute17_20_20261104/) | An exact 19448-point ceiling for the saved finite 20D profile family; searches retained 5730 in 17D and 19448 in 20D. |
| [19D and 21D campaign](compute19_21_20261104/) | Independent reconstruction of Li's 12460- and 31560-point configurations and unsuccessful broader searches. |
| [Structured 21D campaign](attack21_structured_20261104/) | Verified recovery, nonlinear examples, exact local trade barriers, and rational bounds conditional on retaining specified holes. |
| [21D upper-bound certificate](upper21_certificate_20261104/) | An exact interval of 960–971 per component in the full line-profile family, plus an exact feasible pair-moment relaxation at 971. |
| [Nonclassical 13D seeds](coupled13/) | Exact 1,154-point controls and coupled searches; no 1,155-point witness. |
| [Partial-orbit 17D campaign](mixed_tier17_20261104/) | Coupled flat-coset, tier and equator searches, including asymmetric height choices; scoped negative certificates and no new bound. |
| [Character-family campaign](character_families_20261104/) | A robust rational 19D 12,460-point control and coupled 17D/19D/21D searches, without improvement. |
| [Exact orbit controls](orbit_controls_20261104/) | Known 17D, 18D and 22–24D configurations constructed and independently replayed through the engine. |
| [Unified discovery campaign](discovery_campaign_20261104/) | Archived group/ring/support searches and source reviews; no new bound in that campaign. |
| [Construction-search continuation](sol_discovery_20261104/) | Separate 1–24 held-out measurement and finite 13D/14D/31D search pools; no improvement in this earlier continuation. |

The integrated current-source controls cover 19 of the 24 target dimensions.
The separate held-out catalogue matches ten targets using formula controls
and support search; actual support search matches 4D, 7D, 8D and 9D. Replaying
a known configuration is not credited as discovering a new record. See the
[coverage matrix](discovery_campaign_20261104/coverage_1_24/) and
[construction-search contract](../discovery/docs/CONSTRUCTION_SEARCH.md).

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
