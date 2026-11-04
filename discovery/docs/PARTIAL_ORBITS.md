# Coupled partial sign-orbit search

`search.partial_sign_orbits` takes one independently verified explicit integer
code, exact equal-norm candidate representatives, and binary diagonal sign
generators. Its output is a plain integer-coordinate code with a fresh complete
geometry proof. No expected count or dimension-specific construction is loaded.

Each point of an orbit is selectable separately. Internally conflicting orbits
are retained, and sparse representatives with sign stabilizers are deduplicated.
Every incumbent index can be deleted, including equators, roots, height tiers,
axes and other accessory layers. The shared-blocker selector charges a deleted
point once across all insertions and layers.
For each shared blocker, candidate neighborhoods are partitioned into actual
conflict cliques. The valid inequality `sum(selected candidates in clique) <=
deleted blocker` strengthens the numerical relaxation. Every clique edge is
checked against the exact group kernel; these cuts exclude no feasible exchange.

For representatives a,b and sign messages p,q, the product is
`<g_p a, g_q b> = <a, g_(p xor q) b>` in the supplied positive diagonal metric.
A 64-bit table stores the forbidden relative phases for each representative
pair. It retains signed products; an absolute whole-orbit maximum would lose
the individual compatibility needed here. The source is checked to be closed
under the sign group, so blockers can also be transported exactly. An arbitrary
source can use the identity group, represented by an empty generator list.
Only the input source needs this symmetry; the output may break every orbit.

The backend supports at most six independent sign generators. A proved absolute
integer accumulation bound selects exact binary32 or int64 products; larger
values use Python integers. This is an acceleration of integer arithmetic,
without a rounded distance tolerance. Every chosen output pair is then checked
by the existing independent coordinate verifier.

Parameters are `representatives`, `generators`, `blocker_limit`,
`candidate_limit`, optional `pool_seed`, and `seconds`. Screening covers the supplied representative
catalogue; the two explicit limits restrict the optimizer's pool. The receipt
records orbit counts, stabilizers through actual point counts, blocker histogram,
pool truncation, selected points, shared deletions, input hashes and numerical
solver diagnostics. A timeout or gain zero is not an upper-bound certificate.
The operator never labels a result a new record.
By default the cap keeps the lowest-price orbits first. A `pool_seed` instead
selects a reproducible uniform reservoir of individual eligible points. This
explores different partial-orbit frontiers without allocating the full catalogue.

The tests include an orbit whose two points conflict with each other but whose
individual points each clear the source. Selecting one raises the control from
two to three points. Another control requires removing two source points to
insert three; no candidate is a free insertion. Other checks compare all relative
phases and transported blockers directly with integer coordinate products,
exercise stabilizers and duplicate orbits, use a weighted metric and arbitrary
precision, and reject nonclosed sources and resource violations.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  python -m unittest discovery.tests.test_partial_orbits
```

The first research application is the
[17D coupled-layer search](../../research/partial_orbits_20261104/README.md).
It differs from the earlier 16- and 32-point orbit optimizations by selecting
points individually and retaining internally conflicting candidate orbits.
