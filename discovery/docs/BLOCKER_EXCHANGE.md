# Shared-blocker selection

`selection.blocker_exchange` assembles candidates supplied by a geometric or
algebraic compatibility oracle. Each candidate lists the incumbent indices
that must be deleted. Candidate conflicts and optional clique inequalities
control new/new compatibility. No incumbent index is protected.

Deleting a shared blocker is charged once. Optional positive integer
`candidate_weights` and `removal_costs` support complete families with unequal
populations. A cost-two blocker shared by three weight-one candidates gives
a profitable three-for-one-family exchange even though each candidate alone
loses weight. Point-level choices can break an otherwise imposed orbit.

The numerical MILP proposes a selection. Rounded bits are checked against
every supplied exact incidence, conflict and clique; extra deletions are
removed. A compatible free pool is preferred to neutral changes at equal
gain. Numerical dual bounds remain diagnostics. Resource and integer-size
limits apply before allocating the optimization matrix.

The result is an `exchange_proposal` with status `unverified`. Supplied
geometric oracles, family populations and source validity remain assumptions.
Independent coordinate replay is required to establish a lower bound.

The 27D/31D programs in
[`research/algebraic_switches_20261104`](../../research/algebraic_switches_20261104/README.md)
use exact signed-support oracles and full changed-point products. Their final
complete coordinate proofs certify **209,825** and **251,196** points.
