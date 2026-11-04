# Independent pre-handoff review

A fresh Sol agent that did not implement the core reviewed the final package,
including the two native DGLS parent proofs, graph composition, exact/numerical
status flow, integer overflow protection, implicit dimensions, resource guards,
and adversarial modifications to certificates and receipts.

Three medium-priority findings were fixed and independently rechecked:

1. **Unvalidated auxiliary claims.** Extra geometric bounds, dimension-scope
   assertions and novelty flags could survive a successful coordinate replay.
   The shared verifier now rejects protected mathematical fields it cannot
   establish and refuses affirmative SOTA or unrestricted-optimality claims.
   Computed cardinalities and dimensions were not incorrect in the original
   cases, but the extra assertions were not certified.
2. **Late numerical allocation guard.** Numerical coordinates were converted to
   a NumPy array before checking their dimension. JSON shape, width and finite
   scalar checks now run before allocation. The regression test patches the
   allocator to fail if it is called for an invalid input.
3. **Negative-proof replay under optimized Python and empty coverage.** Assertions
   could be disabled, and an empty report could claim success without checking
   the two modules. Replay now uses explicit exceptions, requires both factor
   indices and nonempty witnesses, checks full binary rank to establish ring
   freeness, and rejects surviving extensions. The original malformed cases are
   rejected both normally and under `python -O`.

The reviewer independently confirmed the fixes and found no remaining blocker,
false spherical lower-bound certificate, or graph path promoting numerical
coordinates or an undistanced ring structure into an exact geometric result.
The saved 96D certificate replays with 2^37 points, binary distance at least 24,
and maximum unnormalized inner product at most 48.

## Final-tree checks

- Full package suite: **55 tests passed**, 22.861 seconds.
- Optimized core/coset/negative-replay subset: **25 tests passed**, 0.330 seconds.
- Valid negative replay under optimized Python: both 262144-extension families
  excluded with 11 and 10 witnesses, no survivors.
- All four end-to-end controls replayed exactly after the fixes: 10D 510,
  14D 1932, 32D 352800 and the implicit 96D 2^37 sign layer.
- Repository audit and archival round-trip verification are recorded in the
  handoff validation. No new SOTA or global-optimality claim is made.

Receipt hashes detect modification; they do not authenticate historical
provenance or establish novelty. The work remains uncommitted and unpublished.


Later integration evidence is recorded in `research/discovery_campaign_20261104/FINAL_STATUS.md`, `REVIEW.md`, and `discovery/benchmarks/planner_integration_20261104.json`. The earlier measurements above remain historical snapshots.
