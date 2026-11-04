# Current capability inventory — 2026-11-04

No new certified SOTA has been found. The reviewed foundation and finite scheduler are pushed at e25c2924f377a1061b6c60429c950d34 df 7515e9. The subsequent integrations below remain uncommitted and have not been pushed.

## Executable operators

- `codes.translate_union`
- `cosets.build`
- `cosets.union`
- `geometry.gram_code`
- `geometry.integer_code`
- `geometry.orthogonal_union`
- `geometry.section`
- `geometry.sign_embed`
- `gluing.union`
- `groups.from_generators`
- `groups.ganzhinov10`
- `groups.ganzhinov14`
- `groups.real_stabilizer_sector`
- `groups.search`
- `groups.weil`
- `lattice.root_system`
- `layers.known32`
- `layers.signed_supports`
- `rings.binary_factor`
- `rings.extend`
- `rings.gray`
- `rings.hensel`
- `rings.verify_dgls`
- `search.coset_assembly`
- `selection.coupled_union`
- `selection.isometric_exchange`
- `signed_layers.conditional_swap`
- `signed_layers.coupled_lift`
- `signed_layers.free_zone_join`

All 29 registered operators have executable typed-planner signatures. Detailed new operator-contract schemas/resource estimators remain explicitly partial (3 operators); signature support is not a claim of full mathematical contract coverage. The coset assembly operator performs actual candidate generation but returns UNVERIFIED search results and conditional cardinality forecasts, never a certified sphere.

## Shared infrastructure and proof scope

The allow-listed typed DAG remains intact. Added: finite program enumeration, checkpointed typed DAG synthesis, versioned mathematical-object schemas, exact rational Gram/rank certificates, exact central sections, generic root systems and symbolic real-stabilizer sectors, finite heterogeneous nonorthogonal gluing, point-level isometric exchange/coupled union, and signed-layer switching/free-zone/lift constructions with exact radical comparisons. Numerical coordinates and sampled-distance upper estimates remain candidates only.

Gluing currently needs explicit coordinate/Gram objects and supplied cross-frame Gram data; joint PSD checks their consistency. It supports nonzero rational scales, rational heads, exact overlap deduplication and computed output rank. It does not infer unknown relative frames or handle arbitrary gigantic implicit code/orbit unions. Signed-layer raw coordinates use Q(sqrt 3); individual normalization can introduce further radicals, handled by exact norm-aware comparisons rather than claiming all normalized coordinates stay in that field.

## Actual controls and performance

- Known 10D 510 and 14D 1932 exact controls; the 14D frame construction is formula replay, not autonomous nonabelian representation discovery.
- Generated root controls 1–8, exact 5D Gram reconstruction, E8 central sections, and symbolic 16D 4320 control.
- Full 32D 352800 compressed control remains implemented.
- Full 96D 142162362560 replayed by both existing independent checkers, including all five layers. `run_full96.py` is a specialized benchmark bridge, not generic layer synthesis.
- Native 28D Weil search still yields 6720 numerical points versus 221988. Duplicate phase-block budget waste was fixed, but no competitive improvement resulted.
- Fresh competitive 28D multi-image pools started at 221988 and did not improve it. New 13D/14D coupled-layer experiments retained 1154/1932. New 96D subcode/coset candidates were rejected by stronger searches. No new record is concealed in these controls.

Final integration tests reported 168 discovery tests passing. Coverage is tracked by the separate executable 24-row matrix and fresh-source snapshots; imported/formula-generated/search-discovered paths are distinct. See the individual campaign READMEs for tested versus timed-out and unsearched scope.

## Explicit scaling limits and missing capabilities

Default budgets include 20000 materialized points,10000 orbit points,10000000 pairs,262144 codewords,200 trials,30 solver seconds and 100 DAG nodes. Dimension cap 256 is an input guard, not a tractability guarantee. Binary syndrome backend is limited to 16 syndrome bits and 2^24 transitions. Gram/rank/bit-size and search compiler limits are explicit. Generic orbit searches still materialize bounded numerical orbits; unclosed orbits are not valid whole-orbit candidates and are not promoted.

Missing or restricted: exact generic group-orbit reconstruction/double-coset proofs, nonfree ring modules and broad outer-code composition, interval/algebraic backends beyond implemented cases, arbitrary implicit cross-family gluing and overlap counting, generic high-dimensional shell/height synthesis, automatic continuous family optimization in the planner, and complete contract/estimator metadata. No engine can promise finite-budget completeness over all constructions or many new records.

The overriding next milestone is a candidate exceeding a freshly verified current lower bound, followed by independent exact or rigorous certification—not another known control or another failed finite search.
