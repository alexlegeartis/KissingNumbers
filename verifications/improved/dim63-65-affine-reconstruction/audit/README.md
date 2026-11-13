# Independent reconstruction audit, 63–65D

## Results

- 64D: 335,749,696, gain 38,592.
- 65D: 335,817,922, same gain; owner lift replayed.
- 63D: 184,415,700, gain 21,712 over the published repository value.

`audit.py` imports no discovery or main-verifier logic. It enumerates binary sign kernels directly over all 65,536 local masks and verifies same-support distance >=4. For every high-overlap affected pair it enumerates both actual projected sign sets and checks disjointness; this independently replaces the discovery parity-witness proof. All other support overlaps are <=8 and hence have dot product <=32. It checks identities, no repeated supports, exact class cardinalities, and analytic interactions with dense, quadruple and pole layers. `report.json` pins the reviewed inputs.

`base_full_replay.json` is a successful fresh invocation of the entire current published base package `verify.py`, including its independent finite family review, dense-code enumeration, owner geometry replay, and quadruple designs. This is a replay of existing proof code, not an independently reimplemented complete base theorem.

`section_audit.py` independently counts explicit affine sign words and checks the rank7 projection of the 28D dense binary code. The 63D construction is a section of a valid variant of the 64D construction, not simply a section of every fixed coordinate labeling: the four-support Steiner layer is independently permuted so that its first eight labels correspond to the selected octad. This is allowed because all cross-layer bounds with this layer are support-independent. The explicit permutation is saved in `section_report.json`; direct enumeration gives 102900 retained quadruple points, whereas the unpermuted layer would give 101920. Unconstrained original sign classes are likewise individually relabeled when their overlap with the selected octad has size8; the extended Hamming projection then contributes1120 words. These untouched classes have no high-overlap interactions with changed classes, and so this independent relabeling is valid. Orthogonal restriction preserves norms and inner products and an orthonormal isometry maps to R63.

Reproduce using the four commands in the parent README. `base_full_replay.json`
records the separately executed original base verifier.

These are verified lower bounds, not claims of global optimality.
