# Mathematical objects and operator contracts: first additive increment

This increment adds executable schemas and replay adapters in `discovery/objects.py`, operator metadata in `discovery/contracts.py`, and adversarial tests in `discovery/tests/test_object_contracts.py`. It preserves the existing `Artifact`, allow-listed construction DAG, verifier dispatch, receipt format and three legacy verification statuses. It does not introduce a replacement runner or enable an unimplemented arithmetic backend.

## Object semantics and representation

A version-1 `MathematicalObject` separates these fields:

| Field | Meaning and trust |
|---|---|
| `semantics` | Object type, stated ambient dimension and metric. A dimension is a claim until replay compares it with the representation. No family/dimension dispatch table is imposed. |
| `representation` | Encoding, explicit arithmetic domain/description, and a structurally checked legacy Artifact. Integer coefficient coordinates with a diagonal metric do not promise rational Cartesian unit coordinates. A rational Gram does not promise a rational coordinate factorization. |
| `claimed_facts` | Claims to compare with independently recomputed facts. An unsupported claim such as injectivity is rejected by the current replay adapter rather than accepted as metadata. |
| `provenance` | Claimed attribution and sources; reconstruction, known result and candidate remain distinct descriptions. Provenance cannot establish novelty or worldwide priority. |
| `symmetry` | Described generators and claimed group order. These remain claims; the adapter does not infer finiteness, orbit completeness or an injective action from them. |

Exact span rank belongs in the checked fact `rank_exact`; it is separate from `ambient_dimension`. The coordinate adapter demonstrates two antipodal vectors in 96 coordinate dimensions with exact rank one. A Gram adapter obtains a symbolic real realization in the exact PSD rank; it does not conflate matrix order with geometric dimension.

Descriptors explicitly support rational, algebraic, interval and numeric arithmetic, plus finite-field/ring representation domains. **Interval and numeric descriptions do not implement arithmetic or certify enclosures.** Current reviewed adapters cover only spherical integer-coordinate artifacts, rational Gram artifacts, and the existing exact binary-sign/signed-layer implicit backends. Other semantic/encoding combinations raise `NotImplementedError`; mismatched arithmetic domains are rejected. A numeric input cannot become exact because a descriptor or saved receipt says so.

`MathematicalObject` takes immutable JSON snapshots and returns copies. Its `checked_facts` are always empty. `replay_object(document, context)` calls the existing independent verifier, requires its exact result, compares declared ambient dimension and every new claimed fact, and only then returns a `ReplayAssessment` with checked facts. On deserialization, saved `checked_facts` and saved `proof_strength` are deliberately discarded. A hash identifies data; it is never substituted for replay. These protections address untrusted serialized inputs, not arbitrary hostile Python code running inside the process.

For example:

```python
from discovery.core import Artifact
from discovery.objects import MathematicalObject, replay_object

artifact = Artifact('spherical_code', {
    'coordinates': [[1, 0], [-1, 0]], 'norm_squared': 1,
})
obj = MathematicalObject(
    {'type': 'spherical_configuration', 'ambient_dimension': 2,
     'metric': 'normalized_inner_product'},
    {'encoding': 'integer_coordinates',
     'arithmetic': {'domain': 'rational', 'description': 'integer coefficients and exact diagonal metric'},
     'artifact': artifact.to_dict()},
    {'cardinality': 2, 'rank_exact': 1},
)
checked = replay_object(obj)
assert checked.checked_facts['ambient_dimension'] == 2
assert checked.checked_facts['rank_exact'] == 1
```

The existing verifier still enforces resource budgets and the kissing inequality. The new schema does not widen that verifier into a checker for arbitrary spherical angles or arbitrary implicit descriptions.

## Operator contracts and outcomes

Version-1 `OperatorContract` metadata gives accepted semantic input types, a semantic output type, parameter-domain descriptions, scope (`finite`, `enumerable`, or `heuristic`), a reviewed resource-estimator identifier, proof obligations and limitations. Each input position can list alternative semantic types. These are not inferred from an operator's name, advertised result or legacy kind strings.

The first reviewed contracts are `geometry.integer_code`, `geometry.gram_code` and `geometry.sign_embed`. `contract_coverage()` compares this explicit set with the current runner registry and lists all uncovered operators, including operators added concurrently. Partial coverage is intentional: this increment does not invent contracts for subgroup searches, stabilizers, ring lifts, sections, root systems or general gluing before their exact semantics and bounds are reviewed. Contract metadata is separate from the execution allow list and does not authorize an operation or guarantee success.

Resource estimates use representation shapes, rather than saved cardinality claims. They expose exact matrix/pair shape counts where available, mark exact-rank bit growth and recursive verifier costs unknown, and return `complete_upper_bound=False`. They are planning information; the existing `RunContext` and independent verifiers enforce actual limits. A finite construction contract does not imply enumeration of every construction, and a heuristic contract would not imply exhaustive search.

`ExecutionOutcome` reports completion, timeout, resource limit, error or cancellation independently of proof strength. `UncheckedAttempt` stays unverified with empty checked facts for all of these outcomes. A completed run alone is not a proof. Timeout/resource exhaustion establishes neither invalidity nor exhaustive failure. The old `Verification` enum is unchanged; no interval status or new arithmetic implementation is advertised.

## Critical assessment and next increment

The current DAG is a useful collection of typed constructions and independent checkers, but its legacy Artifact kind alone does not carry enough mathematical semantics for universal composition. Two `spherical_code` objects may use different coordinate bases, diagonal metrics, norms, arithmetic domains and phase conventions. An implicit code can prove a cardinality without presenting an enumerable point set. A generator matrix does not establish independence, and a described group does not establish finiteness. Symbolic cardinality transport additionally requires an injective map or verified fiber counts. Hash equality, stored metadata and numerical near equality cannot discharge those obligations.

Supporting objects in dimensions 1–96 does not establish coverage of every construction in those dimensions. This stage adds no family dimension filters; unavailable representation/checker combinations remain explicitly unavailable. It also makes no claim to have implemented generic interval certification, algebraic-number-field arithmetic, general group enumeration, or arbitrary code/lattice conversion.

The next contract increment should formalize heterogeneous gluing with verified overlap accounting. Its obligations should include:

1. Verify each child object with its own backend and establish a common real realization, norm and metric. Separate Gram blocks need verified cross blocks and a PSD whole, or an exact isometry/embedding certificate; choosing an arbitrary relative coordinate alignment is not a proof.
2. Verify every cross-family angle through a checked pair oracle or a bounded explicit Gram block. A transitive-orbit reduction needs a verified action and actual vector phases, not merely projective line labels.
3. Prove each embedding injective, or provide exact fiber counts. Detect normalized-vector equality separately from line equality and ordinary angle compatibility. Duplicate unit points have inner product one and cannot be counted twice.
4. Compute union cardinality from disjoint ownership or a verified intersection certificate. For two families, `|A union B| = |A|+|B|-|A intersection B|`; for more families, pairwise overlaps alone are insufficient because triple intersections matter. Prefer canonical ownership or a full checked intersection decomposition.
5. Derive ambient dimension and exact rank only when supported by the realized common object. Carry source attribution and the logical scope of each child proof without turning known-result reconstruction into discovery.
6. Predict the cost of child replay and all required cross/overlap certificates before allocation. If those proofs exceed a budget, report an unchecked resource-limited attempt rather than a mathematical rejection or exhausted search.

A currently implemented explicit gluing/union backend can later receive such a reviewed contract; this schema does not claim that a fully general compressed heterogeneous union operator already exists.

## Validation

Run:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discovery.tests.test_object_contracts -v
```

The 11 tests cover malformed numeric claims, rank versus ambient dimension through 96, forged saved exact receipts, noninjective metadata and duplicate coordinates, unverified symmetry, unavailable interval/numeric adapters, rational Gram realization, snapshot aliasing, timeout/resource-limit outcomes, explicit partial coverage and incomplete cost estimates. Existing DAG and Artifact tests remain applicable without modification.
