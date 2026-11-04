"""Additive v1 mathematical-object descriptions and replay-only checked facts.

Descriptions, hashes and saved receipts never establish mathematical facts.
This module does not replace Artifact or change legacy Verification statuses.
"""
from dataclasses import dataclass
import json
from pathlib import Path

from .core import Artifact, RunContext, canonical_json, integer

SEMANTIC_TYPES = {'spherical_configuration', 'binary_code', 'ring_module',
                  'unitary_representation', 'polynomial_factorization', 'syndrome_oracle'}
ARITHMETIC_DOMAINS = {'rational', 'algebraic', 'interval', 'numeric', 'finite_field', 'finite_ring'}
ENCODINGS = {'integer_coordinates', 'rational_gram', 'implicit_recipe', 'complex_matrices',
             'ring_generator', 'polynomial_coefficients', 'syndrome_table'}
PROOF_STRENGTHS = {'unverified', 'numerical', 'exact'}
EXECUTION_OUTCOMES = {'completed', 'timeout', 'resource_limit', 'error', 'cancelled'}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _fields(value, expected, name):
    _require(isinstance(value, dict) and set(value) == set(expected), 'Invalid '+name+' fields')


def _claims(claims, ambient):
    _require(isinstance(claims, dict), 'claimed_facts must be an object')
    # Reject coercions (True, 1.0), malformed ranges and unsupported priority
    # assertions before any expensive mathematical replay.
    for key in ('dimension', 'ambient_dimension', 'rank_exact', 'rank_numerical', 'cardinality', 'length', 'min_distance'):
        if key in claims:
            integer(claims[key], key, 0 if key in ('rank_exact', 'rank_numerical', 'min_distance') else 1)
    for key in ('rank_exact', 'rank_numerical'):
        if ambient is not None and key in claims:
            _require(claims[key] <= ambient, 'Claimed rank exceeds ambient dimension')
    if ambient is not None and 'ambient_dimension' in claims:
        _require(claims['ambient_dimension'] == ambient, 'Conflicting ambient-dimension claims')
    if 'injective' in claims:
        _require(type(claims['injective']) is bool, 'Injectivity claim must be Boolean')
    for key in ('new_sota', 'global_optimality_claim', 'worldwide_priority_claim'):
        if key in claims:
            _require(claims[key] is False, 'Object metadata cannot establish '+key)
    canonical_json(claims)


class MathematicalObject:
    """Immutable JSON snapshot: semantics and representation remain distinct.

    Constructor/load APIs accept claims only. Even a serialized `checked_facts`
    field is discarded on load; call replay_object to obtain checked facts.
    """
    def __init__(self, semantics, representation, claimed_facts=None, provenance=None, symmetry=None):
        _fields(semantics, ('type', 'ambient_dimension', 'metric'), 'semantics')
        _require(isinstance(semantics['type'], str) and semantics['type'] in SEMANTIC_TYPES, 'Unknown mathematical object type')
        ambient = semantics['ambient_dimension']
        if ambient is not None:
            integer(ambient, 'ambient_dimension')
        _require(isinstance(semantics['metric'], str) and bool(semantics['metric']), 'Metric description required')
        _fields(representation, ('encoding', 'arithmetic', 'artifact'), 'representation')
        _require(isinstance(representation['encoding'], str) and representation['encoding'] in ENCODINGS, 'Unknown representation encoding')
        arithmetic = representation['arithmetic']
        _fields(arithmetic, ('domain', 'description'), 'arithmetic')
        _require(isinstance(arithmetic['domain'], str) and arithmetic['domain'] in ARITHMETIC_DOMAINS, 'Unknown arithmetic domain')
        _require(isinstance(arithmetic['description'], str) and bool(arithmetic['description']), 'Arithmetic description required')
        Artifact.from_dict(representation['artifact'])  # Structural validation, not proof.
        claims = {} if claimed_facts is None else claimed_facts
        _claims(claims, ambient)
        provenance = {'classification': 'unspecified', 'sources': [], 'description': ''} if provenance is None else provenance
        _fields(provenance, ('classification', 'sources', 'description'), 'provenance')
        _require(isinstance(provenance['classification'], str) and provenance['classification'] in {'unspecified', 'known_result', 'reconstruction', 'search_candidate', 'new_candidate'}, 'Unknown provenance classification')
        _require(isinstance(provenance['sources'], list) and all(isinstance(x, str) for x in provenance['sources']), 'Provenance sources must be strings')
        _require(isinstance(provenance['description'], str), 'Provenance description must be text')
        symmetry = {'description': '', 'claimed_group_order': None, 'generator_references': []} if symmetry is None else symmetry
        _fields(symmetry, ('description', 'claimed_group_order', 'generator_references'), 'symmetry')
        _require(isinstance(symmetry['description'], str), 'Symmetry description must be text')
        if symmetry['claimed_group_order'] is not None:
            integer(symmetry['claimed_group_order'], 'claimed_group_order')
        _require(isinstance(symmetry['generator_references'], list) and all(isinstance(x, str) for x in symmetry['generator_references']), 'Symmetry references must be strings')
        self._snapshot = canonical_json({'schema_version': 1, 'semantics': semantics,
            'representation': representation, 'claimed_facts': claims,
            'provenance': provenance, 'symmetry': symmetry})

    def to_dict(self):
        return json.loads(self._snapshot)

    @classmethod
    def from_dict(cls, raw):
        _require(isinstance(raw, dict), 'Object document must be an object')
        required = {'schema_version', 'semantics', 'representation', 'claimed_facts', 'provenance', 'symmetry'}
        _require(required <= set(raw) and set(raw) <= required | {'checked_facts', 'proof_strength'}, 'Invalid object document fields')
        _require(type(raw['schema_version']) is int and raw['schema_version'] == 1, 'Unsupported object schema')
        if 'checked_facts' in raw:
            _require(isinstance(raw['checked_facts'], dict), 'Invalid saved checked_facts')
        if 'proof_strength' in raw:
            _require(isinstance(raw['proof_strength'], str) and raw['proof_strength'] in PROOF_STRENGTHS, 'Invalid saved proof strength')
        # Deliberately never copy saved checked facts/status into a proof object.
        return cls(*(raw[key] for key in ('semantics', 'representation', 'claimed_facts', 'provenance', 'symmetry')))

    @property
    def claimed_facts(self):
        return self.to_dict()['claimed_facts']

    @property
    def checked_facts(self):
        return {}


@dataclass(frozen=True)
class ExecutionOutcome:
    status: str
    detail: str = ''

    def __post_init__(self):
        _require(isinstance(self.status, str) and self.status in EXECUTION_OUTCOMES, 'Unknown execution outcome')
        _require(isinstance(self.detail, str), 'Execution detail must be text')


class ReplayAssessment:
    """Not deserializable: created only after the existing verifier has run."""
    def __init__(self, *args, **kwargs):
        raise TypeError('Use replay_object; saved facts and hashes cannot create a proof assessment')

    @classmethod
    def _from_replay(cls, facts, strength, method):
        result = object.__new__(cls)
        result._facts_snapshot = canonical_json(facts)
        result.proof_strength = strength
        result.method = method
        result.outcome = ExecutionOutcome('completed')
        return result

    @property
    def checked_facts(self):
        return json.loads(self._facts_snapshot)

    def to_dict(self):
        return {'checked_facts': self.checked_facts, 'proof_strength': self.proof_strength,
                'method': self.method, 'execution_outcome': self.outcome.status}


@dataclass(frozen=True)
class UncheckedAttempt:
    """A timeout/resource limit says nothing about validity or exhaustion."""
    outcome: ExecutionOutcome

    def __post_init__(self):
        _require(isinstance(self.outcome, ExecutionOutcome), 'ExecutionOutcome required')

    @property
    def proof_strength(self):
        return 'unverified'

    @property
    def checked_facts(self):
        return {}


# Explicitly reviewed adapter coverage. Other semantic/encoding combinations
# can be described, but are not silently promoted by an invented adapter.
REPLAY_ADAPTERS = {
    ('spherical_configuration', 'integer_coordinates'): ('rational', 'spherical_code'),
    ('spherical_configuration', 'rational_gram'): ('rational', 'gram_spherical_code'),
    ('spherical_configuration', 'implicit_recipe'): ('algebraic', 'implicit_spherical_code'),
}


def replay_object(document, context=None):
    """Recompute checked facts from an existing independent artifact verifier.

    No numerical/interval arithmetic backend is supplied by this schema. The
    exact adapters above require an exact legacy verifier result. A numerical
    observation is not relabeled as an exact mathematical fact.
    """
    _require(isinstance(document, MathematicalObject), 'MathematicalObject required')
    raw = document.to_dict()
    semantics, representation = raw['semantics'], raw['representation']
    _require(semantics['metric'] == 'normalized_inner_product', 'No reviewed adapter for this semantic metric')
    key = (semantics['type'], representation['encoding'])
    if key not in REPLAY_ADAPTERS:
        raise NotImplementedError('No reviewed replay adapter for '+str(key))
    domain, kind = REPLAY_ADAPTERS[key]
    _require(representation['arithmetic']['domain'] == domain, 'Arithmetic domain is unsupported by this replay adapter')
    artifact = Artifact.from_dict(representation['artifact'])
    _require(artifact.kind == kind, 'Semantic type and representation artifact disagree')
    if representation['encoding'] == 'integer_coordinates':
        _require('coordinates' in artifact.data, 'Integer-coordinate adapter requires coordinates')
    if representation['encoding'] == 'implicit_recipe':
        _require(artifact.data.get('embedding') == 'binary_sign' or artifact.data.get('backend') == 'signed_support_layers', 'No reviewed implicit-recipe adapter for this backend')
    # Artifacts carry legacy claims too; the existing verifier independently
    # checks those, rather than trusting the old status or a saved digest.
    from .certification import verify_artifact
    report = verify_artifact(artifact, context or RunContext(Path('.')))
    _require(report.status == 'exact', 'This adapter requires exact independent replay')
    computed = report.details.get('facts', report.details)
    _require(isinstance(computed, dict), 'Verifier did not supply checked facts')
    ambient = semantics['ambient_dimension']
    if ambient is not None:
        _require(computed.get('ambient_dimension') == ambient, 'Declared ambient dimension disagrees with replay')
    for key, value in raw['claimed_facts'].items():
        _require(key in computed, 'No checked evidence for claimed '+key)
        _require(canonical_json(value) == canonical_json(computed[key]), 'Claimed '+key+' disagrees with replay')
    return ReplayAssessment._from_replay(computed, 'exact', report.method)
