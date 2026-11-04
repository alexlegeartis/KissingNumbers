"""Additive, explicitly partial operator contracts; no execution or auto-proof.

Contract metadata documents reviewed semantics. The existing allow-listed DAG
and independent checkers remain the execution/certification authorities.
"""
from dataclasses import dataclass
from types import MappingProxyType
import re

from .core import canonical_json
from .objects import SEMANTIC_TYPES

SCOPES = {'finite', 'enumerable', 'heuristic'}
PARAMETER_DOMAINS = {'integer_matrix', 'positive_integer', 'positive_integer_list',
                     'rational_square_matrix', 'positive_rational'}
ESTIMATORS = {'integer_coordinate_shapes', 'gram_shapes', 'recursive_binary_proof'}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class OperatorContract:
    name: str
    input_types: tuple
    output_type: str
    parameter_domains: dict
    scope: str
    resource_estimator: str
    proof_obligations: tuple
    limitations: tuple
    schema_version: int = 1

    def __post_init__(self):
        _require(type(self.schema_version) is int and self.schema_version == 1, 'Unsupported contract schema')
        _require(isinstance(self.name, str) and re.fullmatch(r'[a-z][a-z0-9_]*\.[a-z][a-z0-9_]*', self.name) is not None, 'Invalid contract operator name')
        _require(type(self.input_types) is tuple and all(type(group) is tuple and bool(group) and all(isinstance(kind, str) and kind in SEMANTIC_TYPES for kind in group) and len(set(group)) == len(group) for group in self.input_types), 'Each input needs a tuple of accepted semantic types')
        _require(isinstance(self.output_type, str) and self.output_type in SEMANTIC_TYPES, 'Unknown semantic output type')
        _require(isinstance(self.scope, str) and self.scope in SCOPES, 'Unknown operation scope')
        _require(isinstance(self.resource_estimator, str) and self.resource_estimator in ESTIMATORS, 'No reviewed resource estimator')
        _require(isinstance(self.parameter_domains, dict), 'Parameter domains must be an object')
        for name, domain in self.parameter_domains.items():
            _require(isinstance(name, str) and isinstance(domain, dict) and set(domain) == {'domain', 'required', 'description'}, 'Malformed parameter domain')
            _require(isinstance(domain['domain'], str) and domain['domain'] in PARAMETER_DOMAINS and type(domain['required']) is bool and isinstance(domain['description'], str) and bool(domain['description']), 'Unsupported parameter domain')
        for values in (self.proof_obligations, self.limitations):
            _require(type(values) is tuple and bool(values) and all(isinstance(value, str) and value for value in values), 'Proof obligations and limits must be nonempty text tuples')
        # Snapshot mutable caller-owned domain metadata. Returned metadata is
        # always a copy; contracts are not an execution dispatch registry.
        object.__setattr__(self, '_domains_json', canonical_json(self.parameter_domains))
        object.__setattr__(self, 'parameter_domains', MappingProxyType({name: MappingProxyType(dict(domain)) for name, domain in self.parameter_domains.items()}))

    def to_dict(self):
        import json
        return {'schema_version': self.schema_version, 'name': self.name,
                'input_types': [list(group) for group in self.input_types],
                'output_type': self.output_type, 'parameter_domains': json.loads(self._domains_json),
                'scope': self.scope, 'resource_estimator': self.resource_estimator,
                'proof_obligations': list(self.proof_obligations), 'limitations': list(self.limitations)}


def _domain(kind, required, description):
    return {'domain': kind, 'required': required, 'description': description}


CONTRACTS = MappingProxyType({contract.name: contract for contract in (
    OperatorContract('geometry.integer_code', (), 'spherical_configuration', {
        'coordinates': _domain('integer_matrix', True, 'Nonempty rectangular integer rows; booleans and floats rejected'),
        'weights': _domain('positive_integer_list', False, 'One strictly positive diagonal quadratic-form weight per column'),
        'norm_squared': _domain('positive_integer', True, 'Common squared norm in the stated diagonal form')},
        'finite', 'integer_coordinate_shapes', (
            'Distinct nonzero rows and positive metric', 'Every norm equals the claimed common norm',
            'Every distinct pair product is at most half that norm',
            'Exact rank is independently recomputed; column count is only ambient dimension'), (
            'Explicit finite configuration only; pair and exact-rank budgets apply',
            'No group finiteness, injectivity of an external map, optimality or novelty conclusion')),
    OperatorContract('geometry.gram_code', (), 'spherical_configuration', {
        'gram': _domain('rational_square_matrix', True, 'Square exact rational matrix; entries integers or numerator/denominator pairs'),
        'norm_squared': _domain('positive_rational', False, 'Optional exact positive rational norm claim equal to the common Gram diagonal')},
        'finite', 'gram_shapes', (
            'Exact symmetry and common positive diagonal', 'Exact positive semidefiniteness and rank',
            'Distinctness and every off-diagonal kissing inequality',
            'Existence of a real realization follows from PSD; rational Gram entries do not imply rational Cartesian coordinates'), (
            'Symbolic real realization, with no coordinate factorization supplied',
            'Rational-size, Gram-entry and rank-work budgets apply; no interval arithmetic backend')),
    OperatorContract('geometry.sign_embed', (('binary_code',),), 'spherical_configuration', {},
        'finite', 'recursive_binary_proof', (
            'Independently replay the binary child representation and its cardinality/distance',
            'Positive binary length and distinct words', 'Distance d and length N satisfy 4*d>=N',
            'Sign embedding is injective on binary words; rank remains unknown unless separately proved'), (
            'Input backend must already have a supported independent verifier',
            'Ambient length does not establish full span; recursive proof cost can dominate')),
)})


def contract_coverage(operator_names=None):
    """Report missing metadata rather than inventing contracts for new operators."""
    if operator_names is None:
        from .runner import registry
        operator_names = registry()[0]
    _require(isinstance(operator_names, (dict, list, tuple, set, frozenset)), 'Operator-name collection required')
    names = set(operator_names)
    _require(all(isinstance(name, str) for name in names), 'Operator names must be strings')
    return {'schema_version': 1, 'complete': names <= set(CONTRACTS),
            'covered': sorted(names & set(CONTRACTS)), 'uncovered': sorted(names - set(CONTRACTS)),
            'scope': 'Only explicitly reviewed contracts; coverage is not certification or global construction coverage'}


def estimate_resources(name, params, input_objects=()):
    """Conservative shape information, not a complete execution upper bound.

    All counts come from representation shapes, never saved mathematical claims.
    Exact-rank bit growth and recursive backend costs remain explicitly unknown.
    """
    _require(isinstance(name, str) and name in CONTRACTS, 'No reviewed operator contract')
    _require(isinstance(params, dict), 'Parameters must be an object')
    contract = CONTRACTS[name]
    _require(len(input_objects) == len(contract.input_types), 'Wrong input count for estimate')
    if contract.resource_estimator == 'recursive_binary_proof':
        return {'materialized_spherical_points': 0, 'recursive_verifier_work': None,
                'complete_upper_bound': False, 'reason': 'Binary child must be independently replayed'}
    rows = params.get('coordinates' if contract.resource_estimator == 'integer_coordinate_shapes' else 'gram')
    _require(isinstance(rows, list) and bool(rows) and isinstance(rows[0], list) and bool(rows[0]), 'Nonempty matrix shape required')
    n, d = len(rows), len(rows[0])
    _require(all(isinstance(row, list) and len(row) == d for row in rows), 'Ragged matrix')
    if contract.resource_estimator == 'gram_shapes':
        _require(n == d, 'Gram matrix must be square')
    return {'materialized_spherical_points': n if contract.resource_estimator == 'integer_coordinate_shapes' else 0,
            'matrix_entries': n*d, 'pair_inequalities': n*(n-1)//2,
            'maximum_pair_block_entries': min(128, n)*n if contract.resource_estimator == 'integer_coordinate_shapes' else n*n,
            'rank_bit_work': None, 'complete_upper_bound': False,
            'reason': 'Shape counts only; exact elimination bit growth is not predicted'}
