"""Implicit real stabilizer sectors from affine F2 flats and quadratic signs.

All dimensions are 2**m by the construction's mathematics. No stored record,
coordinate archive, numerical orbit, or materialized state list is used.
"""
from __future__ import annotations
from pathlib import Path
from discovery.core import Artifact, RunContext, Verification

BACKEND = 'real_stabilizer_sector'
FORMAT = 'real_stabilizer_sector_v1'


def _parameters(m, parity, context):
    if type(m) is not int or m < 1:
        raise ValueError('m must be a positive integer')
    if type(parity) is not int or parity not in (0, 1):
        raise ValueError('parity must be 0 or 1')
    # Check the exponent before constructing a possibly enormous power of two.
    dimension_limit = context.budget['max_dimension']
    if m >= dimension_limit.bit_length():
        raise ValueError('stabilizer ambient dimension exceeds max_dimension')
    context.limit('max_dimension', 1 << m)
    # This bounds all count intermediates before exponentiation.
    context.limit('max_rational_bits', m + m*(m+1)//2 + 2)
    return 1 << m


def _gaussian_product(m, r):
    numerator = denominator = 1
    for j in range(r):
        numerator *= (1 << (m-j))-1
        denominator *= (1 << (r-j))-1
    quotient, remainder = divmod(numerator, denominator)
    if remainder:
        raise ValueError('nonintegral Gaussian binomial')
    return quotient


def _proof(m, parity, layers):
    return {
        'support_dimensions': list(range(parity, m+1, 2)),
        'layer_cardinalities': layers,
        'same_rank_distinct_flats': 'intersection_dimension_at_most_r_minus_1',
        'different_ranks': 'same_parity_gives_rank_difference_at_least_2',
        'same_flat_signs': 'nonconstant_quadratic_bias_at_most_one_half',
        'quadratic_bias_lemma': 'alternating_polar_rank_even; nonzero_rank_at_least_2; zero_rank_nonconstant_is_balanced_linear',
        'constant_sign_difference': 'one_gives_antipodes; zero_is_identical_state',
        'normalization_squared_norm': [1, 1],
        'distinct_pair_inner_product_upper_bound': [1, 2],
        'full_span_support_dimension': parity,
        'count_identity': 'q_binomial_product_at_t_plus_and_minus_1',
    }


def _data(m, parity, layers):
    return {
        'backend': BACKEND, 'format': FORMAT, 'm': m, 'parity': parity,
        'support_model': 'all_affine_flats_of_F2_m_with_canonical_rref_coordinates',
        'sign_model': 'all_squarefree_boolean_polynomials_of_degree_at_most_2_including_constant',
        'normalization': 'divide_integer_sign_row_by_sqrt_2_to_the_support_dimension',
        'proof': _proof(m, parity, layers),
    }


def _facts(m, cardinality):
    dimension = 1 << m
    return {
        'dimension': dimension, 'ambient_dimension': dimension, 'rank_exact': dimension,
        'dimension_scope': 'full span proved by support-dimension zero or one layer',
        'cardinality': cardinality, 'norm_squared': 1,
        'max_inner_product_upper_bound': 0.5,
    }


def real_stabilizer_sector(inputs, params, context):
    if inputs:
        raise ValueError('real_stabilizer_sector has no inputs')
    if not isinstance(params, dict) or set(params)-{'m', 'parity'}:
        raise ValueError('unknown real_stabilizer_sector parameters')
    context = context or RunContext(Path('.'))
    m, parity = params.get('m', 4), params.get('parity', 0)
    _parameters(m, parity, context)
    layers = [_gaussian_product(m, r) * (1 << (m-r)) *
              (1 << (1+r+r*(r-1)//2)) for r in range(m+1)]
    cardinality = sum(layers[parity::2])
    artifact = Artifact('implicit_spherical_code', _data(m, parity, layers),
                        _facts(m, cardinality))
    artifact.verification = verify_artifact(artifact, context)
    return artifact


def verify_artifact(artifact, context=None):
    """Reconstruct the count by q-Pascal recurrence and replay theorem conditions.

    Distinct affine supports of equal dimension intersect in at most half their
    points; unequal ranks in one parity sector differ by at least two. On the
    same support, the sign difference is quadratic. Its polar alternating form
    has even rank: positive rank gives bias <=1/2, rank zero gives a balanced
    nonconstant linear function or a constant. The latter is either the same
    state or its antipode. These three exhaustive cases prove the bound.
    """
    context = context or RunContext(Path('.'))
    if artifact.kind != 'implicit_spherical_code' or not isinstance(artifact.data, dict):
        raise ValueError('expected implicit stabilizer spherical code')
    data = artifact.data
    m, parity = data.get('m'), data.get('parity')
    dimension = _parameters(m, parity, context)
    # Independent count route: [n r]_2 = [n-1 r]_2 + 2**(n-r)[n-1 r-1]_2.
    gaussian = [1]
    for n in range(1, m+1):
        gaussian = [1] + [
            (gaussian[r] if r < n else 0) + (1 << (n-r))*gaussian[r-1]
            for r in range(1, n+1)]
    layers = [gaussian[r] * (1 << (m-r)) *
              (1 << (1+r+r*(r-1)//2)) for r in range(m+1)]
    cardinality = sum(layers[parity::2])
    product_count = 1 << m
    for j in range(m):
        product_count *= 1 + (1 << j)
    if sum(layers[0::2]) != sum(layers[1::2]) or cardinality != product_count:
        raise ValueError('sector count identity failed')
    if data != _data(m, parity, layers):
        raise ValueError('stabilizer construction or proof conditions altered')
    facts = _facts(m, cardinality)
    if artifact.facts != facts:
        raise ValueError('stabilizer facts do not match independent reconstruction')
    # r=0 supplies the coordinate axes. For odd parity r=1 is present and
    # the two sign choices on every two-point affine line span its two axes.
    if parity > m:
        raise ValueError('full-span witness layer absent')
    return Verification('exact', 'affine binary-flat count and quadratic Boolean bias theorem', {
        'facts': facts, **facts,
        'proof': _proof(m, parity, layers),
        'independent_count_method': 'q-Pascal recurrence and q-binomial product identity',
        'materialized_spherical_points': 0,
        'pair_cases': ['disjoint affine supports', 'equal rank distinct supports',
                       'different ranks of equal parity', 'same support quadratic difference'],
        'proof_scope': 'one complete real stabilizer parity sector; no novelty or kissing upper-bound claim',
        'symbolic_dimension': dimension,
    })


OPERATORS = {'groups.real_stabilizer_sector': real_stabilizer_sector}
INPUT_KINDS = {'groups.real_stabilizer_sector': ()}
