"""Exact compressed dense-sign / signed-octad / root construction family.

For dense length m>=32 and ambient n>=m, every point has squared norm m:
* dense binary signs on the first m coordinates, zero elsewhere;
* all even signs on each octad, amplitude sqrt(m/8);
* explicit antipodal octad lines at the same amplitude;
* all four signs on every coordinate pair, amplitude sqrt(m/2).
No full point array or quadratic point-pair array is materialized. This is a
specific dimension-parametrized family, not a checker for arbitrary layers.
"""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from ..core import Artifact, RunContext, Verification, check_params, digest, integer

SOURCE32 = 'verifications/improved/dim32-35-signed-layers/source/classic/certificate32.json'
SOURCE32_SHA256 = '3a96dc79d21d07557a221f2cdd6664d2d5270e93d7c98711d7b962a2ba2b41d2'
# The transformed mathematical recipe is pinned independently of saved receipts.
RECIPE32_SHA256 = '8f90efc77509656f83a82661a151054e4ca1c065a737dd800280697c25e659bb'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _rank(rows):
    pivots = {}
    for word in rows:
        while word:
            bit = word.bit_length() - 1
            if bit in pivots:
                word ^= pivots[bit]
            else:
                pivots[bit] = word
                break
    return len(pivots)


def _certificate32(raw):
    return {'format': 'signed_support_layers_v1', 'dimension': raw['dimension'],
            'dense_length': 32, 'dense_generator': raw['dense_generator'],
            'full_even_supports': raw['full_even_supports'],
            'exceptional_antipodal_lines': raw['exceptional_antipodal_lines'],
            'expected_cardinality': raw['points']}


def _provenance32():
    return {'classification': 'known_recipe_benchmark', 'source': SOURCE32,
            'source_sha256': SOURCE32_SHA256,
            'attribution': 'Cheng-Sloane dense code; Qihong Lin support packing; repository signed replacements; known 352800-point benchmark'}


def _check_certificate(certificate, context):
    fields = {'format', 'dimension', 'dense_length', 'dense_generator',
              'full_even_supports', 'exceptional_antipodal_lines', 'expected_cardinality'}
    require(isinstance(certificate, dict) and set(certificate) == fields, 'Invalid layer certificate fields')
    require(certificate['format'] == 'signed_support_layers_v1', 'Unsupported layer certificate format')
    n = integer(certificate['dimension'], 'dimension', 32)
    m = integer(certificate['dense_length'], 'dense_length', 32)
    require(m <= n, 'Dense support exceeds ambient dimension')
    context.limit('max_dimension', n)
    generators = certificate['dense_generator']
    full = certificate['full_even_supports']
    raw_lines = certificate['exceptional_antipodal_lines']
    require(isinstance(generators, list), 'Dense generator must be a list')
    require(isinstance(full, list) and isinstance(raw_lines, list), 'Sparse support lists required')
    # Bound counts before any exponential code or pair work. At most m
    # independent binary generators exist; reject dependencies separately.
    require(len(generators) <= m, 'Too many dense generators')
    context.limit('max_codewords', 1 << len(generators))
    context.limit('max_materialized_points', len(full) + len(raw_lines))
    for word in generators:
        integer(word, 'dense generator mask')
        require(word < 1 << m, 'Dense generator outside dense coordinates')
    require(_rank(generators) == len(generators), 'Dependent dense generators')
    for mask in full:
        integer(mask, 'full support mask')
        require(mask < 1 << n and mask.bit_count() == 8, 'Full support must be an ambient octad')
    require(len(set(full)) == len(full), 'Repeated full support')
    lines = []
    for row in raw_lines:
        require(isinstance(row, list) and len(row) == 2, 'Exceptional line needs support and sign masks')
        mask, signs = row
        integer(mask, 'exceptional support mask')
        integer(signs, 'exceptional sign mask', 0)
        require(mask < 1 << n and mask.bit_count() == 8, 'Exceptional support must be an ambient octad')
        require(signs & mask == signs, 'Exceptional signs outside support')
        require(signs & (mask & -mask) == 0, 'Antipodal line representative must be canonical')
        lines.append((mask, signs))
    require(len(set(lines)) == len(lines), 'Repeated exceptional antipodal line')
    exceptional_supports = {mask for mask, _ in lines}
    require(not exceptional_supports.intersection(full), 'Exceptional and complete support layers overlap')
    compressed_pairs = (len(full) * (len(full)-1)//2 + len(lines) * (len(lines)-1)//2
                        + len(full) * len(exceptional_supports))
    context.limit('max_pairs', compressed_pairs)

    # Enumerate only the binary linear code, never its pair graph. Differences
    # of codewords are codewords, so nonzero weights certify all dense pairs.
    words = {0}
    for generator in generators:
        words |= {word ^ generator for word in words}
    require(len(words) == 1 << len(generators), 'Dense span cardinality mismatch')
    distance = min((word.bit_count() for word in words if word), default=None)
    require(distance is None or 4 * distance >= m, 'Dense code violates kissing distance')
    full_max = max(((a & b).bit_count() for a, b in combinations(full, 2)), default=0)
    require(full_max <= 4, 'Full-support intersection conflict')
    cross_max = max(((a & b).bit_count() for a in exceptional_supports for b in full), default=0)
    require(cross_max <= 4, 'Exceptional/full intersection conflict')
    signed_max = 0
    for (a, signs_a), (b, signs_b) in combinations(lines, 2):
        overlap = a & b
        dot = overlap.bit_count() - 2 * (overlap & (signs_a ^ signs_b)).bit_count()
        signed_max = max(signed_max, abs(dot))
        require(abs(dot) <= 4, 'Exceptional signed conflict')

    # Counts follow from independent generators, even parity (2^7 octad
    # signs), canonical antipodal lines, all coordinate-pair roots, and disjoint
    # support sizes m>=32, 8, 2. Root families span all n coordinates exactly.
    components = {'dense': len(words), 'full_even_octads': 128 * len(full),
                  'exceptional_octad_points': 2 * len(lines), 'roots': 2 * n * (n-1)}
    cardinality = sum(components.values())
    integer(certificate['expected_cardinality'], 'expected_cardinality')
    require(cardinality == certificate['expected_cardinality'], 'Layer cardinality claim mismatch')
    facts = {'dimension': n, 'ambient_dimension': n, 'cardinality': cardinality,
             'norm_squared': m, 'max_inner_product_upper_bound': str(Fraction(m, 2)),
             'dimension_scope': 'exact; complete two-coordinate root families span the ambient space'}
    proof = {
        'components': components, 'dense_code_rank': len(generators),
        'dense_codewords_enumerated': len(words), 'dense_min_distance': distance,
        'compressed_support_pairs_checked': compressed_pairs,
        'full_intersection_max': full_max, 'exceptional_full_intersection_max': cross_max,
        'exceptional_signed_dot_absolute_max': signed_max,
        'pair_bounds': {
            'dense_dense': 'm-2*d <= m/2 from 4*d>=m',
            'full_same_support': 'even sign parity gives at least two flips: m-4*(m/8)=m/2',
            'full_different_supports': '(m/8)*intersection <= m/2',
            'exceptional_full': '(m/8)*intersection <= m/2',
            'exceptional_exceptional': '(m/8)*abs(signed_overlap) <= m/2; antipodes have product -m',
            'roots_roots': 'distinct non-antipodal roots share at most one coordinate: product <= m/2',
            'dense_middle': 'absolute product <= sqrt(8*m) <= m/2 because m>=32',
            'dense_roots': 'absolute product <= 2*sqrt(m/2) <= m/2 because m>=8',
            'middle_roots': 'absolute product <= 2*sqrt(m/8)*sqrt(m/2)=m/2'},
        'distinctness': 'independent binary span; distinct support masks; canonical distinct antipodal lines; disjoint support sizes',
        'arithmetic': 'Python integers; symbolic positive square-root bounds use exact squared inequalities',
        'materialized_spherical_points': 0,
        'scaling': 'O(2^k) integer codewords and O(F^2+E^2+F*S) sparse pair checks; no O(M^2) spherical Gram matrix'}
    return facts, proof


def verify_artifact(artifact, context=None):
    context = context or RunContext(Path('.'))
    require(artifact.kind == 'implicit_spherical_code', 'Layer verification requires an implicit spherical code')
    data = artifact.data
    require(isinstance(data, dict) and set(data) == {'backend', 'certificate', 'provenance'}, 'Invalid layer payload')
    require(data['backend'] == 'signed_support_layers', 'Unknown layer backend')
    provenance = data['provenance']
    require(isinstance(provenance, dict), 'Invalid layer provenance')
    if provenance.get('classification') == 'known_recipe_benchmark':
        require(provenance == _provenance32(), 'Modified known-recipe provenance')
        require(digest(data['certificate']) == RECIPE32_SHA256, 'Modified known32 mathematical payload')
    elif provenance.get('classification') == 'signed_fiber_search':
        require(set(provenance) == {'classification', 'solver_status', 'candidate_lines', 'flexible_fibers'}, 'Invalid search provenance')
        require(isinstance(provenance['solver_status'], str) and len(provenance['solver_status']) <= 500, 'Invalid solver status')
        require(type(provenance['candidate_lines']) is int and provenance['candidate_lines'] >= 0, 'Invalid candidate count')
        require(type(provenance['flexible_fibers']) is int and provenance['flexible_fibers'] >= 0, 'Invalid flexible fiber count')
    else:
        require(provenance == {'classification': 'provided_recipe'}, 'Unsupported provenance claim')
    facts, proof = _check_certificate(data['certificate'], context)
    for key in ('dimension', 'ambient_dimension', 'cardinality', 'norm_squared', 'max_inner_product_upper_bound'):
        if key in artifact.facts:
            require(type(artifact.facts[key]) is type(facts[key]) and artifact.facts[key] == facts[key],
                    'Claimed layer ' + key + ' disagrees with replay')
    return Verification('exact', 'native compressed signed-support layer inequalities and exact binary-span distance',
                        {'facts': facts, 'proof': proof, 'data_sha256': digest(data),
                         'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                         'source': deepcopy(provenance), 'materialized': False})


def signed_supports(inputs, params, context):
    require(not inputs, 'layers.signed_supports takes no input artifacts')
    check_params(params, ('certificate',))
    require('certificate' in params, 'A signed-support certificate is required')
    artifact = Artifact('implicit_spherical_code', {'backend': 'signed_support_layers',
                        'certificate': deepcopy(params['certificate']),
                        'provenance': {'classification': 'provided_recipe'}})
    artifact.verification = verify_artifact(artifact, context)
    artifact.facts = artifact.verification.details['facts']
    return artifact


def known32(inputs, params, context):
    require(not inputs, 'layers.known32 takes no inputs')
    check_params(params, ())
    source = Path(__file__).resolve().parents[2] / SOURCE32
    raw = source.read_bytes()
    require(sha256(raw).hexdigest() == SOURCE32_SHA256, 'Known32 source certificate hash mismatch')
    artifact = Artifact('implicit_spherical_code', {'backend': 'signed_support_layers',
                        'certificate': _certificate32(json.loads(raw)), 'provenance': _provenance32()})
    artifact.verification = verify_artifact(artifact, context)
    artifact.facts = artifact.verification.details['facts']
    return artifact


OPERATORS = {'layers.known32': known32, 'layers.signed_supports': signed_supports}
INPUT_KINDS = {'layers.known32': (), 'layers.signed_supports': ()}

