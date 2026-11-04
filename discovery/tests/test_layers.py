"""Large compressed known control, direct pair oracle and rejection controls."""
from copy import deepcopy
from itertools import combinations, product
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from discovery.core import Artifact, RunContext, Verification
from discovery.engines import layers


def mask(indices):
    return sum(1 << i for i in indices)


def small_certificate(n=32, m=32):
    # Two dense signs, one complete even octad and one disjoint exceptional
    # antipodal line; the root family spans all ambient coordinates.
    return {'format': 'signed_support_layers_v1', 'dimension': n, 'dense_length': m,
            'dense_generator': [(1 << m)-1],
            'full_even_supports': [mask(range(8))],
            'exceptional_antipodal_lines': [[mask(range(8, 16)), 0]],
            'expected_cardinality': 2 + 128 + 2 + 2*n*(n-1)}


def direct_integer_points(certificate):
    """Independent coordinate oracle for m=32 only, used on the small fixture."""
    n = certificate['dimension']
    rows = [[1]*32+[0]*(n-32), [-1]*32+[0]*(n-32)]
    for support in certificate['full_even_supports']:
        indices = [i for i in range(n) if support >> i & 1]
        for signs in product((-1, 1), repeat=8):
            if sum(s < 0 for s in signs) % 2:
                continue
            row = [0]*n
            for index, sign in zip(indices, signs):
                row[index] = 2*sign
            rows.append(row)
    for support, signs in certificate['exceptional_antipodal_lines']:
        row = [(-2 if signs >> i & 1 else 2) if support >> i & 1 else 0 for i in range(n)]
        rows += [row, [-x for x in row]]
    for i, j in combinations(range(n), 2):
        for si, sj in product((-1, 1), repeat=2):
            row = [0]*n
            row[i], row[j] = 4*si, 4*sj
            rows.append(row)
    return np.asarray(rows, dtype=np.int64)


class LayerTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path('.'))

    def generic(self, certificate=None, context=None):
        return layers.signed_supports([], {'certificate': certificate or small_certificate()},
                                      context or self.context)

    def test_actual_352800_control_is_compressed_and_exact(self):
        artifact = layers.known32([], {}, self.context)
        self.assertEqual(artifact.kind, 'implicit_spherical_code')
        self.assertEqual(artifact.data['backend'], 'signed_support_layers')
        self.assertEqual((artifact.facts['dimension'], artifact.facts['cardinality']), (32, 352800))
        self.assertEqual(artifact.verification.status, 'exact')
        proof = artifact.verification.details['proof']
        self.assertEqual(proof['components'], dict(dense=131072, full_even_octads=218112,
                                                  exceptional_octad_points=1632, roots=1984))
        self.assertEqual(proof['dense_min_distance'], 8)
        self.assertEqual(proof['dense_codewords_enumerated'], 131072)
        self.assertEqual(proof['materialized_spherical_points'], 0)
        self.assertEqual(proof['compressed_support_pairs_checked'], 1826076)
        self.assertNotIn('coordinates', artifact.data)
        self.assertEqual(artifact.verification.details['source']['source_sha256'], layers.SOURCE32_SHA256)
        self.assertEqual(len(artifact.verification.details['checker_sha256']), 64)

    def test_independent_direct_all_pairs_for_smaller_layered_fixture(self):
        certificate = small_certificate()
        artifact = self.generic(certificate)
        points = direct_integer_points(certificate)
        self.assertEqual(len(points), artifact.facts['cardinality'])
        self.assertEqual(len(set(map(tuple, points))), len(points))
        self.assertTrue(np.all(np.sum(points*points, axis=1) == 32))
        maximum = -32
        for start in range(0, len(points), 64):
            products = points[start:start+64] @ points.T
            for offset, row in enumerate(products):
                off = row[start+offset+1:]
                if len(off):
                    maximum = max(maximum, int(off.max()))
        self.assertEqual(maximum, 16)

    def test_dimension_parameters_extend_to_96_and_nonintegral_scales(self):
        for n, m in ((32, 32), (48, 40), (72, 72), (96, 32), (96, 96), (96, 33)):
            with self.subTest(n=n, m=m):
                artifact = self.generic(small_certificate(n, m))
                self.assertEqual(artifact.facts['dimension'], n)
                self.assertEqual(artifact.facts['norm_squared'], m)
                self.assertEqual(artifact.facts['cardinality'], 132+2*n*(n-1))
                self.assertEqual(layers.verify_artifact(artifact).status, 'exact')

    def test_saved_verification_receipt_is_not_authority(self):
        artifact = self.generic()
        artifact = Artifact.from_dict(artifact.to_dict())
        artifact.verification = Verification('numerical', 'irrelevant receipt')
        self.assertEqual(layers.verify_artifact(artifact).status, 'exact')
        artifact.facts['cardinality'] += 1
        artifact.verification = Verification('exact', 'forged success')
        with self.assertRaises(ValueError):
            layers.verify_artifact(artifact)

    def test_any_modified_known32_payload_is_rejected(self):
        artifact = layers.known32([], {}, self.context)
        for mutation in ('count', 'generator', 'support', 'provenance'):
            bad = deepcopy(artifact)
            if mutation == 'count':
                bad.data['certificate']['expected_cardinality'] += 1
            elif mutation == 'generator':
                bad.data['certificate']['dense_generator'].reverse()
            elif mutation == 'support':
                bad.data['certificate']['full_even_supports'].reverse()
            else:
                bad.data['provenance']['source_sha256'] = '0'*64
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                layers.verify_artifact(bad)

    def test_generic_dense_distance_dependencies_and_count_are_replayed(self):
        for mutation in ('distance', 'dependent', 'count', 'outside_dense'):
            c = small_certificate()
            if mutation == 'distance':
                c['dense_generator'] = [1]
            elif mutation == 'dependent':
                c['dense_generator'] *= 2
            elif mutation == 'outside_dense':
                c['dense_generator'] = [1 << 32]
            else:
                c['expected_cardinality'] += 1
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.generic(c)

    def test_support_intersection_and_signed_conflicts(self):
        c = small_certificate()
        c['full_even_supports'].append(mask(range(1, 9)))
        c['expected_cardinality'] += 128
        with self.assertRaisesRegex(ValueError, 'intersection'):
            self.generic(c)
        c = small_certificate()
        c['exceptional_antipodal_lines'] = [[mask(range(1, 9)), 0]]
        with self.assertRaisesRegex(ValueError, 'intersection'):
            self.generic(c)
        c = small_certificate()
        c['full_even_supports'] = []
        c['exceptional_antipodal_lines'] = [[255, 0], [255, 2]]
        c['expected_cardinality'] = 2+4+1984
        with self.assertRaisesRegex(ValueError, 'signed conflict'):
            self.generic(c)

    def test_sparse_duplicate_sign_and_field_rejections(self):
        for mutation in ('duplicate_full', 'duplicate_line', 'noncanonical', 'outside_sign',
                         'wrong_support_size', 'overlap_layer', 'unknown', 'boolean'):
            c = small_certificate()
            if mutation == 'duplicate_full': c['full_even_supports'] *= 2
            elif mutation == 'duplicate_line': c['exceptional_antipodal_lines'] *= 2
            elif mutation == 'noncanonical': c['exceptional_antipodal_lines'][0][1] = 1 << 8
            elif mutation == 'outside_sign': c['exceptional_antipodal_lines'][0][1] = 1 << 3
            elif mutation == 'wrong_support_size': c['full_even_supports'] = [1]
            elif mutation == 'overlap_layer': c['exceptional_antipodal_lines'] = [[255, 0]]
            elif mutation == 'unknown': c['extra'] = 1
            else: c['dimension'] = True
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.generic(c)

    def test_limits_precede_exponential_or_pair_work(self):
        with patch.object(layers, '_rank', side_effect=AssertionError('rank reached')):
            with self.assertRaises(ValueError):
                self.generic(context=RunContext(Path('.'), {'max_codewords': 1}))
        c = small_certificate()
        c['full_even_supports'] += [mask(range(16, 24)), mask(range(24, 32))]
        c['expected_cardinality'] += 256
        with patch.object(layers, 'combinations', side_effect=AssertionError('pairs reached')):
            with self.assertRaises(ValueError):
                self.generic(c, RunContext(Path('.'), {'max_pairs': 1}))
        for budget in ({'max_dimension': 31}, {'max_materialized_points': 1}):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                self.generic(context=RunContext(Path('.'), budget))

    def test_operator_inputs_and_payload_types_are_strict(self):
        for operation in (layers.known32, layers.signed_supports):
            with self.assertRaises(ValueError): operation([Artifact('binary_code', {})], {}, self.context)
            with self.assertRaises(ValueError): operation([], {'unknown': 1}, self.context)
        a = self.generic()
        for mutation in ('kind', 'backend', 'provenance', 'extra_data', 'fact_float'):
            bad = deepcopy(a)
            if mutation == 'kind': bad.kind = 'binary_code'
            elif mutation == 'backend': bad.data['backend'] = 'arbitrary'
            elif mutation == 'provenance': bad.data['provenance'] = {'classification': 'new_record'}
            elif mutation == 'extra_data': bad.data['arbitrary'] = 1
            else: bad.facts['dimension'] = 32.0
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): layers.verify_artifact(bad)


if __name__ == '__main__':
    unittest.main()
