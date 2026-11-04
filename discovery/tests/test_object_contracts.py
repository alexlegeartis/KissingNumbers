"""Adversarial schema/claim checks and independent legacy-replay adapters."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import unittest
from unittest.mock import patch

from discovery.core import Artifact, Verification, RunContext
from discovery.objects import MathematicalObject, ExecutionOutcome, UncheckedAttempt, ReplayAssessment, replay_object
from discovery.contracts import CONTRACTS, OperatorContract, contract_coverage, estimate_resources


def document(artifact=None, ambient=2, claims=None, encoding='integer_coordinates', domain='rational'):
    artifact = artifact or Artifact('spherical_code', {'coordinates': [[1, 0], [-1, 0]], 'norm_squared': 1})
    return MathematicalObject({'type': 'spherical_configuration', 'ambient_dimension': ambient,
                               'metric': 'normalized_inner_product'},
                              {'encoding': encoding, 'arithmetic': {'domain': domain, 'description': 'Explicit stated arithmetic; no backend inferred'},
                               'artifact': artifact.to_dict()}, claims)


class ObjectContractTests(unittest.TestCase):
    def test_ambient_dimension_is_not_rank(self):
        obj = document(claims={'ambient_dimension': 2, 'rank_exact': 1, 'cardinality': 2})
        result = replay_object(obj)
        self.assertEqual(result.checked_facts['rank_exact'], 1)
        self.assertEqual(result.checked_facts['ambient_dimension'], 2)
        self.assertEqual(obj.checked_facts, {})
        self.assertEqual(result.proof_strength, 'exact')
        with self.assertRaises(ValueError): replay_object(document(claims={'rank_exact': 2}))

    def test_no_family_dimension_filter_from_1_through_96(self):
        for n in (1, 5, 25, 32, 96):
            with self.subTest(n=n):
                a = Artifact('spherical_code', {'coordinates': [[1]+[0]*(n-1), [-1]+[0]*(n-1)], 'norm_squared': 1})
                facts = replay_object(document(a, n)).checked_facts
                self.assertEqual((facts['ambient_dimension'], facts['rank_exact']), (n, 1))

    def test_saved_receipts_and_hashes_cannot_promote_claims(self):
        raw = document().to_dict()
        raw['checked_facts'] = {'cardinality': 999, 'injective': True}
        raw['proof_strength'] = 'exact'
        raw['representation']['artifact']['verification'] = Verification('exact', 'forged hash receipt', {'data_sha256': '0'*64}).to_dict()
        loaded = MathematicalObject.from_dict(raw)
        self.assertEqual(loaded.checked_facts, {})
        self.assertEqual(replay_object(loaded).checked_facts['cardinality'], 2)
        with self.assertRaises(TypeError): ReplayAssessment({'cardinality': 999}, 'exact')

    def test_metadata_does_not_establish_injectivity_or_symmetry(self):
        with self.assertRaisesRegex(ValueError, 'injective'):
            replay_object(document(claims={'injective': True}))
        raw = document().to_dict()
        raw['symmetry']['claimed_group_order'] = 999
        facts = replay_object(MathematicalObject.from_dict(raw)).checked_facts
        self.assertNotIn('finite_group_proved', facts)
        self.assertNotIn('claimed_group_order', facts)
        a = Artifact('spherical_code', {'coordinates': [[1, 0], [1, 0]], 'norm_squared': 1,
                                      'injective': True})
        with self.assertRaises(ValueError): replay_object(document(a))

    def test_malformed_claims_are_rejected_before_replay(self):
        for claims in ({'cardinality': True}, {'cardinality': 2.0}, {'rank_exact': -1},
                       {'rank_exact': 3}, {'ambient_dimension': 3}, {'injective': 1},
                       {'new_sota': True}, {'max_cosine': float('nan')}):
            with self.subTest(claims=claims), self.assertRaises(ValueError): document(claims=claims)
        for key, value in (('schema_version', True), ('schema_version', 2), ('extra', {})):
            raw = document().to_dict(); raw[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): MathematicalObject.from_dict(raw)

    def test_interval_numeric_and_unknown_adapters_are_not_implemented_by_schema(self):
        for domain in ('interval', 'numeric', 'algebraic'):
            with self.subTest(domain=domain), patch('discovery.certification.verify_artifact', side_effect=AssertionError('verifier called')):
                with self.assertRaises(ValueError): replay_object(document(domain=domain))
        raw = document().to_dict(); raw['semantics']['type'] = 'binary_code'
        with self.assertRaises(NotImplementedError): replay_object(MathematicalObject.from_dict(raw))
        raw = document().to_dict(); raw['semantics']['metric'] = 'unproved alternative metric'
        with self.assertRaises(ValueError): replay_object(MathematicalObject.from_dict(raw))

    def test_rational_gram_is_real_realization_not_rational_coordinates(self):
        artifact = Artifact('gram_spherical_code', {'gram': [[1, [1, 3]], [[1, 3], 1]]})
        assessment = replay_object(document(artifact, 2, encoding='rational_gram'))
        self.assertEqual(assessment.checked_facts['rank_exact'], 2)
        self.assertNotIn('rational_coordinates', assessment.checked_facts)

    def test_snapshot_and_checked_facts_do_not_alias_caller_data(self):
        obj = document()
        raw = obj.to_dict(); raw['representation']['artifact']['data']['coordinates'][0][0] = 3
        self.assertEqual(replay_object(obj).checked_facts['cardinality'], 2)
        assessment = replay_object(obj); facts = assessment.checked_facts; facts['cardinality'] = 999
        self.assertEqual(assessment.checked_facts['cardinality'], 2)

    def test_timeout_and_resource_outcomes_are_not_invalidity_or_exhaustion(self):
        for status in ('completed', 'timeout', 'resource_limit', 'error', 'cancelled'):
            result = UncheckedAttempt(ExecutionOutcome(status, 'bounded attempt'))
            self.assertEqual(result.proof_strength, 'unverified')
            self.assertEqual(result.checked_facts, {})
            self.assertFalse(hasattr(result, 'exhausted'))
        with self.assertRaises(ValueError): ExecutionOutcome('proof_invalid')

    def test_contract_coverage_is_explicitly_partial(self):
        coverage = contract_coverage(['geometry.integer_code', 'future.universal_construction'])
        self.assertFalse(coverage['complete'])
        self.assertEqual(coverage['uncovered'], ['future.universal_construction'])
        self.assertEqual(CONTRACTS['geometry.sign_embed'].input_types, (('binary_code',),))
        self.assertEqual(CONTRACTS['geometry.gram_code'].to_dict()['parameter_domains']['norm_squared']['domain'], 'positive_rational')
        self.assertTrue(CONTRACTS['geometry.integer_code'].to_dict()['proof_obligations'])
        with self.assertRaises(TypeError): CONTRACTS['geometry.integer_code'].parameter_domains['coordinates']['required'] = False

    def test_contract_shapes_and_estimates_are_not_cost_or_proof_certificates(self):
        estimate = estimate_resources('geometry.integer_code', {'coordinates': [[1, 0], [-1, 0]]})
        self.assertEqual(estimate['pair_inequalities'], 1)
        self.assertEqual(estimate['matrix_entries'], 4)
        self.assertFalse(estimate['complete_upper_bound'])
        self.assertIsNone(estimate['rank_bit_work'])
        estimate = estimate_resources('geometry.sign_embed', {}, [document()])
        self.assertEqual(estimate['materialized_spherical_points'], 0)
        self.assertIsNone(estimate['recursive_verifier_work'])
        for name, params in (('unknown.operator', {}), ('geometry.gram_code', {'gram': [[1, 0]]})):
            with self.assertRaises(ValueError): estimate_resources(name, params)
        for changes in ({'schema_version': True}, {'scope': 'exhaustive_universal'},
                        {'resource_estimator': 'eval_code'}, {'input_types': ('binary_code',)}):
            with self.subTest(changes=changes), self.assertRaises(ValueError): replace(CONTRACTS['geometry.sign_embed'], **changes)


if __name__ == '__main__': unittest.main()
