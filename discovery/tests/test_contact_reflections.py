from pathlib import Path
from unittest.mock import patch
import unittest
import numpy as np
from discovery.core import Artifact,RunContext,Verification,digest
from discovery.certification import verify_artifact
from discovery.engines.contact_reflections import contact_reflections
from discovery.engines.construction_search import contact_exchange
from discovery.runner import run_program,output_artifact


class ReflectionTests(unittest.TestCase):
    def setUp(self):
        self.context=RunContext(Path('.'))
        self.rows=[[-1,-2,1],[-1,1,-2],[1,1,2],[-1,2,1],[2,-1,-1],[-2,-1,-1]]
    def source(self):return Artifact('spherical_code',{'coordinates':self.rows,'norm_squared':6})
    def test_strictly_beats_negative_contact_sums_and_checks_full_output(self):
        old=contact_exchange([self.source()],{'candidate_limit':256,'seconds':1},self.context)
        result=contact_reflections([self.source()],{'seconds':2},self.context)
        self.assertEqual(old.facts['cardinality'],6)
        self.assertEqual(result.facts['cardinality'],8)
        self.assertEqual(result.data['search']['removed_source_indices'],[2])
        self.assertEqual(result.verification.details['data_sha256'],digest(result.data))
        X=np.array(result.data['coordinates'],np.int64);G=X@X.T;np.fill_diagonal(G,-6)
        self.assertLessEqual(int(G.max()),3)
        self.assertEqual(verify_artifact(result,self.context).status,'exact')
    def test_weighted_hexagon_and_arbitrary_precision(self):
        for k in (1,2**80):
            source=Artifact('spherical_code',{'coordinates':[[2*k,0],[k,k],[-k,k],[-2*k,0],[-k,-k]],
                'weights':[1,3],'norm_squared':4*k*k})
            result=contact_reflections([source],{},self.context)
            self.assertEqual(result.facts['cardinality'],6)
        source=self.source();source.data['weights']=None
        self.assertEqual(contact_reflections([source],{},self.context).facts['cardinality'],8)
    def test_invalid_source_cannot_borrow_a_forged_receipt(self):
        bad=Artifact('spherical_code',{'coordinates':[[2,1,1],[2,1,-1]],'norm_squared':6})
        bad.verification=Verification('exact','forged cached receipt')
        with patch('research.contact_reflection_20261104.search',side_effect=AssertionError('called')):
            with self.assertRaises(ValueError):contact_reflections([bad],{},self.context)
    def test_resource_and_seed_preflights(self):
        for params in ({'candidate_limit':True},{'pair_limit':0},{'seconds':float('nan')},
                       {'seeds':[0,0]},{'seeds':[True]},{'blocker_limit':-1}):
            with self.assertRaises(ValueError):contact_reflections([self.source()],params,self.context)
        with patch('research.contact_reflection_20261104.search',side_effect=AssertionError('called')):
            with self.assertRaisesRegex(ValueError,'max_materialized_points'):
                contact_reflections([self.source()],{},RunContext(Path('.'),{'max_materialized_points':6}))
    def test_registered_program_retains_generation_provenance(self):
        program={'schema_version':1,'nodes':[{'id':'source','op':'geometry.integer_code','params':self.source().data},
          {'id':'search','op':'search.contact_reflections','inputs':['source'],'params':{'seconds':1}}],'outputs':['search']}
        result=output_artifact(run_program(program,self.context))
        self.assertEqual(result.facts['cardinality'],8)
        self.assertEqual(result.data['search']['pairs_examined'],30)
    def test_checked_matching_avoids_the_numerical_solver(self):
        rows=[[-1,2,-1],[-2,1,1],[1,-2,-1],[2,1,-1],[1,-1,2],[-1,-2,1],[-1,-1,-2],[1,2,1]]
        source=Artifact('spherical_code',{'coordinates':rows,'norm_squared':6})
        with patch('scipy.optimize.milp',side_effect=AssertionError('solver called')):
            result=contact_reflections([source],{},self.context)
        assigned=result.data['search']['matching_source_indices']
        self.assertEqual(len(set(assigned)),result.data['search']['candidates_after_filter'])
        self.assertEqual(result.facts['cardinality'],len(rows))
    def test_generation_cap_is_reported_with_exact_output(self):
        result=contact_reflections([self.source()],{'seeds':[0],'candidate_limit':1},self.context)
        self.assertEqual(result.data['search']['truncated'],'candidate_limit')
        self.assertEqual(result.data['search']['seed_count'],1)
        self.assertEqual(verify_artifact(result,self.context).status,'exact')
    def test_empty_seeds_is_an_explicit_exact_noop(self):
        result=contact_reflections([self.source()],{'seeds':[]},self.context)
        self.assertEqual(result.data['search']['pairs_examined'],0)
        self.assertEqual(result.data['coordinates'],self.rows)
    def test_wrong_kinds_and_unknown_parameters_are_rejected(self):
        for inputs in ([],[Artifact('binary_code',{})]):
            with self.assertRaises(ValueError):contact_reflections(inputs,{},self.context)
        with self.assertRaises(ValueError):contact_reflections([self.source()],{'expected_count':8},self.context)


if __name__=='__main__':unittest.main()
