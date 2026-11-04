from pathlib import Path
from copy import deepcopy
import unittest

from discovery.core import RunContext
from discovery.engines.body_heads import search
from discovery.runner import run_program,output_artifact


class BodyHeadExchangeTests(unittest.TestCase):
    def setUp(self):
        self.context=RunContext(Path('.'))
        self.params={'bodies':[[1]],'body_norm_squared':1,'heads':[[-1],[0],[1]],
                     'incumbent_bodies':[0],'incumbent_heads':[[0]]}

    def test_one_for_two_exact_geometry_and_unverified_global_status(self):
        artifact=search([],self.params,self.context)
        self.assertEqual(artifact.kind,'exchange_proposal')
        self.assertEqual(artifact.facts['gain'],1)
        self.assertEqual(artifact.data['removed_indices'],[0])
        self.assertEqual(artifact.data['added_head_numerators'],[[-1],[1]])
        self.assertEqual(artifact.verification.status,'unverified')
        self.assertFalse(artifact.facts['new_sota'])

    def test_direct_insertion_and_empty_search(self):
        p=deepcopy(self.params);p['heads']=[[2]]
        a=search([],p,self.context)
        self.assertEqual((a.facts['added_points'],a.facts['removed_points']),(1,0))
        p['heads']=[[1]]
        self.assertEqual(search([],p,self.context).facts['gain'],0)

    def test_fixed_families_and_pure_heads_gate_replacements(self):
        p=deepcopy(self.params);p['pure_heads']=[[1]]
        self.assertEqual(search([],p,self.context).facts['gain'],0)
        p.pop('pure_heads');p.update(fixed_body_norms=[0],fixed_heads=[[1]],
            fixed_double_body_bounds=[[0]])
        self.assertEqual(search([],p,self.context).facts['gain'],0)

    def test_runner_exposes_a_proposal_without_promoting_it(self):
        program={'schema_version':1,'nodes':[{'id':'exchange','op':'search.body_head_exchange',
                  'params':self.params}],'outputs':['exchange']}
        a=output_artifact(run_program(program,self.context))
        self.assertEqual(a.facts['gain'],1)
        self.assertEqual(a.verification.status,'unverified')

    def test_invalid_shapes_norms_and_budgets(self):
        for key,value in [('bodies',[[True]]),('bodies',[[2]]),('head_denominator',True),
                          ('incumbent_bodies',[True]),('heads',[[0],[0]]),
                          ('incumbent_heads',[[1],[2]]),('body_classes',[1]),
                          ('screen_epsilon',float('nan'))]:
            p=deepcopy(self.params);p[key]=value
            with self.assertRaises(ValueError):search([],p,self.context)
        with self.assertRaisesRegex(ValueError,'max_codewords'):
            search([],self.params,RunContext(Path('.'),{'max_codewords':1}))

    def test_scaled_body_products_and_large_rational_heads(self):
        p=deepcopy(self.params);scale=10**13
        p['bodies']=[[12]];p['body_squared_scale']=144
        p['heads']=[[-scale],[0],[scale]];p['head_denominator']=scale
        self.assertEqual(search([],p,self.context).facts['gain'],1)

    def test_two_for_three_is_a_coupled_triangle(self):
        p={'bodies':[[1]],'body_norm_squared':1,
           'heads':[[2,2,2,-1,-1],[2,2,-1,2,-1],[2,2,-1,-1,2]],
           'incumbent_bodies':[0,0],
           'incumbent_heads':[[3,0,0,0,0],[0,3,0,0,0]]}
        self.assertEqual(search([],p,self.context).facts['gain'],0)
        p['max_removed']=2
        a=search([],p,self.context)
        self.assertEqual(a.data['removed_indices'],[0,1])
        self.assertEqual((a.facts['added_points'],a.facts['gain']),(3,1))

    def test_replacement_limit_is_explicit(self):
        p=deepcopy(self.params);p['heads']=[[-1],[1]]
        p['replacement_test_limit']=1
        self.assertFalse(search([],p,self.context).data['search']['replacement_search_truncated'])


if __name__=='__main__':unittest.main()
