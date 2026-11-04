from itertools import combinations
from pathlib import Path
from unittest.mock import patch
import unittest

from discovery.core import Artifact, RunContext
from discovery.certification import verify_artifact
from discovery.engines.construction_search import (integer_pool, support_packing,
                                                   contact_exchange, select_integer_pool)
from discovery.runner import run_program, output_artifact


class ConstructionSearchTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path('.'))

    def test_blind_nine_support_discovery(self):
        result = support_packing([], {'dimension': 9, 'cycles': [8,1], 'seconds': 2}, self.context)
        self.assertEqual(result.facts['cardinality'], 306)
        self.assertEqual(result.facts['rank_exact'], 9)
        self.assertFalse(result.data['search']['target_counts_used'])
        supports = result.data['search']['supports']
        self.assertTrue(all(len(set(a)&set(b)) <= 2 for a,b in combinations(supports,2)))
        self.assertEqual(verify_artifact(result, self.context).status, 'exact')

    def test_twelve_recovers_only_historical_840(self):
        result = support_packing([], {'dimension':12,'cycles':[2,2,4,4],'seconds':2},self.context)
        self.assertEqual(result.facts['cardinality'], 840)
        self.assertLess(result.facts['cardinality'], 841)

    def test_arbitrary_pool_trades_initial_points_and_deduplicates(self):
        # Three initial points, three points of a second copy; optimum is four.
        rows = [[-25,0],[15,20],[20,-15], [-15,20],[-7,-24],[24,-7],[-25,0]]
        result = select_integer_pool(rows, None, 625, self.context, initial_indices=[0,1,2], seconds=2)
        self.assertEqual(result.facts['cardinality'],4)
        self.assertEqual(result.data['search']['distinct_pool_points'],6)
        self.assertTrue(result.data['search']['all_initial_points_removable'])
        self.assertEqual(verify_artifact(result,self.context).status,'exact')

    def test_three_cohorts_not_two_copy_constraint(self):
        # Six equally spaced exact weighted integer vectors, supplied as three
        # mutually conflicting cohorts plus duplicates and an extra bad point.
        rows=[[2,0],[1,1],[-1,1],[-2,0],[-1,-1],[1,-1],[2,0]]
        result=integer_pool([],{'coordinates':rows,'weights':[1,3], 'norm_squared':4},self.context)
        self.assertEqual(result.facts['cardinality'],6)

    def test_complete_tiny_pool_matches_bruteforce(self):
        rows=[[5,0],[3,4],[0,5],[-3,4],[-5,0],[-3,-4],[0,-5],[3,-4]]
        optimum=max(len(c) for k in range(9) for c in combinations(rows,k)
                    if all(2*sum(a*b for a,b in zip(x,y))<=25 for x,y in combinations(c,2)))
        result=select_integer_pool(rows,None,25,self.context,seconds=2)
        self.assertEqual(result.facts['cardinality'],optimum)

    def test_contact_generation_and_runner(self):
        program={'schema_version':1,'nodes':[
            {'id':'roots','op':'lattice.root_system','params':{'family':'D','rank':4}},
            {'id':'exchange','op':'search.contact_exchange','inputs':['roots'],
             'params':{'candidate_limit':100,'seconds':1}}], 'outputs':['exchange']}
        result=output_artifact(run_program(program,self.context))
        self.assertEqual(result.facts['cardinality'],24)
        self.assertTrue(result.data['search']['generation_complete'])

    def test_preallocation_guards(self):
        with patch('discovery.engines.construction_search.combinations',side_effect=AssertionError('allocated')):
            with self.assertRaisesRegex(ValueError,'max_codewords'):
                support_packing([],{'dimension':96},self.context)
        with patch('numpy.array',side_effect=AssertionError('allocated')):
            with self.assertRaisesRegex(ValueError,'max_pairs'):
                integer_pool([],{'coordinates':[[2,0],[0,2]],'norm_squared':4},
                             RunContext(Path('.'),{'max_pairs': 0.5}))

    def test_bad_metrics_norms_cycles_and_initial_conflicts(self):
        for rows,weights,norm in [([[1,0]],None,2), ([[True,0]],None,1),
                                  ([[1,0]],[1,0],1), ([[1,0],[1]],None,1)]:
            with self.assertRaises(ValueError):
                select_integer_pool(rows,weights,norm,self.context)
        with self.assertRaisesRegex(ValueError,'Cycle lengths'):
            support_packing([],{'dimension':8,'cycles':[7]},self.context)
        with self.assertRaisesRegex(ValueError,'Initial subset'):
            select_integer_pool([[5,0],[3,4]],None,25,self.context,initial_indices=[0,1])

    def test_arbitrary_precision_pool_products(self):
        scale=2**40
        result=select_integer_pool([[scale,0],[-scale,0],[0,scale]],None,scale**2,self.context)
        self.assertEqual(result.facts['cardinality'],3)
        self.assertEqual(result.verification.details['integer_backend'],'Python integers')


if __name__=='__main__':
    unittest.main()
