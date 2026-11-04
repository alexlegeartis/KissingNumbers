from copy import deepcopy
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from discovery.core import RunContext
from discovery.engines.construction_search import blocker_exchange
from discovery.runner import run_program,output_artifact

class BlockerExchangeTests(unittest.TestCase):
    def setUp(self):
        self.context=RunContext(Path('.'))
        self.params={'blockers':[[0],[0],[0]],'removal_costs':[2]}
    def test_shared_family_deletion_requires_three_candidates(self):
        result=blocker_exchange([],self.params,self.context)
        self.assertEqual(result.data['selected_candidate_indices'],[0,1,2])
        self.assertEqual(result.data['removed_source_indices'],[0])
        self.assertEqual(result.facts['gain'],1)
        self.assertEqual(result.verification.status,'unverified')
        self.assertFalse(result.facts['new_sota'])
    def test_candidate_conflicts_and_cliques_prevent_false_gains(self):
        for extra in ({'candidate_conflicts':[[0,1],[0,2],[1,2]]},{'cliques':[[0,1,2]]}):
            p=deepcopy(self.params);p.update(extra)
            self.assertEqual(blocker_exchange([],p,self.context).facts['gain'],0)
    def test_free_candidates_and_sparse_large_source_indices(self):
        p={'blockers':[[],[1000000],[1000000]],'candidate_weights':[1,2,2]}
        result=blocker_exchange([],p,self.context)
        self.assertEqual(result.facts['gain'],4)
        self.assertEqual(result.data['removed_source_indices'],[1000000])
    def test_malformed_oracles_and_budgets(self):
        for change in ({'blockers':[[True]]},{'blockers':[[[]]]},{'blockers':[[0,0]]},
                       {'candidate_weights':[1,2]}, {'candidate_weights':[True,1,1]},
                       {'candidate_conflicts':[[0,0]]},{'candidate_conflicts':[[0,3]]},
                       {'cliques':[[0,[]]]},{'removal_costs':[]},{'seconds':float('nan')}):
            p=deepcopy(self.params);p.update(change)
            with self.subTest(change=change),self.assertRaises(ValueError):
                blocker_exchange([],p,self.context)
        for budget in ({'max_codewords':3},{'max_gram_entries':5}):
            with self.assertRaises(ValueError):blocker_exchange([],self.params,RunContext(Path('.'),budget))
    def test_forged_solver_feasibility_is_rejected(self):
        class Forged:
            x=np.array([1,1,1,0]);status=0;message='forged';mip_dual_bound=0;mip_gap=0
        with patch('scipy.optimize.milp',return_value=Forged()):
            with self.assertRaisesRegex(ValueError,'exact oracle'):
                blocker_exchange([],self.params,self.context)
    def test_registered_dag_emits_a_proposal(self):
        program={'schema_version':1,'nodes':[{'id':'exchange','op':'selection.blocker_exchange',
            'params':self.params}],'outputs':['exchange']}
        self.assertEqual(output_artifact(run_program(program,self.context)).facts['gain'],1)

if __name__=='__main__':unittest.main()
