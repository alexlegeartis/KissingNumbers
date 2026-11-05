from copy import deepcopy
from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.runner import registry
from discovery.engines.mixed_group_search import mixed_search
from discovery.engines.group_seed_search import certify_group_seeds
from discovery.certification import verify_artifact

class MixedGroupTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));swap={'matrix':[[0,1],[1,0]],'denominator':1};neg={'matrix':[[-1,0],[0,-1]],'denominator':1};nswap={'matrix':[[0,-1],[-1,0]],'denominator':1}
  self.p={'dimension':2,'generators':[swap,neg],'orbit_specs':[{'stabilizer_generators':[swap]},{'stabilizer_generators':[nswap]}],'initial_seeds':[[1,1],[1,-1]],'keep_initial':True,'attempts':1,'max_iterations':8,'seconds':5,'scale':100}
 def test_mixed_stabilizers_promote_exactly(self):
  self.assertIn('search.mixed_stabilizer_orbits',registry()[0]);r=mixed_search([],self.p,self.ctx);self.assertEqual(r.data['result']['stabilizer_orders'],[2,2]);self.assertEqual(r.data['result']['orbit_cardinality_upper_forecast'],4)
  a=certify_group_seeds([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],4);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_initial_seed_must_obey_its_stabilizer(self):
  p=deepcopy(self.p);p['initial_seeds'][0]=[1,0]
  with self.assertRaisesRegex(ValueError,'not fixed'):mixed_search([],p,self.ctx)
 def test_subgroup_must_be_inside_group(self):
  p=deepcopy(self.p);p['orbit_specs'][0]['stabilizer_generators']=[{'matrix':[[-1,0],[0,1]],'denominator':1}]
  with self.assertRaisesRegex(ValueError,'not a subgroup'):mixed_search([],p,self.ctx)
 def test_keeping_initial_never_bypasses_exact_promotion(self):
  p=deepcopy(self.p);p['initial_seeds'][1]=[1,1];p['orbit_specs'][1]=deepcopy(p['orbit_specs'][0]);r=mixed_search([],p,self.ctx)
  with self.assertRaises(ValueError):certify_group_seeds([r],{},self.ctx)
if __name__=='__main__':unittest.main()
