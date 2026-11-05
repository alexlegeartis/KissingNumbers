from copy import deepcopy
from pathlib import Path
import unittest
from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
from discovery.engines.group_seed_search import group_seed_search,certify_group_seeds,verify_multi_artifact
from discovery.runner import registry

class GroupSeedTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));self.generators=[{'matrix':[[-1,0],[0,-1]],'denominator':1}]
 def test_numerical_search_and_exact_promotion(self):
  ops,kinds=registry();self.assertIn('search.group_orbit_seeds',ops);self.assertIn('groups.certify_multi_orbit',ops)
  p={'dimension':2,'orbit_count':2,'generators':self.generators,'attempts':1,'max_iterations':10,'seconds':4,'initial_seeds':[[1,0],[0,1]],'scale':100}
  proposal=group_seed_search([],p,self.ctx);self.assertEqual(proposal.verification.status,'unverified');self.assertLess(proposal.data['result']['best_max_correlation'],.001)
  a=certify_group_seeds([proposal],{},self.ctx);self.assertEqual(a.facts['cardinality'],4);self.assertEqual(a.facts['dimension'],2);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_unequal_integer_norms_are_checked_exactly(self):
  a=Artifact('implicit_spherical_code',{'backend':'exact_multi_group_orbits','generators':self.generators,'seed_vectors':[[2,0],[0,7]]})
  r=verify_multi_artifact(a,self.ctx);self.assertEqual(r.status,'exact');self.assertEqual(r.details['facts']['cardinality'],4)
  bad=deepcopy(a);bad.data['seed_vectors'][1]=[2,1]
  with self.assertRaises(ValueError):verify_multi_artifact(bad,self.ctx)
 def test_duplicate_orbit_rejected(self):
  a=Artifact('implicit_spherical_code',{'backend':'exact_multi_group_orbits','generators':self.generators,'seed_vectors':[[1,0],[2,0]]})
  with self.assertRaises(ValueError):verify_multi_artifact(a,self.ctx)
 def test_orthogonality_and_integral_group_lattice_checked(self):
  p={'dimension':2,'orbit_count':1,'generators':[{'matrix':[[2,0],[0,1]],'denominator':1}]}
  with self.assertRaises(ValueError):group_seed_search([],p,self.ctx)
 def test_half_integer_group_and_boundary_are_exact(self):
  H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
  a=Artifact('implicit_spherical_code',{'backend':'exact_multi_group_orbits','generators':[{'matrix':H,'denominator':2}],'seed_vectors':[[1,0,0,0]]})
  r=verify_multi_artifact(a,self.ctx);self.assertEqual(r.status,'exact');self.assertEqual(r.details['facts']['cardinality'],2)
 def test_common_denominator_closure_is_enforced(self):
  p={'dimension':2,'orbit_count':1,'generators':[{'matrix':[[3,-4],[4,3]],'denominator':5}]}
  with self.assertRaisesRegex(ValueError,'common-denominator'):group_seed_search([],p,self.ctx)
 def test_common_stabilizer_search_and_exact_promotion(self):
  swap={'matrix':[[0,1],[1,0]],'denominator':1};neg={'matrix':[[-1,0],[0,-1]],'denominator':1}
  p={'dimension':2,'orbit_count':1,'generators':[swap,neg],'stabilizer_generators':[swap], 'initial_seeds':[[1,1]],'attempts':1,'max_iterations':10,'seconds':4,'scale':100}
  proposal=group_seed_search([],p,self.ctx);r=proposal.data['result'];self.assertEqual(r['group_order'],4);self.assertEqual(r['stabilizer_order'],2);self.assertEqual(r['orbit_cardinality_upper_forecast'],2);self.assertEqual(r['seed_vectors'][0][0],r['seed_vectors'][0][1]);a=certify_group_seeds([proposal],{},self.ctx);self.assertEqual(a.facts['cardinality'],2);self.assertEqual(a.facts['dimension'],1)
 def test_stabilizer_must_be_in_group_and_have_fixed_vectors(self):
  swap={'matrix':[[0,1],[1,0]],'denominator':1};neg=self.generators[0];p={'dimension':2,'orbit_count':1,'generators':[neg],'stabilizer_generators':[swap]}
  with self.assertRaisesRegex(ValueError,'not a subgroup'):group_seed_search([],p,self.ctx)
  p['stabilizer_generators']=[neg]
  with self.assertRaisesRegex(ValueError,'no nonzero fixed'):group_seed_search([],p,self.ctx)
 def test_forged_search_seed_must_pass_promotion(self):
  p={'dimension':2,'orbit_count':2,'generators':self.generators,'attempts':1,'max_iterations':10,'seconds':4,'initial_seeds':[[1,0],[0,1]]}
  a=group_seed_search([],p,self.ctx);a.data['result']['seed_vectors']=[[1,0],[1,0]]
  with self.assertRaises(ValueError):certify_group_seeds([a],{},self.ctx)
if __name__=='__main__':unittest.main()
