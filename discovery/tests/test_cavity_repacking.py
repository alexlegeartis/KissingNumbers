from copy import deepcopy
from pathlib import Path
import unittest
from discovery.core import RunContext,Artifact
from discovery.runner import registry
from discovery.certification import verify_artifact
from discovery.engines.cavity_repacking import search_cavity,certify_cavity,verify_exact
class CavityTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));self.p={'source':[[1,0],[-1,0],[0,1],[0,-1]],'remove_indices':[2],'initial_replacements':[[0,1]],'attempts':1,'max_iterations':10,'seconds':4,'scale':100}
 def test_repacking_search_and_certification(self):
  self.assertIn('search.cavity_repacking',registry()[0]);r=search_cavity([],self.p,self.ctx);self.assertEqual(r.verification.status,'unverified');a=certify_cavity([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],4);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_exact_checks_fixed_fixed_pairs(self):
  data={'backend':'exact_cavity_repacking','source':[[1,0],[1,0],[0,1]],'metric':[1,1],'remove_indices':[],'replacement_vectors':[[0,-1]]}
  with self.assertRaisesRegex(ValueError,'separation'):verify_exact(Artifact('implicit_spherical_code',data),self.ctx)
 def test_unequal_norms_and_weighted_metric(self):
  data={'backend':'exact_cavity_repacking','source':[[2,0],[-1,0]],'metric':[1,2],'remove_indices':[],'replacement_vectors':[[0,1000000]]};a=Artifact('implicit_spherical_code',data);self.assertEqual(verify_exact(a,self.ctx).details['facts']['cardinality'],3)
 def test_invalid_remove_and_zero_seed(self):
  p=deepcopy(self.p);p['remove_indices']=[2,2]
  with self.assertRaisesRegex(ValueError,'removal'):search_cavity([],p,self.ctx)
  p=deepcopy(self.p);p['initial_replacements']=[[0,0]]
  with self.assertRaisesRegex(ValueError,'zero'):search_cavity([],p,self.ctx)
if __name__=='__main__':unittest.main()
