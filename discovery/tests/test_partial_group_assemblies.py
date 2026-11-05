from copy import deepcopy
from pathlib import Path
import unittest
from discovery.core import RunContext,Artifact
from discovery.runner import registry
from discovery.certification import verify_artifact
from discovery.engines.partial_group_assemblies import partial_search,certify_partial,verify_exact
class PartialGroupTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));neg={'matrix':[[-1,0],[0,-1]],'denominator':1}
  self.p={'dimension':2,'orbit_specs':[{'generators':[neg]},{'generators':[]}], 'initial_seeds':[[2,0],[0,7]],'keep_initial':True,'attempts':1,'max_iterations':8,'seconds':5,'scale':100}
 def test_group_orbit_plus_independent_point(self):
  self.assertIn('search.partial_group_assemblies',registry()[0]);r=partial_search([],self.p,self.ctx);self.assertEqual(r.data['result']['orbit_cardinality_upper_forecast'],3);a=certify_partial([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],3);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_duplicate_across_different_groups_rejected(self):
  p=deepcopy(self.p);p['initial_seeds'][1]=[3,0];r=partial_search([],p,self.ctx)
  with self.assertRaisesRegex(ValueError,'separation'):certify_partial([r],{},self.ctx)
 def test_noncontained_stabilizer_rejected(self):
  p=deepcopy(self.p);p['orbit_specs'][0]['stabilizer_generators']=[{'matrix':[[0,1],[1,0]],'denominator':1}]
  with self.assertRaisesRegex(ValueError,'not a subgroup'):partial_search([],p,self.ctx)
 def test_seed_not_fixed_is_rejected(self):
  p=deepcopy(self.p);swap={'matrix':[[0,1],[1,0]],'denominator':1};p['orbit_specs'][0]={'generators':[swap],'stabilizer_generators':[swap]}
  with self.assertRaisesRegex(ValueError,'not fixed'):partial_search([],p,self.ctx)
 def test_exact_separation_on_boundary_and_unequal_norms(self):
  data={'backend':'exact_partial_group_orbits','groups':[[],[]],'seed_vectors':[[1,1,0],[1,0,1]],'metric':[1,1,1]};self.assertEqual(verify_exact(Artifact('implicit_spherical_code',data),self.ctx).status,'exact')
  data['seed_vectors']=[[2,1,0],[1,1,0]]
  with self.assertRaisesRegex(ValueError,'separation'):verify_exact(Artifact('implicit_spherical_code',data),self.ctx)
 def test_weighted_metric_and_metric_preservation(self):
  p=deepcopy(self.p);p['metric']=[1,2];r=partial_search([],p,self.ctx);a=certify_partial([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],3);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
  p['orbit_specs'][0]['generators']=[{'matrix':[[0,1],[1,0]],'denominator':1}]
  with self.assertRaisesRegex(ValueError,'preserve the metric'):partial_search([],p,self.ctx)
 def test_beta_schedule_validation(self):
  p=deepcopy(self.p);p['betas']=[0]
  with self.assertRaisesRegex(ValueError,'beta'):partial_search([],p,self.ctx)
if __name__=='__main__':unittest.main()
