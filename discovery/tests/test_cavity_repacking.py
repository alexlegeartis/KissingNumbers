from copy import deepcopy
from pathlib import Path
import unittest
import numpy as np
from discovery.core import RunContext,Artifact
from discovery.runner import registry
from discovery.certification import verify_artifact
from discovery.engines.cavity_repacking import search_cavity,certify_cavity,verify_exact,_smooth_value_gradient,adaptive_cavity
class CavityTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));self.p={'source':[[1,0],[-1,0],[0,1],[0,-1]],'remove_indices':[2],'initial_replacements':[[0,1]],'attempts':1,'max_iterations':10,'seconds':4,'scale':100}
 def test_smooth_objective_and_gradient(self):
  X=np.array([[.2,.5],[.6,-.4],[-.3,.7]]);F=np.array([[.8,.1],[-.4,-.6]]);beta=4.
  value,grad,mx=_smooth_value_gradient(X,F,beta)
  terms=[X[i]@F[j] for i in range(3) for j in range(2)]+[X[i]@X[j] for i in range(3) for j in range(i+1,3)]
  expected=np.log(np.exp(beta*np.array(terms)).sum())/beta;self.assertAlmostEqual(value,expected,places=12);self.assertAlmostEqual(mx,max(terms),places=12)
  for i in range(3):
   for j in range(2):
    h=1e-6;Y=X.copy();Y[i,j]+=h;a=_smooth_value_gradient(Y,F,beta)[0];Y[i,j]-=2*h;b=_smooth_value_gradient(Y,F,beta)[0];self.assertAlmostEqual(grad[i,j],(a-b)/(2*h),places=7)
  value,grad,mx=_smooth_value_gradient(X,np.empty((0,2)),beta);self.assertTrue(np.isfinite(value));self.assertEqual(grad.shape,X.shape)
 def test_repacking_search_and_certification(self):
  self.assertIn('search.cavity_repacking',registry()[0]);r=search_cavity([],self.p,self.ctx);self.assertEqual(r.verification.status,'unverified');a=certify_cavity([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],4);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_exact_checks_fixed_fixed_pairs(self):
  data={'backend':'exact_cavity_repacking','source':[[1,0],[1,0],[0,1]],'metric':[1,1],'remove_indices':[],'replacement_vectors':[[0,-1]]}
  with self.assertRaisesRegex(ValueError,'separation'):verify_exact(Artifact('implicit_spherical_code',data),self.ctx)
 def test_unequal_norms_and_weighted_metric(self):
  data={'backend':'exact_cavity_repacking','source':[[2,0],[-1,0]],'metric':[1,2],'remove_indices':[],'replacement_vectors':[[0,1000000]]};a=Artifact('implicit_spherical_code',data);self.assertEqual(verify_exact(a,self.ctx).details['facts']['cardinality'],3)
 def test_adaptive_expansion_and_exact_replay(self):
  p={'source':[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]],'seed_candidate':[0,0,1],'expansion_rounds':2,'add_per_round':1,'max_iterations':2,'seconds':5,'round_seconds':2,'scale':100}
  r=adaptive_cavity([],p,self.ctx);h=r.data['result']['adaptive_history'];self.assertEqual(len(h),2);self.assertEqual(h[0]['removed'],0);self.assertEqual(h[1]['removed'],1);a=certify_cavity([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],5)
 def test_irrational_norm_threshold_boundary(self):
  data={'backend':'exact_cavity_repacking','source':[[1,1,0]],'metric':[1,1,1],'remove_indices':[],'replacement_vectors':[[1,0,1]]};self.assertEqual(verify_exact(Artifact('implicit_spherical_code',data),self.ctx).status,'exact')
  data['replacement_vectors']=[[2,1,0]]
  with self.assertRaisesRegex(ValueError,'separation'):verify_exact(Artifact('implicit_spherical_code',data),self.ctx)
 def test_invalid_remove_and_zero_seed(self):
  p=deepcopy(self.p);p['remove_indices']=[2,2]
  with self.assertRaisesRegex(ValueError,'removal'):search_cavity([],p,self.ctx)
  p=deepcopy(self.p);p['initial_replacements']=[[0,0]]
  with self.assertRaisesRegex(ValueError,'zero'):search_cavity([],p,self.ctx)
if __name__=='__main__':unittest.main()
