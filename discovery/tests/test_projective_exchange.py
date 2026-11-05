from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.projective_exchange import search,certify
from discovery.runner import registry

class ProjectiveExchangeTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_unequal_norm_extension_and_exact_replay(self):
  p={'source':[[1,0,0],[-1,0,0],[0,1,0]],'candidates':[[0,-2,0],[2,0,0]],'seconds':3};a=search([],p,self.ctx);self.assertEqual(a.verification.status,'unverified');b=certify([a],{},self.ctx);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(b.facts['dimension'],2);self.assertEqual(b.facts['ambient_dimension'],3);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
 def test_no_gain_returns_verified_source(self):
  p={'source':[[1,0],[-1,0],[0,1],[0,-1]],'candidates':[[2,0]],'seconds':3};a=search([],p,self.ctx);self.assertEqual(a.data['selected_indices'],[]);b=certify([a],{},self.ctx);self.assertEqual(b.facts['cardinality'],4)
 def test_invalid_source_and_candidate_rejected(self):
  with self.assertRaises(ValueError):search([],{'source':[[1,0],[1,0]],'candidates':[[0,1]],'seconds':1},self.ctx)
  with self.assertRaises(ValueError):search([],{'source':[[1,0]],'candidates':[[0,0]],'seconds':1},self.ctx)
 def test_proposal_replay_and_tampering(self):
  p={'source':[[1,0,0],[-1,0,0],[0,1,0]],'candidates':[[0,-2,0]],'seconds':3}
  a=search([],p,self.ctx)
  self.assertEqual(verify_artifact(a,self.ctx).status,'unverified')
  a.data['minimum_gain']=2
  with self.assertRaises(ValueError):certify([a],{},self.ctx)
 def test_registry(self):
  r,_=registry();self.assertIn('search.projective_exchange',r);self.assertIn('geometry.certify_projective_exchange',r)
if __name__=='__main__':unittest.main()
