from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.engines.rational_reflection_search import search,quadratic_search
from discovery.engines.projective_exchange import certify
from discovery.certification import verify_artifact
from discovery.runner import registry

class RationalReflectionSearchTests(unittest.TestCase):
 def test_real_search_and_certification(self):
  c=RunContext(Path('.'));p={'source':[[1,0,0],[-1,0,0],[0,1,0]],'axes':[[0,1,0]],'candidate_limit':10,'seconds':5};a=search([],p,c);b=certify([a],{},c);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(verify_artifact(b,c).status,'exact')
 def test_weighted_unequal_norms(self):
  c=RunContext(Path('.'));p={'source':[[1,0,0],[-1,0,0],[0,2,0]],'axes':[[0,1,0]],'metric':[1,2,3],'candidate_limit':10,'seconds':5};a=search([],p,c);b=certify([a],{},c);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(verify_artifact(b,c).status,'exact')
 def test_reject_zero_axis_and_no_new_candidates(self):
  c=RunContext(Path('.'));p={'source':[[1,0],[-1,0]],'axes':[[0,0]]}
  with self.assertRaises(ValueError):search([],p,c)
  p['axes']=[[1,0]]
  with self.assertRaises(ValueError):search([],p,c)
 def test_reject_invalid_source(self):
  c=RunContext(Path('.'));p={'source':[[1,0],[2,0]],'axes':[[1,1]]}
  with self.assertRaises(ValueError):search([],p,c)
 def test_quadratic_reflection_certification(self):
  from discovery.engines.quadratic_exchange import certify as qcertify
  c=RunContext(Path('.'));p={'source_a':[[1,0,0],[-1,0,0],[0,1,0]],'source_b':[[0,0,0]]*3,'axes_a':[[0,1,0]],'axes_b':[[0,0,0]],'candidate_limit':10,'seconds':5};a=quadratic_search([],p,c);b=qcertify([a],{},c);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(verify_artifact(b,c).status,'exact')
 def test_quadratic_weighted_metric(self):
  from discovery.engines.rational_reflection_search import certify_weighted
  c=RunContext(Path('.'));p={'source_a':[[1,0,0],[-1,0,0],[0,1,0]],'source_b':[[0,0,0]]*3,'axes_a':[[0,1,0]],'axes_b':[[0,0,0]],'metric':[1,2,3],'candidate_limit':10,'seconds':5};a=quadratic_search([],p,c);b=certify_weighted([a],{},c);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(b.facts['dimension'],3);self.assertEqual(verify_artifact(b,c).status,'exact')
 def test_quadratic_safety_uses_actual_norm(self):
  from discovery.engines.quadratic_exchange import _points
  c=RunContext(Path('.'));a=[16,1]+[0]*18;b=[0]*20;p={'source_a':[a],'source_b':[b],'candidate_a':[[1,0]+[0]*18],'candidate_b':[b]};A,B,C,D=_points(p,c);self.assertEqual(A.shape,(1,20))
 def test_quadratic_overflow_rejected(self):
  from discovery.engines.quadratic_exchange import _points
  c=RunContext(Path('.'));a=[100000]+[0]*19;b=[0]*20;p={'source_a':[a],'source_b':[b],'candidate_a':[[1]+[0]*19],'candidate_b':[b]}
  with self.assertRaises(ValueError):_points(p,c)
 def test_weighted_tampering_rejected(self):
  from discovery.engines.rational_reflection_search import certify_weighted
  c=RunContext(Path('.'));p={'source_a':[[1,0,0],[-1,0,0],[0,1,0]],'source_b':[[0,0,0]]*3,'axes_a':[[0,1,0]],'axes_b':[[0,0,0]],'metric':[1,2,3],'candidate_limit':10,'seconds':5};a=quadratic_search([],p,c);a.data['proposal']['source_a'][0][2]=1
  with self.assertRaises(ValueError):certify_weighted([a],{},c)
 def test_registry(self):
  r,_=registry();self.assertIn('search.rational_reflection_exchange',r);self.assertIn('search.quadratic_reflection_exchange',r)
if __name__=='__main__':unittest.main()
