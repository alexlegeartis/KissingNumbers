from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.binary_coset_discovery import search,certify
from discovery.runner import registry

class BinaryCosetDiscoveryTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_subcode_coset_expansion(self):
  p={'length':4,'parent_generator':[3,12],'codimension':1,'minimum_distance':2,'trials':3,'seconds':5,'max_candidates':16}
  a=search([],p,self.ctx);self.assertIsNotNone(a.data['proposal']);self.assertEqual(verify_artifact(a,self.ctx).status,'unverified');b=certify([a],{},self.ctx);self.assertGreater(b.facts['cardinality'],4);self.assertGreaterEqual(b.facts['min_distance'],2);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
 def test_invalid_parent_and_budget(self):
  p={'length':4,'parent_generator':[3,3],'codimension':1,'minimum_distance':2,'seconds':2,'max_candidates':16}
  with self.assertRaises(ValueError):search([],p,self.ctx)
  p['parent_generator']=[3,12];p['max_candidates']=2000
  with self.assertRaises(ValueError):search([],p,self.ctx)
 def test_coordinate_equalities(self):
  p={'length':4,'parent_generator':[3,12],'codimension':1,'minimum_distance':2,'coordinate_equalities':[[0,2]],'seconds':5,'max_candidates':16}
  a=search([],p,self.ctx);self.assertIsNotNone(a.data['proposal']);self.assertGreater(certify([a],{},self.ctx).facts['cardinality'],4)
  p['coordinate_equalities']=[[0,1]]
  with self.assertRaises(ValueError):search([],p,self.ctx)
 def test_registry(self):
  r,_=registry();self.assertIn('search.binary_subcode_cosets',r);self.assertIn('codes.certify_binary_subcode_cosets',r)
if __name__=='__main__':unittest.main()
