from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.support_recombination import search,certify
from discovery.runner import registry

class SupportRecombinationTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_new_family_and_replay(self):
  p={'dimension':12,'source':[(1<<8)-1],'permutations':[list(range(4,12))+list(range(4))],'seconds':5}
  a=search([],p,self.ctx);self.assertEqual(a.verification.status,'unverified');self.assertEqual(verify_artifact(a,self.ctx).status,'unverified');b=certify([a],{},self.ctx);self.assertEqual(b.facts['cardinality'],2);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
 def test_protected_family_and_tampering(self):
  p={'dimension':12,'source':[(1<<8)-1],'protected_supports':[(1<<8)-1],'permutations':[],'seconds':1}
  with self.assertRaises(ValueError):search([],p,self.ctx)
  p={'dimension':12,'source':[(1<<8)-1],'permutations':[list(range(4,12))+list(range(4))],'seconds':5};a=search([],p,self.ctx);a.data['selected']=[255,127|256]
  with self.assertRaises(ValueError):certify([a],{},self.ctx)
 def test_registry(self):
  r,_=registry();self.assertIn('search.support_recombination',r);self.assertIn('codes.certify_support_recombination',r)
if __name__=='__main__':unittest.main()
