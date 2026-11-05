from pathlib import Path
import unittest
import numpy as np
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.quadratic_exchange import search,certify,_conflicts
from discovery.engines.signed_layers import compatible_dot
from discovery.runner import registry

class QuadraticExchangeTests(unittest.TestCase):
 def setUp(self): self.ctx=RunContext(Path('.'))
 def test_radical_exchange(self):
  p={'source_a':[[1,0],[-1,0],[0,1]],'source_b':[[0,0]]*3,'candidate_a':[[0,0]],'candidate_b':[[0,-1]],'seconds':3}
  a=search([],p,self.ctx);self.assertEqual(a.verification.status,'unverified');self.assertEqual(verify_artifact(a,self.ctx).status,'unverified');b=certify([a],{},self.ctx);self.assertEqual(b.facts['cardinality'],4);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
  a.data['minimum_gain']=2
  with self.assertRaises(ValueError):certify([a],{},self.ctx)
 def test_conflict_signs(self):
  rng=np.random.default_rng(4);a=rng.integers(-4,5,(30,4));b=rng.integers(-4,5,(30,4));c=rng.integers(-4,5,(25,4));d=rng.integers(-4,5,(25,4));n1=(a*a+3*b*b).sum(1);n2=(c*c+3*d*d).sum(1);K=_conflicts(a,b,c,d,n1,n2)
  for i in range(len(a)):
   for j in range(len(c)):
    p=int(a[i]@c[j]+3*b[i]@d[j]);q=int(a[i]@d[j]+b[i]@c[j]);self.assertEqual(K[i,j],not compatible_dot(p,q,int(n1[i]*n2[j])))
 def test_registry(self):
  r,_=registry();self.assertIn('search.quadratic_exchange',r);self.assertIn('geometry.certify_quadratic_exchange',r)
if __name__=='__main__':unittest.main()
