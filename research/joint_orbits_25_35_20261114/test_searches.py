"""Small exact geometry/matching regressions for the combined campaign."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('union',HERE/'union_certificates.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class MatchingGeometry(unittest.TestCase):
 def setUp(self):m.CACHE.clear()
 def test_shared_blocker_produces_a_real_gain(self):
  A=np.array([[8,8,0,0,0]],np.int16);B=np.zeros_like(A)
  C=np.array([[8,4,4,4,4],[4,8,-4,-4,-4]],np.int16)
  RA,RB,r=m.union(A,B,C,np.zeros_like(C),3)
  self.assertEqual(r['gain'],1);self.assertEqual((r['added'],r['removed']),(2,1))
  self.assertEqual(len(RA),2)
  G=RA.astype(np.int64)@RA.astype(np.int64).T
  self.assertTrue(np.array_equal(np.diag(G),[128,128]));self.assertLessEqual(G[0,1],64)
 def test_complete_matching_preserves_the_incumbent(self):
  A=np.array([[8,8,0,0,0]],np.int16);C=np.array([[8,4,4,4,4]],np.int16)
  RA,RB,r=m.union(A,np.zeros_like(A),C,np.zeros_like(C),3)
  self.assertEqual(r['gain'],0);self.assertTrue(np.array_equal(RA,A))
 def test_cached_and_uncached_graphs_agree(self):
  A=np.array([[8,8,0,0,0]],np.int16);C=np.array([[8,4,4,4,4]],np.int16)
  first=m.union(A,np.zeros_like(A),C,np.zeros_like(C),3)[2]
  second=m.union(A,np.zeros_like(A),C,np.zeros_like(C),3)[2]
  self.assertEqual(first,second)
 def test_four_sign_orbit_conflict_formula(self):
  import itertools
  a=np.array([8,4,4,4,4,0,0,0,4,4,0])
  b=np.array([4,8,-4,-4,0,0,0,0,4,0,4])
  actual=max(int(np.r_[s*a[:8],t*a[8:]]@b) for s,t in itertools.product((-1,1),repeat=2))
  self.assertEqual(actual,abs(int(a[:8]@b[:8]))+abs(int(a[8:]@b[8:])))
 def test_quadratic_sign_at_the_boundary(self):
  p=np.array([64,64,64,65,63]);q=np.array([0,1,-1,-1,1])
  self.assertEqual(m.above(p,q,3).tolist(),[False,True,False,False,True])
if __name__=='__main__':unittest.main()
