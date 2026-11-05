import unittest
import numpy as np
from discovery.engines.exchange_components import search_components

class ExchangeComponentTests(unittest.TestCase):
 def test_random_small_graphs_against_exhaustion(self):
  rng=np.random.default_rng(22)
  for trial in range(20):
   m,n=5,4;B=rng.random((m,n))<.4;C=np.triu(rng.random((m,m))<.3,1);a,b=np.where(B);i,j=np.where(C)
   picked,removed,status=search_components(m,n,a,b,i,j,10);exists=False
   for mask in range(1,1<<m):
    S=np.array([(mask>>k)&1 for k in range(m)],bool)
    if C[np.ix_(S,S)].any():continue
    exists |= S.sum()>B[S].any(0).sum()
   self.assertEqual(bool(picked),exists,msg=f'trial {trial}: {status}')
   if picked:
    self.assertGreater(len(picked),len(removed));self.assertFalse(C[np.ix_(picked,picked)].any());self.assertTrue(set(np.flatnonzero(B[picked].any(0)))<=set(removed))
if __name__=='__main__':unittest.main()
