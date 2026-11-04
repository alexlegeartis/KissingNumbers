import unittest
import numpy as np
from code_search import binary_factors,hensel,generator,gray,enumerate_small,translate_search

class CodeSearchTests(unittest.TestCase):
    def test_hensel_and_direct_binary_octacode(self):
        g,h=binary_factors(7)[0];g,h,steps=hensel(g,h,4)
        self.assertTrue(np.array_equal(np.convolve(g,h)%4,np.ones(7)))
        G=generator(g,7,4);d,hist=enumerate_small(G,4)
        self.assertEqual(d,6);self.assertEqual(sum(hist),256)
        digits=(np.arange(256)[:,None]//4**np.arange(4))%4
        table,_=gray(4);B=table[digits@G%4].reshape(256,16)
        self.assertEqual(len(np.unique(B,axis=0)),256)
        D=np.count_nonzero(B[:,None,:]!=B[None,:,:],axis=2);np.fill_diagonal(D,16)
        self.assertEqual(int(D.min()),6)
    def test_all_gray_symbol_pairs(self):
        for q in (4,8,16):
            T,w=gray(q)
            self.assertEqual(len(np.unique(T,axis=0)),q)
            for a in range(q):
                for b in range(q):self.assertEqual(int(np.count_nonzero(T[a]!=T[b])),int(w[(a-b)%q]))
    def test_dgls_translate_discovery(self):
        T,_=gray(8);result=translate_search(T,24,24)
        self.assertEqual(result['multiplier_certified'],2)
        self.assertEqual(result['cross_distance_lower_bound'],24)
        self.assertEqual(result['local_quotient_cosets'],2)
    def test_invalid_factor_rejected(self):
        g,h=binary_factors(7)[0];g[0]^=1
        with self.assertRaises(AssertionError):hensel(g,h,8)
    def test_distance_objective_changes_assembly(self):
        T,_=gray(8)
        self.assertEqual(translate_search(T,8,10)['multiplier_certified'],1)
        self.assertEqual(translate_search(T,8,8)['multiplier_certified'],2)
if __name__=='__main__':unittest.main()
