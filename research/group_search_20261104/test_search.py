import unittest
import numpy as np
from group_search import weil,orbit,phase_lift,maxcos
from verify import build10,build14,check

class SearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gs=weil(); cls.seed=np.eye(5,dtype=complex)[0]
        cls.A=orbit(cls.seed,cls.gs,projective=False)
    def test_40_lines_80_vectors(self):
        self.assertEqual(len(orbit(self.seed,self.gs)),40)
        self.assertEqual(len(self.A),80)
    def test_projective_gauge_is_not_vector_phase(self):
        # Prevent regression: arbitrary projective gauge destroys this code.
        lines=orbit(self.seed,self.gs)
        self.assertGreater(maxcos(phase_lift(lines,6)),.57)
        self.assertLessEqual(maxcos(phase_lift(self.A,3)),.5+1e-12)
    def test_compressed_cross_orbit_checks(self):
        # Full g-orbits, rotated by scalar phases, stay g-invariant.
        X=phase_lift(self.A,3)
        for theta in [0,.137,np.pi/6,np.pi/3,np.pi/2]:
            Y=X*np.exp(1j*theta)
            reps=phase_lift(self.A[:1],3)
            self.assertAlmostEqual(maxcos(reps,Y),maxcos(X,Y),places=11)
    def test_exact_controls_and_invalid_mutation(self):
        for dim,builder,N in [(10,build10,510),(14,build14,1932)]:
            C,w,n,_=builder();self.assertEqual(check(C,w,n,N,dim)['violations'],0)
            C[0]=C[1]
            with self.assertRaises(AssertionError):check(C,w,n,N,dim)
    def test_exact_inequality_rejection(self):
        C,w,n,_=build14()
        # Integer norm-preserving coordinate swap, independent of duplicates.
        for row in range(0,len(C),17):
            for i in range(14):
                for j in range(i):
                    v=C[row].copy();v[i],v[j]=v[j],v[i]
                    dots=C@v;dots[row]=-n
                    if not np.any(np.all(C==v,axis=1)) and max(dots)>n/2:
                        C[row]=v
                        with self.assertRaises(AssertionError):check(C,w,n,1932,14)
                        return
        self.fail('no violating norm-preserving mutation located')
if __name__=='__main__':unittest.main()
