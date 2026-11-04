from itertools import combinations
from pathlib import Path
import unittest

import numpy as np

from discovery.core import RunContext
from discovery.numerical_shapes import orbit_shape_objective


class NumericalShapeTests(unittest.TestCase):
    def setUp(self):
        self.context=RunContext(Path('.'))

    def test_analytic_gradient_against_finite_differences(self):
        shapes=[[3.,4.],[4.,-3.]]
        objective,audit=orbit_shape_objective([[1,1]],shapes,[[0,5],[0,-5]],self.context)
        x=np.array(shapes).ravel();loss,gradient=objective(x)
        self.assertGreater(loss,0)
        for i in range(len(x)):
            plus=x.copy();minus=x.copy();plus[i]+=1e-6;minus[i]-=1e-6
            finite=(objective(plus)[0]-objective(minus)[0])/2e-6
            self.assertAlmostEqual(gradient[i],finite,places=7)
        self.assertEqual(audit(x)['verification_status'],'numerical')
        self.assertFalse(audit(x)['new_sota'])

    def test_all_explicit_orbit_pairs_and_base_pairs_are_covered(self):
        shapes=np.array([[1.,2.,3.],[3.,1.,2.],[2.,3.,-1.]])
        base=[[2,0,0],[-2,0,0]]
        objective,audit=orbit_shape_objective([[1,1,1]],shapes.tolist(),base,self.context)
        report=audit(shapes.ravel())
        points=np.r_[shapes,-shapes];points/=np.linalg.norm(points,axis=1)[:,None]
        maximum=max(float(a@b) for a,b in combinations(points,2))
        self.assertAlmostEqual(report['max_pair_cosine'],maximum)
        self.assertAlmostEqual(report['max_base_cosine'],float(abs(points[:,0]).max()))

    def test_noninvariant_base_and_stabilizers_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'not invariant'):
            orbit_shape_objective([[1,1]],[[3.,4.]],[[5,0]],self.context)
        with self.assertRaisesRegex(ValueError,'free sign'):
            orbit_shape_objective([[1,0]],[[0.,4.]],[],self.context)
        with self.assertRaisesRegex(ValueError,'max_gram_entries'):
            orbit_shape_objective([[1,1]],[[3.,4.]],[],RunContext(Path('.'),{'max_gram_entries':1}))


if __name__=='__main__':unittest.main()
