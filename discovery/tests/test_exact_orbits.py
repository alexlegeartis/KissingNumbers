from copy import deepcopy
from pathlib import Path
from unittest.mock import patch
import unittest

from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
from discovery.engines.exact_orbits import integer_orbit


ROT={'matrix':[[0,-1],[1,0]],'denominator':1}


class ExactOrbitTests(unittest.TestCase):
    def setUp(self):self.context=RunContext(Path('.'))

    def square(self):return integer_orbit([],{'seed':[1,0],'generators':[ROT]},self.context)

    def test_exact_closed_transitive_square(self):
        a=self.square();r=verify_artifact(a,self.context)
        self.assertEqual(a.facts['cardinality'],4)
        self.assertEqual(a.facts['rank_exact'],2)
        self.assertTrue(a.facts['finite_group_proved'])
        self.assertEqual(r.details['seed_pairs_checked'],3)
        self.assertEqual(r.details['unordered_pairs_covered'],6)

    def test_stabilized_lower_dimensional_orbit_does_not_prove_group_finite(self):
        g={'matrix':[[5,0,0],[0,3,-4],[0,4,3]],'denominator':5}
        a=integer_orbit([],{'seed':[1,0,0],'generators':[g]},self.context)
        self.assertEqual(a.facts['rank_exact'],1)
        self.assertFalse(a.facts['finite_group_proved'])
        a.facts['finite_group_proved']=True
        with self.assertRaisesRegex(ValueError,'finite_group_proved'):verify_artifact(a,self.context)

    def test_closed_union_of_two_orbits_is_rejected(self):
        data=deepcopy(self.square().data)
        data['coordinates']=[[5,0],[0,5],[-5,0],[0,-5],[3,4],[-4,3],[-3,-4],[4,-3]]
        data['norm_squared']=25
        with self.assertRaisesRegex(ValueError,'one transitive orbit'):
            verify_artifact(Artifact('spherical_code',data),self.context)

    def test_missing_member_wrong_norm_and_duplicate_reject(self):
        for change in ('missing','norm','duplicate'):
            a=self.square();a.facts={}
            if change=='missing':a.data['coordinates'].pop()
            elif change=='norm':a.data['norm_squared']=2
            else:a.data['coordinates'].append(a.data['coordinates'][0][:])
            with self.assertRaises(ValueError):verify_artifact(a,self.context)

    def test_nonorthogonal_generator_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'exactly orthogonal'):
            integer_orbit([],{'seed':[1,0],'generators':[{'matrix':[[2,0],[0,1]],'denominator':1}]},self.context)

    def test_violating_whole_orbit_and_integral_scale_are_scoped(self):
        g={'matrix':[[3,-4],[4,3]],'denominator':5}
        with self.assertRaisesRegex(ValueError,'partial orbits remain open'):
            integer_orbit([],{'seed':[5,0],'generators':[g]},self.context)
        with self.assertRaisesRegex(ValueError,'integral coordinate scale'):
            integer_orbit([],{'seed':[1,0],'generators':[g]},self.context)

    def test_budget_rejections_precede_graph_allocation(self):
        a=self.square()
        with patch('discovery.engines.exact_orbits.np.empty',side_effect=AssertionError('allocated')):
            with self.assertRaisesRegex(ValueError,'max_orbit_points'):
                verify_artifact(a,RunContext(Path('.'),{'max_orbit_points':3}))
        with self.assertRaisesRegex(ValueError,'max_orbit_points'):
            integer_orbit([],{'seed':[1,0],'generators':[ROT]},RunContext(Path('.'),{'max_orbit_points':2}))

    def test_large_integer_coordinates_and_generator_fallback(self):
        scale=2**40
        a=integer_orbit([],{'seed':[scale,0],'generators':[ROT]},self.context)
        self.assertEqual(a.facts['norm_squared'],scale**2)
        big=2**80
        a=integer_orbit([],{'seed':[1],'generators':[{'matrix':[[big]],'denominator':big}]},self.context)
        self.assertEqual(a.facts['cardinality'],1)
        self.assertEqual(verify_artifact(a,self.context).status,'exact')


if __name__=='__main__':unittest.main()
