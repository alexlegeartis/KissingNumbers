from fractions import Fraction
from itertools import combinations
from pathlib import Path
from unittest.mock import patch
import unittest

import numpy as np

from discovery.core import Artifact,RunContext
from discovery.engines import groups
from discovery.group_constraints import analyze_orthogonal_actions
from discovery.planner import operator_catalogue
from discovery.runner import run_program,output_artifact


def sign_change_wedges(n):
    pairs=list(combinations(range(n),2));d=len(pairs)
    return [[[int(i==j)*(-1 if axis in pair else 1) for j in range(d)]
             for i,pair in enumerate(pairs)] for axis in range(n)]


class GroupConstraintTests(unittest.TestCase):
    def setUp(self):self.context=RunContext(Path('.'))

    def test_real_complex_type_needs_only_symmetric_average(self):
        # C4 on R2 is real irreducible; its singleton class sum is not scalar,
        # but the symmetric part is zero. No absolute-irreducibility shortcut.
        result=analyze_orthogonal_actions([1,0],[[[0,-1],[1,0]]])
        self.assertTrue(result['symmetric_average_is_scalar'])
        self.assertEqual(result['scalar_average'],0)
        self.assertFalse(result['whole_generated_orbit_ruled_out'])

    def test_stabilizers_prevent_false_positive_above_half(self):
        actions=[]
        for axis in range(3):
            actions.append([[int(i==j)*(-1 if i==axis else 1) for j in range(3)] for i in range(3)])
        # Standard sign changes in dimension3 average1/3. Use dimension9
        # to make the average7/9; e0 still has only two distinct orbit points.
        actions=[[[int(i==j)*(-1 if i==axis else 1) for j in range(9)] for i in range(9)] for axis in range(9)]
        result=analyze_orthogonal_actions([1]+[0]*8,actions)
        self.assertEqual(result['scalar_average'],[7,9])
        self.assertEqual(result['fixed_action_count'],8)
        self.assertEqual(result['necessary_fixed_actions_for_kissing'],5)
        self.assertFalse(result['whole_generated_orbit_ruled_out'])

    def test_b9_exterior_reflection_class_requires_stabilizer(self):
        result=analyze_orthogonal_actions([1]*36,sign_change_wedges(9))
        self.assertEqual(result['scalar_average'],[5,9])
        self.assertEqual(result['necessary_fixed_actions_for_kissing'],1)
        self.assertTrue(result['whole_generated_orbit_ruled_out'])
        self.assertFalse(result['partial_orbits_pruned'])
        self.assertEqual(len(result['violating_images']),9)

    def test_equality_tightness_and_stabilizer_exception(self):
        pairs=list(combinations(range(8),2))
        seed=[int(i%2==0 and j==i+1) for i,j in pairs]
        result=analyze_orthogonal_actions(seed,sign_change_wedges(8))
        self.assertEqual(result['scalar_average'],[1,2])
        self.assertTrue(result['equality_constraints_apply'])
        self.assertTrue(result['equality_constraints_satisfied'])
        rank_two=analyze_orthogonal_actions([1]+[0]*27,sign_change_wedges(8))
        self.assertEqual(rank_two['fixed_action_count'],6)
        self.assertFalse(rank_two['equality_constraints_apply'])
        self.assertFalse(rank_two['whole_generated_orbit_ruled_out'])

    def test_reducible_non_scalar_average_is_not_a_character_prune(self):
        result=analyze_orthogonal_actions([1,1],[[[-1,0],[0,1]]])
        self.assertFalse(result['symmetric_average_is_scalar'])
        self.assertNotIn('necessary_fixed_actions_for_kissing',result)

    def test_exact_rationals_and_preallocation_guards(self):
        result=analyze_orthogonal_actions([1,0],[[[[3,5],[-4,5]],[[4,5],[3,5]]]])
        self.assertTrue(result['whole_generated_orbit_ruled_out'])
        with patch('discovery.group_constraints.bounded_fraction',side_effect=AssertionError('allocated')):
            with self.assertRaisesRegex(ValueError,'max_dimension'):
                analyze_orthogonal_actions([1,0],[[[1,0],[0,1]]],RunContext(Path('.'),{'max_dimension':1}))
        for actions in ([[[2,0],[0,1]]],[[[True,0],[0,1]]],[[[1.,0],[0,1]]]):
            with self.assertRaises(ValueError):analyze_orthogonal_actions([1,0],actions)

    def test_exterior_functor_and_unitary_replay(self):
        a=np.array([[0,1,0],[-1,0,0],[0,0,1]],float)
        b=np.diag([1,-1,1])
        rep=groups.from_generators([],{'generators_real':[a.tolist(),b.tolist(),(a@b).tolist()]},self.context)
        lifted=groups.exterior_square([rep],{},self.context)
        value=groups._generators(lifted)
        self.assertTrue(np.array_equal(value[0]@value[1],value[2]))
        self.assertEqual(lifted.verification.status,'numerical')
        self.assertFalse(lifted.verification.details['finite_group_proved'])
        self.assertTrue(operator_catalogue()['groups.exterior_square']['planner_supported'])

    def test_exterior_preflight_before_source_allocation(self):
        source=Artifact('unitary_representation',{'generators_real':[np.eye(4).tolist()]})
        with patch.object(groups.np,'asarray',side_effect=AssertionError('allocated')):
            with self.assertRaisesRegex(ValueError,'max_dimension'):
                groups.exterior_square([source],{},RunContext(Path('.'),{'max_dimension':5}))
            with self.assertRaisesRegex(ValueError,'max_materialized_points'):
                groups.exterior_square([source],{},RunContext(Path('.'),{'max_materialized_points':71}))

    def test_exterior_through_shared_dag(self):
        program={'schema_version':1,'nodes':[
            {'id':'base','op':'groups.from_generators','params':{'generators_real':[[[-1,0,0],[0,1,0],[0,0,1]]]}},
            {'id':'wedge','op':'groups.exterior_square','inputs':['base']}],'outputs':['wedge']}
        artifact=output_artifact(run_program(program,self.context))
        self.assertEqual(artifact.facts['complex_dimension'],3)
        self.assertEqual(artifact.verification.status,'numerical')


if __name__=='__main__':unittest.main()
