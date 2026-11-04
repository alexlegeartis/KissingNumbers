from itertools import combinations
from pathlib import Path
import unittest

import numpy as np

from discovery.core import Artifact, RunContext, digest
from discovery.certification import verify_artifact
from discovery.engines.partial_orbits import (partial_sign_orbits, phase_kernel,
                                              screen_orbits, sign_group)
from discovery.runner import registry


class PartialOrbitTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path('.'))

    def search(self, source, representatives, generators=(), weights=None, norm=25, **kw):
        return partial_sign_orbits([Artifact('spherical_code',{'coordinates':source,
            'weights':weights or [1]*len(source[0]),'norm_squared':norm})],
            {'representatives':representatives,'generators':list(generators),'seconds':2,**kw},self.context)

    def test_breaks_internally_conflicting_orbit(self):
        result = self.search([[2,1],[-2,1]],[[1,-2]],[[1,0]],norm=5)
        self.assertEqual(result.facts['cardinality'],3)
        self.assertEqual(result.data['search']['selected_pool_points'],2)
        self.assertEqual(len(result.data['search']['selected_candidate_indices']),1)
        self.assertFalse(result.data['search']['whole_orbits_required'])
        self.assertEqual(result.verification.details['data_sha256'],digest(result.data))
        self.assertEqual(verify_artifact(result,self.context).status,'exact')

    def test_stabilizer_and_duplicate_orbits_not_overcounted(self):
        result = self.search([[3,4],[-3,4]],[[0,-5],[0,-5],[3,4]],[[1,0]])
        self.assertEqual(result.facts['cardinality'],3)
        self.assertEqual(result.data['search']['selected_pool_points'],1)
        self.assertEqual(result.data['search']['distinct_new_orbits'],1)

    def test_coupled_two_for_three_all_old_points_removable(self):
        result = self.search([[65,0],[-16,63],[-60,-25]],
                             [[25,60],[-56,33],[-25,-60]],norm=4225)
        self.assertEqual(result.facts['cardinality'],4)
        self.assertEqual(len(result.data['search']['removed_source_indices']),2)
        self.assertEqual(len(result.data['search']['selected_candidate_indices']),3)

    def test_phase_kernel_equals_all_direct_signed_products(self):
        reps = [[3,4,0],[0,3,4],[4,0,3]]
        group = sign_group([[1,0,1],[0,1,1]],3,self.context)
        kernel,_ = phase_kernel(reps,group,[1,1,1],25,self.context)
        for a,b in combinations(range(3),2):
            for p in range(4):
                for q in range(4):
                    exact = sum(int(x*y*s*t) for x,y,s,t in zip(reps[a],reps[b],group[p],group[q]))
                    self.assertEqual(bool(int(kernel[a,b]) >> (p^q) & 1),2*exact>25)

    def test_phase_mask_handles_the_high_unsigned_bit(self):
        generators=[[int(i==j) for i in range(6)] for j in range(6)]
        group=sign_group(generators,6,self.context)
        kernel,_=phase_kernel([[1]*6,[-1]*6],group,[1]*6,6,self.context)
        self.assertTrue(int(kernel[0,1]) >> 63 & 1)
        self.assertFalse(int(kernel[0,1]) & 1)

    def test_blocker_transport_matches_direct_geometry(self):
        source = [[3,4],[-3,4]]
        active,blockers,group,perm,_ = screen_orbits(source,[1,1],25,[[4,3]],[[1,0]],9,self.context)
        for p,signs in enumerate(group):
            direct = {i for i,row in enumerate(source)
                      if 2*sum(int(x*y*s) for x,y,s in zip(active[0],row,signs))>25}
            self.assertEqual(direct,set(map(int,perm[p,blockers[0]])))

    def test_limits_are_reported(self):
        result = self.search([[2,1],[-2,1]],[[1,-2]],[[1,0]],norm=5,candidate_limit=1)
        self.assertTrue(result.data['search']['pool_truncated'])
        self.assertEqual(result.data['search']['eligible_points'],2)

    def test_point_reservoir_explores_reproducible_different_frontiers(self):
        results=[]
        for seed in range(6):
            result=self.search([[5,0]],[[-5,0],[0,5],[0,-5]],candidate_limit=1,pool_seed=seed)
            self.assertEqual(result.facts['cardinality'],2)
            results.append(result.data['coordinates'])
        self.assertGreater(len({tuple(map(tuple,r)) for r in results}),1)
        again=self.search([[5,0]],[[-5,0],[0,5],[0,-5]],candidate_limit=1,pool_seed=0)
        self.assertEqual(again.data['coordinates'],results[0])

    def test_weighted_metric_and_large_integer_fallback(self):
        result = self.search([[2,0],[-2,0]],[[0,1]],[[1,0]],weights=[1,4],norm=4)
        self.assertEqual(result.facts['cardinality'],3)
        scale=2**80
        result = self.search([[scale,0],[-scale,0]],[[0,scale]],[[1,0]],norm=scale**2)
        self.assertEqual(result.facts['cardinality'],3)
        self.assertEqual(result.data['search']['kernel_integer_backend'],'Python integers')

    def test_rejects_nonclosed_source_bad_norm_and_dependent_generators(self):
        with self.assertRaisesRegex(ValueError,'not closed'):
            self.search([[5,0]],[[0,5]],[[1,0]])
        with self.assertRaisesRegex(ValueError,'wrong norm'):
            self.search([[5,0]],[[0,4]])
        for generators in ([[1,0],[1,0]],[[True,0]],[[1]]):
            with self.assertRaises(ValueError):
                self.search([[5,0],[-5,0]],[[0,5]],generators)

    def test_registered(self):
        operators,kinds=registry()
        self.assertIn('search.partial_sign_orbits',operators)
        self.assertEqual(kinds['search.partial_sign_orbits'],('spherical_code',))

    def test_budget_rejection_and_empty_catalogue(self):
        result = self.search([[5,0],[-5,0]],[],[[1,0]])
        self.assertEqual(result.facts['cardinality'],2)
        self.context=RunContext(Path('.'),{'max_materialized_points':2})
        with self.assertRaisesRegex(ValueError,'max_materialized_points'):
            self.search([[5,0],[-5,0]],[[0,5]],[[1,0]])
        with self.assertRaises(ValueError):
            self.search([[5,0],[-5,0]],[],[[1,0]],seconds=float('nan'))


if __name__=='__main__':
    unittest.main()
