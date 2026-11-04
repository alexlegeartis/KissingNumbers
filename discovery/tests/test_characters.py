from itertools import product
from pathlib import Path
from unittest.mock import patch
import unittest

import numpy as np

from discovery.core import RunContext
from discovery.engines.characters import (walsh_products, family_problem, check_assignment,
    optimize_families, materialize_assignment, search, assemble_orbits)
from discovery.runner import run_program, output_artifact


class CharacterFamilyTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path('.'))
        self.g = [[1,0,0], [0,1,0], [0,0,1]]
        self.shapes = [[1,2,2], [2,1,2], [2,2,1]]
        self.quotient = [[1,1,0], [0,1,1]]

    def test_walsh_matches_direct_weighted_products_and_bigints(self):
        for scale in (1, 2**40):
            shapes = [[scale*x for x in row] for row in self.shapes]
            # Walsh itself is useful for arbitrary norms; these permuted shapes
            # need equal weights for the common-norm construction contract.
            dots, weights, norm, columns = walsh_products(self.g, shapes, [3]*3, self.context)
            for a, b, m in product(range(3), range(3), range(8)):
                direct = sum(w*x*y*(-1 if (m & c).bit_count()%2 else 1)
                             for w,x,y,c in zip(weights, shapes[a], shapes[b], columns))
                self.assertEqual(int(dots[a,b,m]), direct)

    def test_joint_family_changes_beat_fixed_shape_and_match_exhaustion(self):
        p = family_problem(self.g, self.shapes, None, self.quotient, self.context)
        feasible = []
        for assignment in product(range(-1,3), repeat=4):
            try:
                count = check_assignment(p, list(assignment))
            except ValueError:
                continue
            feasible.append((count, assignment))
        optimum = max(count for count, a in feasible)
        fixed = max(count for count, a in feasible if set(a) <= {-1,0})
        self.assertEqual((optimum, fixed), (4,2))
        initial = next(list(a) for count,a in feasible if count==fixed and set(a)<={-1,0})
        assignment, diagnostic = optimize_families(p, self.context, 2, initial)
        self.assertEqual(check_assignment(p, assignment), optimum)
        self.assertEqual(len(materialize_assignment(p,self.shapes,assignment)),8)
        self.assertEqual(diagnostic['initial_families'],2)
        result = search([], {'generators':self.g,'shapes':self.shapes,'quotient':self.quotient,
                             'initial_assignment':initial,'seconds':2}, self.context)
        self.assertEqual(result.facts['cardinality'],8)
        self.assertEqual(result.verification.status,'exact')

    def test_quotient_kernel_detects_internal_conflicts(self):
        p=family_problem(self.g,self.shapes,None,[],self.context)
        self.assertFalse(any(p['admissible']))
        with self.assertRaisesRegex(ValueError,'internal conflict'):
            check_assignment(p,[0])

    def test_every_cross_table_entry_matches_materialized_families(self):
        p=family_problem(self.g,self.shapes,None,self.quotient,self.context)
        for q,r,a,b in product(range(4),range(4),range(3),range(3)):
            if q==r:continue
            first=[[(1-2*((m&c).bit_count()%2))*x for c,x in zip(p['columns'],self.shapes[a])]
                   for m in p['cosets'][q]]
            second=[[(1-2*((m&c).bit_count()%2))*x for c,x in zip(p['columns'],self.shapes[b])]
                    for m in p['cosets'][r]]
            bad=any(2*sum(x*y for x,y in zip(v,w))>9 for v in first for w in second)
            self.assertEqual(bool(p['forbidden'][q^r,a,b]),bad)

    def test_active_neighborhood_preserves_frozen_families(self):
        p=family_problem(self.g,self.shapes,None,self.quotient,self.context)
        a,_=optimize_families(p,self.context,2,initial=[0,-1,-1,-1],active=[1,2,3])
        self.assertEqual(a[0],0)
        self.assertGreaterEqual(check_assignment(p,a),1)
        with self.assertRaisesRegex(ValueError,'distinct coset'):
            optimize_families(p,self.context,2,active=[True])

    def test_runner_and_bad_schema(self):
        program={'schema_version':1,'nodes':[{'id':'sphere','op':'search.character_families',
                  'params':{'generators':self.g,'shapes':self.shapes,'quotient':self.quotient}}],
                 'outputs':['sphere']}
        self.assertEqual(output_artifact(run_program(program,self.context)).facts['cardinality'],8)
        for generators,shapes,quotient in [([[True,0,0]],self.shapes,[]),
                ([[1,0,0],[1,0,0]],self.shapes,[]), (self.g,[[1.,2,2]],[]),
                (self.g,self.shapes,[[1,1,0],[1,1,0]]), (self.g,[[1,0,0]],[])]:
            with self.assertRaises(ValueError):
                family_problem(generators,shapes,None,quotient,self.context)
        with self.assertRaisesRegex(ValueError,'common positive'):
            family_problem(self.g,[[1,2,2],[1,1,1]],None,[],self.context)

    def test_preallocation_guard_and_invalid_incumbent(self):
        with patch('numpy.zeros',side_effect=AssertionError('allocation')):
            with self.assertRaisesRegex(ValueError,'max_gram_entries'):
                walsh_products(self.g,self.shapes,None,
                               RunContext(Path('.'),{'max_gram_entries':10}))
        p=family_problem(self.g,self.shapes,None,self.quotient,self.context)
        with self.assertRaisesRegex(ValueError,'cross conflict'):
            optimize_families(p,self.context,1,[0]*4)

    def test_sparse_stabilizers_and_replacement_of_entire_incumbent(self):
        # The old axis orbit has two points and a nontrivial sign stabilizer.
        # The new four-point orbit conflicts with both old points. Joint
        # weighted assembly must delete the entire old family to improve.
        result=assemble_orbits([],{'generators':[[1,0],[0,1]],
            'shapes':[[25,0],[15,20]],'initial_shapes':[0]},self.context)
        self.assertEqual(result.facts['cardinality'],4)
        self.assertEqual(result.data['search']['selected_shapes'],[1])
        self.assertEqual(result.data['search']['initial_points'],2)
        self.assertTrue(result.data['search']['all_initial_families_removable'])

    def test_overlapping_orbits_cannot_double_count(self):
        result=assemble_orbits([],{'generators':[[1,0],[0,1]],
            'shapes':[[15,20],[-15,20]]},self.context)
        self.assertEqual(result.facts['cardinality'],4)
        with self.assertRaisesRegex(ValueError,'conflict or overlap'):
            assemble_orbits([],{'generators':[[1,0],[0,1]],
                'shapes':[[15,20],[-15,20]],'initial_shapes':[0,1]},self.context)


if __name__=='__main__':unittest.main()
