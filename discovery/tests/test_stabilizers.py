"""Independent small-state enumeration and compressed theorem-backend negatives."""
from copy import deepcopy
from itertools import combinations
from pathlib import Path
import unittest

from discovery.core import RunContext
from discovery.engines import stabilizers


def explicit_small_sector(m, parity):
    """Enumerate flats independently as XOR-closed subsets, then quadratic signs."""
    n = 1 << m
    spaces = []
    for subset in range(1 << n):
        if not subset & 1:
            continue
        V = [x for x in range(n) if subset >> x & 1]
        if any(a ^ b not in V for a in V for b in V):
            continue
        r = len(V).bit_length()-1
        if len(V) != 1 << r or r % 2 != parity:
            continue
        basis, span = [], {0}
        for x in V:
            if x not in span:
                basis.append(x)
                span |= {y ^ x for y in tuple(span)}
        spaces.append((r, basis))
    points = []
    for r, basis in spaces:
        V = [0]
        for b in basis:
            V += [x ^ b for x in V]
        covered = set()
        for a in range(n):
            if a in covered:
                continue
            support = [a ^ v for v in V]
            covered.update(support)
            pairs = list(combinations(range(r), 2))
            for polynomial in range(1 << (1+r+len(pairs))):
                row = [0]*n
                for label, coordinate in enumerate(support):
                    value = polynomial & 1
                    for i in range(r):
                        value ^= ((polynomial >> (1+i)) & 1)*((label >> i) & 1)
                    for k, (i, j) in enumerate(pairs):
                        value ^= ((polynomial >> (1+r+k)) & 1)*((label >> i) & 1)*((label >> j) & 1)
                    row[coordinate] = 1 if value == 0 else -1
                points.append((tuple(row), 1 << r))
    return points


class StabilizerTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path('.'))

    def test_small_independent_all_pairs(self):
        for m in (1, 2, 3):
            for parity in (0, 1):
                with self.subTest(m=m, parity=parity):
                    artifact = stabilizers.real_stabilizer_sector([], {'m': m, 'parity': parity}, self.context)
                    points = explicit_small_sector(m, parity)
                    self.assertEqual(len(points), artifact.facts['cardinality'])
                    self.assertEqual(len(set(row for row, _ in points)), len(points))
                    for i, (a, anorm) in enumerate(points):
                        self.assertEqual(sum(x*x for x in a), anorm)
                        for b, bnorm in points[i+1:]:
                            dot = sum(x*y for x, y in zip(a, b))
                            if dot > 0:
                                self.assertLessEqual(4*dot*dot, anorm*bnorm)

    def test_generated_16_control_and_family_counts(self):
        expected = {4: 4320, 5: 146880, 6: 9694080}
        for m, count in expected.items():
            for parity in (0, 1):
                artifact = stabilizers.real_stabilizer_sector([], {'m': m, 'parity': parity}, self.context)
                self.assertEqual(artifact.facts['cardinality'], count)
                self.assertEqual(artifact.facts['dimension'], 1 << m)
                self.assertEqual(artifact.facts['ambient_dimension'], 1 << m)
                self.assertEqual(artifact.facts['rank_exact'], 1 << m)
                self.assertEqual(artifact.verification.details['materialized_spherical_points'], 0)
                self.assertNotIn('coordinates', artifact.data)
                self.assertEqual(stabilizers.verify_artifact(artifact).status, 'exact')

    def test_m_parity_and_budget_validation(self):
        for m, parity in ((0, 0), (-1, 0), (True, 0), (4.0, 0), (4, True), (4, -1), (4, 2)):
            with self.subTest(m=m, parity=parity), self.assertRaises(ValueError):
                stabilizers.real_stabilizer_sector([], {'m': m, 'parity': parity}, self.context)
        for m in (9, 10**30):
            with self.assertRaises(ValueError):
                stabilizers.real_stabilizer_sector([], {'m': m}, self.context)
        tiny = RunContext(Path('.'), {'max_rational_bits': 3})
        with self.assertRaises(ValueError):
            stabilizers.real_stabilizer_sector([], {'m': 4}, tiny)
        artifact = stabilizers.real_stabilizer_sector([], {'m': 8, 'parity': 1}, self.context)
        self.assertEqual(artifact.facts['dimension'], 256)

    def test_corrupt_m_count_proof_and_definition_rejected(self):
        artifact = stabilizers.real_stabilizer_sector([], {'m': 4}, self.context)
        mutations = [
            lambda a: a.data.update(m=5),
            lambda a: a.facts.update(cardinality=1),
            lambda a: a.facts.update(ambient_dimension=15),
            lambda a: a.data['proof']['layer_cardinalities'].__setitem__(0, 0),
            lambda a: a.data['proof'].update(distinct_pair_inner_product_upper_bound=[1, 1]),
            lambda a: a.data.update(sign_model='arbitrary_signs'),
            lambda a: a.data.update(parity=True),
        ]
        for mutate in mutations:
            damaged = deepcopy(artifact)
            mutate(damaged)
            with self.assertRaises(ValueError):
                stabilizers.verify_artifact(damaged)

    def test_registry_dispatch_and_recipe(self):
        import json
        from discovery.runner import run_program, output_artifact
        from discovery.certification import verify_artifact
        recipe = json.loads(Path('discovery/recipes/stabilizer16.json').read_text())
        receipt = run_program(recipe, self.context)
        artifact = output_artifact(receipt)
        self.assertEqual(artifact.facts['cardinality'], 4320)
        self.assertEqual(verify_artifact(artifact, self.context).status, 'exact')


if __name__ == '__main__':
    unittest.main()
