from copy import deepcopy
from pathlib import Path
import unittest
import json
from pathlib import Path

from discovery.benchmarks.root_controls import root_recipe
from discovery.certification import verify_integer_spherical
from discovery.core import RunContext
from discovery.runner import output_artifact, run_program


class RootControlTests(unittest.TestCase):
    def test_committed_fixtures_regenerate_exactly(self):
        from discovery.benchmarks.root_controls import ROOT_NAMES
        root = Path(__file__).resolve().parents[1]
        for dimension, name in ROOT_NAMES.items():
            with self.subTest(dimension=dimension):
                self.assertEqual(json.loads((root / 'recipes' / f'roots_{name}.json').read_text()), root_recipe(dimension))

    def test_all_low_dimension_controls_through_runner(self):
        # Independently known ranks, cardinalities and inner products. These
        # expectations are not read from the generator's benchmark constants.
        expected = [(1, 2, -1), (2, 6, 1), (3, 12, 1), (4, 24, 1),
                    (5, 40, 1), (6, 72, 4), (7, 126, 4), (8, 240, 4)]
        for dimension, count, maximum in expected:
            with self.subTest(dimension=dimension):
                artifact = output_artifact(run_program(root_recipe(dimension), RunContext(Path('.'))))
                self.assertEqual(artifact.verification.status, 'exact')
                self.assertEqual(artifact.facts['dimension'], dimension)
                self.assertEqual(artifact.facts['cardinality'], count)
                self.assertEqual(artifact.facts['max_inner_product'], maximum)
                self.assertEqual(artifact.facts['pairs_checked'], count * (count - 1) // 2)

    def test_wrong_e8_sign_parity_fails_exact_separation(self):
        data = deepcopy(root_recipe(8)['nodes'][0]['params'])
        row = next(x for x in data['coordinates'] if all(abs(c) == 1 for c in x))
        row[0] *= -1  # Norm and row count survive; the kissing inequality fails.
        with self.assertRaisesRegex(ValueError, 'Pairwise kissing inequality violated'):
            verify_integer_spherical(data)

    def test_e8_respects_pair_budget(self):
        with self.assertRaisesRegex(ValueError, 'max_pairs budget exceeded'):
            run_program(root_recipe(8), RunContext(Path('.'), {'max_pairs': 28679}))


if __name__ == '__main__':
    unittest.main()
