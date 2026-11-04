"""Independent count, norm, metric and resource controls for root formulas."""
from copy import deepcopy
from itertools import permutations
from pathlib import Path
import unittest
from unittest.mock import patch

from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines import lattices
from discovery.runner import output_artifact, run_program


class RootSystemTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def build(self, family, rank, context=None):
        return lattices.root_system([], {"family": family, "rank": rank}, context or self.context)

    def test_classical_counts_and_dimensions_1_to_8(self):
        for family, rank, count in (("A", 1, 2), ("A", 2, 6), ("A", 3, 12),
                                    ("D", 4, 24), ("D", 5, 40), ("E", 6, 72),
                                    ("E", 7, 126), ("E", 8, 240)):
            with self.subTest(family=family, rank=rank):
                artifact = self.build(family, rank)
                self.assertEqual(artifact.verification.status, "exact")
                self.assertEqual(artifact.facts["cardinality"], count)
                self.assertEqual(artifact.facts["dimension"], rank)
                self.assertEqual(artifact.facts["ambient_dimension"], rank)
                self.assertEqual(artifact.verification.details["pairs_checked"], count * (count - 1) // 2)
                self.assertEqual(verify_artifact(artifact, self.context).status, "exact")

    def test_helmert_metric_matches_original_A_root_gram(self):
        # Verify the full weighted Gram identity against ambient e_i-e_j roots,
        # independently of the spherical acceptance threshold or shared checker.
        rank = 5
        artifact = self.build("A", rank)
        ambient = []
        for first, second in permutations(range(rank + 1), 2):
            vector = [0] * (rank + 1)
            vector[first], vector[second] = 1, -1
            ambient.append(vector)
        coordinates, weights = artifact.data["coordinates"], artifact.data["weights"]
        scale = artifact.data["norm_squared"] // 2
        for i, first in enumerate(coordinates):
            for j, second in enumerate(coordinates):
                weighted = sum(w * x * y for w, x, y in zip(weights, first, second))
                original = sum(x * y for x, y in zip(ambient[i], ambient[j]))
                self.assertEqual(weighted, scale * original)

    def test_exceptional_sections_preserve_integer_metric(self):
        e8 = list(lattices._type_e8())
        self.assertEqual(len(e8), 240)
        self.assertEqual(len(set(map(tuple, e8))), 240)
        for rank in (6, 7):
            selected = [row for row in e8 if row[6] + row[7] == 0 and (rank == 7 or row[5] + row[7] == 0)]
            artifact = self.build("E", rank)
            weights = artifact.data["weights"]
            for original, projected in zip(selected, artifact.data["coordinates"]):
                self.assertEqual(sum(x * x for x in original), sum(w * x * x for w, x in zip(weights, projected)))
            self.assertEqual(len(selected), {6: 72, 7: 126}[rank])

    def test_invalid_requests_and_preallocation_budgets(self):
        for params in ({"family": "B", "rank": 3}, {"family": "D", "rank": 1},
                       {"family": "E", "rank": 5}, {"family": "A", "rank": True},
                       {"family": "A", "rank": 3, "typo": 1}):
            with self.subTest(params=params), self.assertRaises(ValueError):
                lattices.root_system([], params, self.context)
        for family, rank, budget in (("A", 96, {"max_pairs": 1}),
                                     ("D", 96, {"max_materialized_points": 10}),
                                     ("E", 8, {"max_dimension": 7})):
            with patch.object(lattices, "_type_a", side_effect=AssertionError("allocated")), patch.object(
                    lattices, "_type_d", side_effect=AssertionError("allocated")), patch.object(
                    lattices, "_type_e8", side_effect=AssertionError("allocated")):
                with self.subTest(family=family, budget=budget), self.assertRaises(ValueError):
                    self.build(family, rank, RunContext(Path("/tmp"), budget))

    def test_generic_families_have_no_rank_8_cutoff(self):
        for family, rank, count in (("A", 9, 90), ("D", 12, 264)):
            artifact = self.build(family, rank)
            self.assertEqual(artifact.facts["cardinality"], count)
            self.assertEqual(artifact.facts["dimension"], rank)

    def test_runner_integration_and_modified_geometry_rejection(self):
        program = {"schema_version": 1, "nodes": [{"id": "roots", "op": "lattice.root_system",
                   "params": {"family": "E", "rank": 6}}], "outputs": ["roots"]}
        artifact = output_artifact(run_program(program, self.context))
        self.assertEqual(artifact.facts["dimension"], 6)
        bad = deepcopy(artifact)
        bad.data["weights"][-1] = 1
        with self.assertRaises(ValueError):
            verify_artifact(bad, self.context)


if __name__ == "__main__":
    unittest.main()
