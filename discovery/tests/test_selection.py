"""Exact coupled exchanges, checked matching bounds and independent geometry."""
from itertools import combinations, product
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from discovery.core import Artifact, RunContext, Verification
from discovery.certification import verify_artifact
from discovery.engines.geometry import integer_code
from discovery.engines import selection


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))
        self.source = integer_code([], {"coordinates": [[-5, 0], [3, 4], [4, -3]],
                                      "norm_squared": 25}, self.context)

    def test_real_exchange_removes_two_and_adds_three_with_exact_geometry(self):
        result = selection.isometric_exchange([self.source], {"normals": [[1, 2]]}, self.context)
        self.assertEqual(result.facts["cardinality"], 4)
        self.assertEqual(result.facts["rank_exact"], 2)
        certificate = result.verification.details["isometric_exchange"]
        self.assertEqual(certificate["left_indices"], [1])
        self.assertEqual(certificate["right_indices"], [0, 1, 2])
        self.assertEqual(len(certificate["matching"]), 2)
        self.assertEqual(certificate["conflict_edges"], 2)
        self.assertEqual(certificate["cross_pairs_checked"], 9)
        result.verification = Verification("unverified", "saved receipts have no authority")
        self.assertEqual(verify_artifact(result, self.context).status, "exact")
        result.data["coordinates"][0] = result.data["coordinates"][1][:]
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            verify_artifact(result, self.context)

    def test_matching_cover_bound_matches_exhaustive_small_graphs(self):
        # Independent brute force over every subset, for every 2-by-3 graph.
        for bits in product((0, 1), repeat=6):
            graph = [[v for v in range(3) if bits[3*u + v]] for u in range(2)]
            result = selection.maximum_independent_bipartite(graph, 3)
            best = 0
            for selected in product((False, True), repeat=5):
                if all(not (selected[u] and selected[2+v]) for u, row in enumerate(graph) for v in row):
                    best = max(best, sum(selected))
            self.assertEqual(len(result["left_indices"]) + len(result["right_indices"]), best)
        with patch("scipy.sparse.csgraph.maximum_bipartite_matching", return_value=np.array([-1])):
            with self.assertRaisesRegex(ValueError, "augmenting"):
                selection.maximum_independent_bipartite([[0]], 1)

    def test_overlap_identity_weighted_metric_and_matrix_word_agree(self):
        identity = selection.isometric_exchange([self.source], {"normals": []}, self.context)
        self.assertEqual(identity.facts["cardinality"], 3)
        self.assertEqual(identity.verification.details["isometric_exchange"]["overlap_points"], 3)
        weighted = integer_code([], {"coordinates": [[1, 0], [-1, 0]], "weights": [1, 2], "norm_squared": 1}, self.context)
        first = selection.isometric_exchange([weighted], {"normals": [[1, 1]]}, self.context)
        matrix = [[[1, 3], [-2, 3]], [[-4, 3], [-1, 3]]]
        second = selection.isometric_exchange([weighted], {"matrix": matrix}, self.context)
        self.assertEqual(first.data, second.data)
        self.assertEqual(first.facts["cardinality"], 4)
        self.assertEqual(first.data["weights"], [1, 2])
        self.assertEqual(verify_artifact(second, self.context).status, "exact")

    def test_large_integer_graph_uses_exact_fallback(self):
        scale = 10**10
        source = integer_code([], {"coordinates": [[scale*x for x in row] for row in self.source.data["coordinates"]],
                                  "norm_squared": 25*scale**2}, self.context)
        result = selection.isometric_exchange([source], {"normals": [[1, 2]]}, self.context)
        self.assertEqual(result.facts["cardinality"], 4)
        self.assertEqual(result.verification.details["isometric_exchange"]["integer_backend"], "Python integers")

    def test_two_different_codes_and_norm_alignment_are_coupled(self):
        # The second independently valid source is the changed-alphabet image
        # of the first, at a different integer-coordinate scale.
        image = integer_code([], {"coordinates": [[-15, 20], [-7, -24], [24, -7]],
                                 "norm_squared": 625}, self.context)
        result = selection.coupled_union([self.source, image], {}, self.context)
        self.assertEqual(result.facts["cardinality"], 4)
        self.assertEqual(result.verification.details["coupled_union"]["source_norm_alignment_scales"], [5, 1])
        self.assertEqual(verify_artifact(result, self.context).status, "exact")
        irrational_ratio = integer_code([], {"coordinates": [[1, 1]], "norm_squared": 2}, self.context)
        with self.assertRaisesRegex(ValueError, "rational square"):
            selection.coupled_union([self.source, irrational_ratio], {}, self.context)
        wrong_metric = integer_code([], {"coordinates": [[1, 0]], "weights": [25, 2], "norm_squared": 25}, self.context)
        with self.assertRaisesRegex(ValueError, "same explicit diagonal"):
            selection.coupled_union([self.source, wrong_metric], {}, self.context)

    def test_invalid_maps_inputs_and_budgets_reject_before_dense_allocation(self):
        cases = [{}, {"normals": [[0, 0]]}, {"normals": [[1]]}, {"normals": [[True, 0]]},
                 {"normals": [[0.5, 1]]}, {"matrix": [[1, 1], [0, 1]]},
                 {"normals": [], "matrix": [[1, 0], [0, 1]]}, {"normals": [], "typo": 1}]
        for params in cases:
            with self.subTest(params=params), self.assertRaises(ValueError):
                selection.isometric_exchange([self.source], params, self.context)
        for budget in ({"max_dimension": 1}, {"max_pairs": 14}, {"max_materialized_points": 5}):
            with patch.object(selection.np, "asarray", side_effect=AssertionError("dense allocation")):
                with self.subTest(budget=budget), self.assertRaises(ValueError):
                    selection.isometric_exchange([self.source], {"normals": [[1, 2]]}, RunContext(Path("/tmp"), budget))
        with self.assertRaises(ValueError):
            selection.isometric_exchange([Artifact("gram_spherical_code", {"gram": [[1]]})], {"normals": []}, self.context)
        bad = Artifact("spherical_code", {"coordinates": [[5, 0], [4, 3]], "norm_squared": 25},
                       verification=Verification("exact", "forged"))
        with self.assertRaisesRegex(ValueError, "inequality"):
            selection.isometric_exchange([bad], {"normals": []}, self.context)

    def test_registered_program_can_compose_coupled_selection(self):
        from discovery.runner import output_artifact, run_program
        program = {"schema_version": 1, "nodes": [
            {"id": "base", "op": "geometry.integer_code", "params": self.source.data},
            {"id": "exchange", "op": "selection.isometric_exchange", "inputs": ["base"], "params": {"normals": [[1, 2]]}}
        ], "outputs": ["exchange"]}
        result = output_artifact(run_program(program, self.context))
        self.assertEqual(result.facts["cardinality"], 4)
        self.assertEqual(verify_artifact(result, self.context).status, "exact")


if __name__ == "__main__":
    unittest.main()
