"""Exact joint-realization, nonorthogonal maps, overlaps and replay controls."""
from copy import deepcopy
from pathlib import Path
import unittest
from unittest.mock import patch

from discovery.core import Artifact, RunContext, Verification
from discovery.certification import verify_artifact
from discovery.compatibility import within, between, certified_union
from discovery.engines.geometry import integer_code, gram_code
from discovery.runner import input_kinds_match, output_artifact, run_program


class GluingTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))
        self.antipodes = integer_code([], {"coordinates": [[1, 0], [-1, 0]], "norm_squared": 1}, self.context)
        self.singleton = gram_code([], {"gram": [[1]]}, self.context)

    def test_mixed_nonorthogonal_two_dimensional_union(self):
        result = certified_union(self.antipodes, self.singleton, [[[1, 2]], [[-1, 2]]], self.context)
        self.assertEqual(result.kind, "gram_spherical_code")
        self.assertEqual(result.verification.status, "exact")
        self.assertEqual(result.facts["cardinality"], 3)
        self.assertEqual(result.facts["dimension"], 2)
        self.assertEqual(result.facts["max_inner_product"], [1, 2])
        self.assertEqual(verify_artifact(result, self.context).status, "exact")
        self.assertEqual(between(self.antipodes, self.singleton, [[[1, 2]], [[-1, 2]]], self.context).status, "exact")

    def test_rational_head_changes_union_rank_and_negative_scale_supported(self):
        second = gram_code([], {"gram": [[1, -1], [-1, 1]]}, self.context)
        cross = [[[1, 2], [-1, 2]], [[-1, 2], [1, 2]]]
        result = certified_union(self.antipodes, second, cross, self.context,
                                 heads=[[[1, 2]], [[-1, 2]]])
        self.assertEqual(result.facts["dimension"], 3)
        self.assertEqual(result.facts["norm_squared"], [5, 4])
        self.assertEqual(result.facts["cardinality"], 4)
        self.assertEqual(result.verification.details["input_ambient_dimensions"], [2, 1])
        with self.assertRaises(ValueError):
            certified_union(self.antipodes, second, cross, RunContext(Path("/tmp"), {"max_dimension": 2}),
                            heads=[[[1, 2]], [[-1, 2]]])
        negative = certified_union(self.antipodes, self.singleton, [[1], [-1]], self.context, scales=[-1, 1])
        self.assertEqual(negative.facts["cardinality"], 2)
        self.assertEqual(negative.facts["duplicates_removed"], 1)

    def test_rational_scales_match_different_source_norms(self):
        source = integer_code([], {"coordinates": [[2, 0], [-2, 0]], "norm_squared": 4}, self.context)
        result = certified_union(source, self.singleton, [[1], [-1]], self.context, scales=[[1, 2], 1])
        self.assertEqual(result.facts["norm_squared"], 1)
        self.assertEqual(result.facts["cardinality"], 3)
        self.assertEqual(result.data["gram"][0][-1], [1, 2])

    def test_partial_and_complete_overlaps_have_exact_cardinality(self):
        axes = integer_code([], {"coordinates": [[1, 0], [0, 1]], "norm_squared": 1}, self.context)
        other = gram_code([], {"gram": [[1, 0], [0, 1]]}, self.context)
        result = certified_union(axes, other, [[0, -1], [1, 0]], self.context)
        self.assertEqual(result.facts["cardinality"], 3)
        self.assertEqual(result.facts["duplicates_removed"], 1)
        self.assertEqual(result.data["gluing"]["classes"], [[[0, 0]], [[0, 1], [1, 0]], [[1, 1]]])
        identical = certified_union(axes, other, [[1, 0], [0, 1]], self.context)
        self.assertEqual(identical.facts["cardinality"], 2)
        self.assertEqual(identical.facts["duplicates_removed"], 2)

    def test_inconsistent_cross_psd_cannot_hide_behind_bounds_or_overlap(self):
        opposite = gram_code([], {"gram": [[1, -1], [-1, 1]]}, self.context)
        # Each cross entry satisfies the angle bound, but two antipodal rows
        # cannot both have inner product +1/2 against the same vector.
        with self.assertRaisesRegex(ValueError, "semidefinite|indefinite"):
            certified_union(self.antipodes, opposite,
                            [[[1, 2], [1, 2]], [[1, 2], [1, 2]]], self.context)
        axes = integer_code([], {"coordinates": [[1, 0], [0, 1]], "norm_squared": 1}, self.context)
        # A fake duplicate at the first coordinate must have zero product
        # with the second coordinate. Deduplication may not hide that conflict.
        with self.assertRaisesRegex(ValueError, "semidefinite|indefinite"):
            certified_union(axes, self.singleton, [[1], [[1, 2]]], self.context)

    def test_wrong_norm_threshold_shapes_and_noninjectivity_reject(self):
        cases = [({"scales": [0, 1]}, "noninjective"),
                 ({"scales": [1, 2]}, "same norm"),
                 ({"heads": [[], [1]]}, "equal width")]
        for params, reason in cases:
            with self.subTest(params=params), self.assertRaisesRegex(ValueError, reason):
                certified_union(self.antipodes, self.singleton, [[0], [0]], self.context, **params)
        with self.assertRaisesRegex(ValueError, "inequality"):
            certified_union(self.antipodes, self.singleton, [[[3, 4]], [[-3, 4]]], self.context)
        for cross in ([[0]], [[0, 0], [0, 0]], [[0.5], [-0.5]]):
            with self.subTest(cross=cross), self.assertRaises(ValueError):
                certified_union(self.antipodes, self.singleton, cross, self.context)
        with self.assertRaises(ValueError):
            within(self.antipodes, self.context, scale=0)
        with self.assertRaises(ValueError):
            within(self.antipodes, self.context, head=0)

    def test_stored_receipts_provenance_and_mutations(self):
        self.antipodes.data["provenance"] = {"note": "retained historical annotation"}
        self.antipodes.verification = Verification("unverified", "ignored saved status")
        result = certified_union(self.antipodes, self.singleton, [[[1, 2]], [[-1, 2]]], self.context)
        self.assertEqual(result.data["gluing"]["sources"][0]["data"]["provenance"], self.antipodes.data["provenance"])
        result.verification = Verification("unverified", "no receipt authority")
        self.assertEqual(verify_artifact(result, self.context).status, "exact")
        for field in ("gram", "classes", "source", "claim"):
            bad = deepcopy(result)
            if field == "gram": bad.data["gram"][0][1] = 0
            elif field == "classes": bad.data["gluing"]["classes"] = []
            elif field == "source": bad.data["gluing"]["sources"][0]["facts"]["cardinality"] = 99
            else: bad.facts["new_sota"] = True
            with self.subTest(field=field), self.assertRaises(ValueError):
                verify_artifact(bad, self.context)

    def test_nested_gluing_and_backend_resource_guards(self):
        glued = certified_union(self.antipodes, self.singleton, [[[1, 2]], [[-1, 2]]], self.context)
        repeated = certified_union(glued, self.singleton, [[[1, 2]], [[-1, 2]], [1]], self.context)
        self.assertEqual(repeated.facts["cardinality"], 3)
        self.assertEqual(verify_artifact(repeated, self.context).status, "exact")
        with self.assertRaises(ValueError):
            verify_artifact(repeated, RunContext(Path("/tmp"), {"max_nodes": 2}))
        for budget in ({"max_gram_entries": 8}, {"max_pairs": 2}, {"max_materialized_points": 2}):
            with patch("discovery.compatibility._source_gram", side_effect=AssertionError("source work started")):
                with self.subTest(budget=budget), self.assertRaises(ValueError):
                    certified_union(self.antipodes, self.singleton, [[0], [0]], RunContext(Path("/tmp"), budget))
        for source in (Artifact("implicit_spherical_code", {}), Artifact("interval_spherical_code", {}),
                       Artifact("spherical_code", {"real_coordinates": [[1.0]]})):
            with self.subTest(kind=source.kind), self.assertRaisesRegex(ValueError, "Unsupported|explicit"):
                within(source, self.context)

    def test_runner_accepts_only_explicit_declared_kind_alternatives(self):
        self.assertTrue(input_kinds_match([self.antipodes, self.singleton],
                        (("spherical_code", "gram_spherical_code"), ("spherical_code", "gram_spherical_code"))))
        self.assertFalse(input_kinds_match([self.singleton], ("spherical_code",)))
        program = {"schema_version": 1, "nodes": [
            {"id": "a", "op": "geometry.integer_code", "params": {"coordinates": [[1, 0], [-1, 0]], "norm_squared": 1}},
            {"id": "b", "op": "geometry.gram_code", "params": {"gram": [[1]]}},
            {"id": "union", "op": "gluing.union", "inputs": ["a", "b"],
             "params": {"cross_gram": [[[1, 2]], [[-1, 2]]]}}], "outputs": ["union"]}
        result = output_artifact(run_program(program, self.context))
        self.assertEqual((result.facts["cardinality"], result.facts["dimension"]), (3, 2))
        program["nodes"][-1]["params"] = {}
        with self.assertRaisesRegex(ValueError, "explicit.*cross_gram"):
            run_program(program, self.context)


if __name__ == "__main__":
    unittest.main()
