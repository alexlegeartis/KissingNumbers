"""Controls for bounded group discovery; numerical screening is never exact."""
from copy import deepcopy
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from discovery.core import Artifact, RunContext, Verification
from discovery.engines import groups


class GroupTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(workdir=Path("/tmp"), budget={})

    def test_weil_shape_unitarity_and_parameter_validation(self):
        representation = groups.weil([], {}, self.context)
        self.assertEqual(representation.facts["complex_dimension"], 5)
        self.assertEqual(representation.verification.status, "numerical")
        self.assertEqual(groups._generators(representation).shape, (5, 5, 5))
        for params in ({"p": 9}, {"parity": 0}, {"p": 5, "rank": 5}, {"typo": 1}):
            with self.assertRaises(ValueError):
                groups.weil([], params, self.context)

    def test_supplied_generators_support_other_dimensions(self):
        generators = [[[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                      [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]]
        representation = groups.from_generators([], {"generators_real": generators, "label": "signed permutations"}, self.context)
        self.assertEqual(representation.facts["complex_dimension"], 3)
        self.assertEqual(representation.verification.status, "numerical")
        self.assertFalse(representation.verification.details["finite_group_proved"])
        vectors, status = groups.orbit([1, 0, 0], groups._generators(representation), 6)
        self.assertEqual((len(vectors), status), (6, "closed_numerically"))
        complex_rep = groups.from_generators([], {"generators_real": [[[0]]], "generators_imag": [[[1]]]}, self.context)
        self.assertEqual(complex_rep.verification.status, "numerical")

    def test_supplied_generator_rejections_precede_allocation(self):
        oversized = {"generators_real": [[[1, 0], [0, 1]]]}
        for budget in ({"max_dimension": 1}, {"max_materialized_points": 7}):
            with patch.object(groups.np, "asarray", side_effect=AssertionError("allocated")):
                with self.assertRaises(ValueError):
                    groups.from_generators([], oversized, RunContext(Path("/tmp"), budget))
        for params in ({"generators_real": [[[2]]]}, {"generators_real": [[[True]]]},
                       {"generators_real": [[[1, 0]]]}, {"generators_real": [[[1]]], "typo": 1},
                       {"generators_real": [[[1]]], "generators_imag": []}):
            with self.assertRaises(ValueError):
                groups.from_generators([], params, self.context)

    def test_numerical_coordinate_guards_precede_allocation(self):
        cases=[([[1.,0.,0.]], {'max_dimension':2}),
               ([[1.,0.],[0.]], {}), ([[True,0.]], {}),
               ([[float('inf'),0.]], {})]
        for rows,budget in cases:
            a=Artifact('spherical_code',{'real_coordinates':rows})
            with self.subTest(rows=rows),patch.object(groups.np,'asarray',side_effect=AssertionError('allocated')):
                report=groups.verify_artifact(a,RunContext(Path('/tmp'),budget))
                self.assertEqual(report.status,'unverified')
    def test_orbit_limit_is_explicit_and_never_closed(self):
        generator = np.array([[1j]])
        vectors, status = groups.orbit([1], [generator], limit=3)
        self.assertEqual((len(vectors), status), (3, "limit_reached"))
        vectors, status = groups.orbit([1], [generator], limit=4)
        self.assertEqual((len(vectors), status), (4, "closed_numerically"))

    def test_compressed_checks_match_full_for_closed_vector_orbits(self):
        # Nontrivial coordinate mixing and phases exercise representative gauge.
        generators = [np.array([[0, 1], [-1, 0]], complex), np.diag([1j, -1j])]
        first, status = groups.orbit([1, 0], generators, 100)
        self.assertEqual(status, "closed_numerically")
        second, status = groups.orbit([1 / np.sqrt(2), 1j / np.sqrt(2)], generators, 100)
        self.assertEqual(status, "closed_numerically")
        for offset in (0, 0.13, 0.31):
            a, reps_a = groups.phase_block(first, 3, offset)
            b, _ = groups.phase_block(second, 2, 0.17)
            self.assertAlmostEqual(groups.compressed_max_cosine(reps_a, b), groups.max_cosine(a, b), places=10)
            self.assertAlmostEqual(groups.compressed_max_cosine(reps_a, a, True), groups.max_cosine(a), places=10)

    def test_phase_block_duplicates_are_not_charged_against_budget(self):
        full = np.array([[1 + 0j], [-1 + 0j]])
        one, _ = groups.phase_block(full, 1)
        two, _ = groups.phase_block(full, 2)
        self.assertEqual(groups._block_fingerprint(one), groups._block_fingerprint(two))
        self.assertEqual(groups._block_fingerprint(one), groups._block_fingerprint(one[::-1]))
        representation = groups.from_generators([], {"generators_real": [[[-1]]]}, self.context)
        result = groups.search([representation], {"trials": 0, "phase_orders": [1, 2, 3, 6],
                               "phase_steps": 1, "max_materialized_points": 8}, self.context)
        self.assertEqual(result.facts["cardinality"], 6)
        self.assertEqual(result.verification.status, "numerical")
        self.assertEqual(result.data["block_coverage"]["unique_blocks_admitted"], 2)
        self.assertEqual(result.data["block_coverage"]["admitted_point_rows"], 8)
        self.assertEqual(result.data["block_coverage"]["duplicate_blocks_skipped"], 2)
        self.assertFalse(result.data["materialization_limit_reached"])

    def test_numerical_verifier_never_promotes_and_rejects_forgery(self):
        artifact = Artifact(kind="spherical_code", data={"real_coordinates": [[1., 0.], [-1., 0.]]},
                            facts={"cardinality": 2},
                            verification=Verification(status="exact", method="forged", details={}))
        self.assertEqual(groups.verify_artifact(artifact).status, "numerical")
        artifact.facts["cardinality"] = 3
        self.assertEqual(groups.verify_artifact(artifact).status, "unverified")

    def test_search_is_deterministic_bounded_and_numerical(self):
        representation = groups.weil([], {"rank": 1}, self.context)
        params = {"trials": 2, "linecap": 20, "vectorcap": 80, "phase_steps": 2,
                  "phase_orders": [1, 2], "max_blocks": 20, "solver_seconds": 2}
        first = groups.search([representation], params, self.context)
        second = groups.search([representation], params, self.context)
        self.assertEqual(first.facts, second.facts)
        self.assertEqual(first.data["selected_blocks"], second.data["selected_blocks"])
        self.assertEqual(first.kind, "spherical_code")
        self.assertEqual(first.verification.status, "numerical")
        self.assertFalse(first.facts["global_optimality_claim"])
        params["vectorcap"] = 1
        limited = groups.search([representation], params, self.context)
        self.assertEqual(limited.kind, "search_result")
        self.assertEqual(limited.verification.status, "unverified")
        self.assertTrue(all(row["vector_status"] == "limit_reached" for row in limited.data["inventory"]))
        with self.assertRaises(ValueError):
            groups.search([representation], {"trial": 2}, self.context)

    def test_exact_control_reconstruction_and_mutation(self):
        artifact = groups.ganzhinov10([], {}, self.context)
        self.assertEqual(artifact.verification.status, "exact")
        self.assertEqual(artifact.verification.details["cardinality"], 510)
        self.assertEqual(artifact.verification.details["dimension"], 10)
        self.assertEqual(artifact.verification.details["pairs_checked"], 129795)
        bad = deepcopy(artifact)
        bad.data["coordinates"][0] = bad.data["coordinates"][1][:]
        self.assertEqual(groups.verify_artifact(bad).status, "unverified")
        bad = deepcopy(artifact)
        bad.facts["dimension"] = 9
        self.assertEqual(groups.verify_artifact(bad).status, "unverified")
        self.assertEqual(groups.verify_artifact(artifact, RunContext(workdir=Path("/tmp"), budget={"max_pairs": 10})).status, "unverified")

    def test_exact_rejects_pair_violation_and_false_facts(self):
        bad = Artifact(kind="spherical_code", data={"coordinates": [[2, 1], [2, -1]],
                       "weights": [1, 1], "norm_squared": 5})
        self.assertEqual(groups.verify_artifact(bad).status, "unverified")
        good = Artifact(kind="spherical_code", data={"coordinates": [[1, 0], [-1, 0]],
                        "weights": [1, 1], "norm_squared": 1}, facts={"max_inner_product": 0})
        self.assertEqual(groups.verify_artifact(good).status, "unverified")
        good.facts["max_inner_product"] = -1
        self.assertEqual(groups.verify_artifact(good).status, "exact")

    def test_exact_14_recipe_control(self):
        artifact = groups.ganzhinov14([], {}, self.context)
        self.assertEqual(artifact.verification.status, "exact")
        self.assertEqual(artifact.verification.details["cardinality"], 1932)
        self.assertEqual(artifact.verification.details["dimension"], 14)
        self.assertEqual(artifact.verification.details["pairs_checked"], 1865346)
        self.assertIn("recipe replay", artifact.data["provenance"])
        artifact.to_dict()

    def test_search_honors_verification_pair_budget(self):
        representation = groups.weil([], {"rank": 1}, self.context)
        context = RunContext(workdir=Path("/tmp"), budget={"max_pairs": 1, "max_trials": 1})
        result = groups.search([representation], {"trials": 10, "phase_steps": 1, "phase_orders": [1]}, context)
        self.assertEqual(result.kind, "search_result")
        self.assertEqual(result.data["output_status"], "verification_pair_limit_reached")
        self.assertEqual(result.data["limits"]["trials"], 1)
        self.assertEqual(result.verification.status, "unverified")


if __name__ == "__main__":
    unittest.main()
