"""Small independent controls for the competitive research exchange machinery."""
import importlib.util
from pathlib import Path
import time
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location("coupled28", Path(__file__).with_name("coupled_search.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
proof_spec = importlib.util.spec_from_file_location("blocker_proof", Path(__file__).with_name("verify_pool_obstruction.py"))
proof_module = importlib.util.module_from_spec(proof_spec)
proof_spec.loader.exec_module(proof_module)


class CoupledSearchTests(unittest.TestCase):
    def test_full_sextet_family_count_and_partitions(self):
        frames, outside = module.all_sextets()
        self.assertEqual((len(frames), len(outside)), (1771, 1736))
        self.assertEqual(len(set(frames)), 1771)
        for frame in frames:
            self.assertEqual(sorted(sum((list(row) for row in frame), [])), list(range(24)))

    def test_fractional_image_is_retained_by_exact_common_scaling(self):
        frame = tuple(tuple(range(i, i + 4)) for i in range(0, 24, 4))
        point = np.zeros((1, 28), np.int16)
        point[0, :5] = [11, 2, 1, 1, 1]
        self.assertEqual(int(np.sum(point.astype(np.int64)**2)), 128)
        original, transformed, norm, factor, matrix, denominator = module.image(point, [frame], [0])
        self.assertEqual((factor, norm, denominator), (2, 512, 2))
        self.assertTrue(np.any(transformed % 2))
        self.assertTrue(np.array_equal(transformed, point.astype(np.int64) @ matrix.T))
        self.assertTrue(np.array_equal(original, point * 2))

    def test_positive_exchange_and_independent_integer_replay(self):
        original = np.array([[8, 8, 0], [-8, 0, 8], [0, -8, -8]], np.int16)
        transformed = -original
        result, report, witness = module.exchange(original, transformed, 128, blocker_limit=8,
                                                  max_checks=100, max_edges=100, deadline=time.monotonic() + 10)
        self.assertEqual(report["state"], "complete")
        self.assertEqual((report["gain"], len(result)), (3, 6))
        added = transformed[witness["added_image_indices"]]
        direct = module.direct_delta(original, result, added)
        self.assertEqual(direct["status"], "PASS")
        self.assertEqual(direct["direct_integer_pairs"], 15)
        gram = result.astype(np.int64) @ result.T.astype(np.int64)
        np.fill_diagonal(gram, -128)
        self.assertLessEqual(int(gram.max()), 64)

    def test_work_cutoff_is_incomplete_not_zero_gain(self):
        original = np.array([[8, 8, 0], [-8, 0, 8], [0, -8, -8]], np.int16)
        result, report, witness = module.exchange(original, -original, 128, blocker_limit=8,
                                                  max_checks=1, max_edges=100, deadline=time.monotonic() + 10)
        self.assertEqual(report["state"], "incomplete")
        self.assertNotIn("gain", report)
        self.assertIsNone(witness)
        self.assertTrue(np.array_equal(result, original))

    def test_edge_cutoff_is_incomplete_before_accumulation(self):
        source = np.array([[8, 8, 0, 0, 0], [-8, -8, 0, 0, 0]], np.int16)
        candidates = np.array([[10, 2, 4, 2, 2], [-10, -2, -4, -2, -2]], np.int16)
        _, report, witness = module.exchange(source, candidates, 128, blocker_limit=8,
                                             max_checks=100, max_edges=1, deadline=time.monotonic() + 10)
        self.assertEqual(report["state"], "incomplete")
        self.assertEqual(report["reason"], "edge_cutoff")
        self.assertIsNone(witness)

    def test_two_images_share_one_deletion_for_positive_gain(self):
        baseline = np.array([[8, 8, 0, 0, 0]], np.int16)
        first = np.array([[10, 2, 4, 2, 2]], np.int16)
        second = np.array([[2, 10, -4, -2, -2]], np.int16)
        self.assertEqual(int(first[0].astype(np.int64) @ second[0].astype(np.int64)), 16)
        pools = [{"points": points, "indptr": np.array([0, 1]), "blockers": np.array([0]), "factor": 1}
                 for points in (first, second)]
        result = module.shared_blocker_search(baseline, pools)
        self.assertEqual((result["state"], result["gain"]), ("gain", 1))
        self.assertEqual(result["proof"]["status"], "PASS")
        self.assertEqual(result["proof"]["cardinality"], 2)

    def test_certified_core_skipping_matches_full_conflict_graph(self):
        small = np.array([[8, 8, 0], [-8, 0, 8], [0, -8, -8]], np.int16)
        original = np.pad(small, ((0, 0), (0, 25)))
        transformed = -original
        lookup = {row[:24].tobytes() for row in np.concatenate((original, transformed))}
        result, report, witness = module.exchange(original, transformed, 128, blocker_limit=8,
                                                  max_checks=1, max_edges=100, deadline=time.monotonic() + 10,
                                                  core_lookup=lookup)
        self.assertEqual(report["state"], "complete")
        self.assertEqual(report["gain"], 3)
        self.assertEqual(report["dot_products"], 0)
        self.assertEqual(report["certified_core_pairs_skipped"], 9)
        self.assertEqual(module.direct_delta(original, result, transformed[witness["added_image_indices"]])["status"], "PASS")

    def test_global_blocker_injection_and_negative_controls(self):
        source = np.array([[8, 8, 0, 0, 0], [-8, -8, 0, 0, 0]], np.int64)
        candidates = np.array([[10, 2, 4, 2, 2], [-10, -2, -4, -2, -2]], np.int64)
        report = proof_module.check(source, candidates, np.array([0, 1]), 128)
        self.assertEqual(report["maximum_possible_gain"], 0)
        for matching in (np.array([0, 0]), np.array([1, 0])):
            with self.assertRaises(ValueError):
                proof_module.check(source, candidates, matching, 128)
        bad = candidates.copy()
        bad[0, -1] = -2**63
        with self.assertRaisesRegex(ValueError, "magnitude"):
            proof_module.check(source, bad, np.array([0, 1]), 128)


if __name__ == "__main__":
    unittest.main()
