"""Exact syndrome quotient checks against independent small binary enumerations."""
from copy import deepcopy
from itertools import product
from pathlib import Path
import random
import unittest
from unittest.mock import patch

from discovery.core import Artifact, RunContext, Verification
from discovery.engines import cosets


def extended_hamming_check():
    # Columns are (1,x) for the eight x in F2^3: an [8,4,4] kernel.
    return [[1]*8] + [[(j >> i) & 1 for j in range(8)] for i in range(3)]


def brute_cosets(H):
    groups = {s: [] for s in range(1 << len(H))}
    for word in product((0, 1), repeat=len(H[0])):
        syndrome = sum((sum(a*b for a, b in zip(row, word)) % 2) << i
                       for i, row in enumerate(H))
        groups[syndrome].append(word)
    return groups


class CosetTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def build(self, H=None, context=None):
        return cosets.build([], {"parity_check": H or extended_hamming_check()},
                            context or self.context)

    def test_extended_hamming_kernel_exact_distance_and_replay(self):
        oracle = self.build()
        self.assertEqual(oracle.kind, "syndrome_oracle")
        self.assertEqual(oracle.facts, dict(length=8, rank=4, kernel_dimension=4,
                                           cardinality=16, min_distance=4))
        self.assertEqual(oracle.data["proof"]["syndrome_min_weights"][0], 0)
        self.assertEqual(oracle.data["proof"]["kernel_min_nonzero_weight"], 4)
        self.assertEqual(len(oracle.data["proof"]["syndrome_min_weights"]), 16)
        oracle.verification = Verification()  # Serialized status is never proof.
        self.assertEqual(cosets.verify_artifact(oracle).status, "exact")

    def test_all_syndrome_distances_match_independent_enumeration(self):
        H = extended_hamming_check()
        oracle = self.build(H)
        groups = brute_cosets(H)
        weights = [min(sum(w) for w in groups[s]) for s in range(16)]
        self.assertEqual(oracle.data["proof"]["syndrome_min_weights"], weights)
        for left in range(16):
            for right in range(left):
                observed = min(sum(x != y for x, y in zip(a, b))
                               for a in groups[left] for b in groups[right])
                self.assertEqual(observed, weights[left ^ right])

    def test_fixed_coset_unions_match_full_binary_pair_distance(self):
        H = extended_hamming_check()
        oracle = self.build(H)
        groups = brute_cosets(H)
        for labels, expected in (([0], 4), ([0, 2], 2), ([0, 1], 1), ([3, 5, 9], 2)):
            with self.subTest(labels=labels):
                code = cosets.union([oracle], {"syndrome_labels": labels}, self.context)
                words = [word for label in labels for word in groups[label]]
                distance = min(sum(x != y for x, y in zip(words[i], words[j]))
                               for i in range(len(words)) for j in range(i))
                self.assertEqual(distance, expected)
                self.assertEqual(code.facts, dict(length=8, cardinality=16*len(labels),
                                                 min_distance=expected))
                self.assertEqual(cosets.verify_artifact(code).details["facts"], code.facts)
                self.assertNotIn("words", code.data)

    def test_nonempty_dependency_dp_does_not_use_empty_zero_syndrome(self):
        # No weight-one kernel word; duplicate nonzero columns give weight two.
        oracle = self.build([[1, 1]])
        self.assertEqual(oracle.data["proof"]["syndrome_min_weights"], [0, 1])
        self.assertEqual(oracle.facts["min_distance"], 2)
        # An actual zero column must, separately, certify distance one.
        oracle = self.build([[1, 0]])
        self.assertEqual(oracle.facts["min_distance"], 1)

    def test_trivial_kernel_sentinel_is_not_claimed_as_actual_distance(self):
        oracle = self.build([[1, 0], [0, 1]])
        self.assertEqual(oracle.facts["cardinality"], 1)
        self.assertIsNone(oracle.facts["min_distance"])
        self.assertEqual(oracle.data["proof"]["kernel_min_nonzero_weight"], 3)
        one = cosets.union([oracle], {"syndrome_labels": [2]}, self.context)
        self.assertIsNone(one.facts["min_distance"])
        two = cosets.union([oracle], {"syndrome_labels": [0, 3]}, self.context)
        self.assertEqual(two.facts, dict(length=2, cardinality=2, min_distance=2))
        self.assertEqual(cosets.verify_artifact(two).status, "exact")

    def test_random_full_rank_matrices_against_enumeration(self):
        rng = random.Random(21063)
        for n in range(2, 9):
            r = min(3, n)
            # Identity minor certifies rank; remaining columns are independent draws.
            H = [[int(i == j) if j < r else rng.randrange(2) for j in range(n)]
                 for i in range(r)]
            oracle = self.build(H)
            groups = brute_cosets(H)
            self.assertEqual(oracle.data["proof"]["syndrome_min_weights"],
                             [min(map(sum, groups[s])) for s in range(1 << r)])
            nonempty = [sum(w) for w in groups[0] if any(w)]
            self.assertEqual(oracle.data["proof"]["kernel_min_nonzero_weight"],
                             min(nonempty, default=n+1))

    def test_implicit_kernel_can_be_much_larger_than_quotient_budget(self):
        oracle = self.build([[1]*100], RunContext(Path("/tmp"), {"max_codewords": 2}))
        self.assertEqual(oracle.facts["cardinality"], 2**99)
        self.assertEqual(oracle.facts["min_distance"], 2)
        self.assertEqual(len(oracle.data["proof"]["syndrome_min_weights"]), 2)

    def test_matrix_rejections_and_rank_checked_before_tables(self):
        invalid = ([], [[]], [[1, 0], [1]], [[1, 2]], [[True]], [[1.0]],
                   [[0]], [[1, 0], [1, 0]], [[1], [1]])
        for H in invalid:
            with self.subTest(H=H), patch.object(cosets, "array", side_effect=AssertionError("allocated")):
                with self.assertRaises(ValueError):
                    cosets.build([], {"parity_check": H}, self.context)

    def test_caps_precede_exponential_storage(self):
        oversized = [[int(i == j) for j in range(17)] for i in range(17)]
        limited = RunContext(Path("/tmp"), {"max_codewords": 15})
        dimensions = RunContext(Path("/tmp"), {"max_dimension": 7})
        for H, context in ((oversized, self.context), (extended_hamming_check(), limited),
                           (extended_hamming_check(), dimensions)):
            with self.subTest(rows=len(H), budget=context.budget), patch.object(
                    cosets, "array", side_effect=AssertionError("must not allocate")):
                with self.assertRaises(ValueError):
                    self.build(H, context)

    def test_modified_certificates_facts_and_receipts(self):
        oracle = self.build()
        code = cosets.union([oracle], {"syndrome_labels": [0, 2]}, self.context)
        changed = []
        for artifact in (oracle, code):
            bad = deepcopy(artifact)
            bad.facts["min_distance"] = 99
            bad.verification = Verification("exact", "forged")
            changed.append(bad)
        for key, value in (("kernel_min_nonzero_weight", 8),
                           ("syndrome_min_weights", [0]*16), ("method", "forged proof")):
            bad = deepcopy(oracle)
            bad.data["proof"][key] = value
            changed.append(bad)
        bad = deepcopy(oracle)
        bad.data["parity_check"][1] = bad.data["parity_check"][0][:]
        changed.append(bad)
        bad = deepcopy(code)
        bad.data["syndrome_labels"] = [0, 0]
        changed.append(bad)
        bad = deepcopy(code)
        bad.data["proof"]["cross_min_distance"] = 4
        changed.append(bad)
        bad = deepcopy(oracle)
        bad.facts["rank"] = 4.0
        changed.append(bad)
        for artifact in changed:
            with self.subTest(artifact=artifact), self.assertRaises(ValueError):
                cosets.verify_artifact(artifact)

    def test_labels_pair_budget_and_parameter_validation(self):
        oracle = self.build()
        for labels in ([], [0, 0], [-1], [16], [True], [1.0], "01"):
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                cosets.union([oracle], {"syndrome_labels": labels}, self.context)
        limited = RunContext(Path("/tmp"), {"max_pairs": 1})
        with self.assertRaises(ValueError):
            cosets.union([oracle], {"syndrome_labels": [0, 2, 4]}, limited)
        for operator, inputs in ((cosets.build, []), (cosets.union, [oracle])):
            with self.assertRaises(ValueError):
                operator(inputs, {"unknown": True}, self.context)
            with self.assertRaises(ValueError):
                operator(inputs, {}, self.context)
        with self.assertRaises(ValueError):
            cosets.build([oracle], {"parity_check": [[1]]}, self.context)
        with self.assertRaises(ValueError):
            cosets.union([Artifact("binary_code", {})], {"syndrome_labels": [0]}, self.context)
        forged = deepcopy(oracle)
        forged.facts["cardinality"] *= 2
        with self.assertRaises(ValueError):
            cosets.union([forged], {"syndrome_labels": [0]}, self.context)

    def test_json_roundtrip_and_saved_receipt_has_no_authority(self):
        oracle = self.build()
        code = cosets.union([oracle], {"syndrome_labels": [0, 2]}, self.context)
        code = Artifact.from_dict(code.to_dict())
        code.verification = Verification("numerical", "irrelevant saved receipt")
        self.assertEqual(cosets.verify_artifact(code).status, "exact")


if __name__ == "__main__":
    unittest.main()
