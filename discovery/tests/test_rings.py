"""Exact finite-ring regression checks; no code coordinates are test inputs."""
from copy import deepcopy
from itertools import product
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from discovery.core import Artifact, RunContext, Verification
from discovery.engines import rings


def ring_recipe(n=7, q=4, index=0, context=None):
    context = context or RunContext(Path("/tmp"))
    pair = rings.binary_factor([], {"n": n, "index": index}, context)
    pair = rings.hensel([pair], {"modulus": q}, context)
    return rings.extend([pair], {"extension": -1}, context)


class RingTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def test_gray_isometry_all_symbol_pairs(self):
        for q in (4, 8, 16):
            table, weight = rings.gray_table(q)
            self.assertEqual(table.shape, (q, q//2))
            self.assertEqual(len({tuple(x) for x in table}), q)
            for a, b in product(range(q), repeat=2):
                self.assertEqual(int(np.count_nonzero(table[a] != table[b])), int(weight[(a-b) % q]))

    def test_octacode_exhaustive_distance_and_replay(self):
        ring = ring_recipe()
        code = rings.gray([ring], {}, self.context)
        self.assertEqual(code.facts, dict(length=16, cardinality=256, min_distance=6))
        # Independent full binary pair check, rather than only weights at zero.
        G = np.asarray(ring.data["generator"], dtype=np.int64)
        coefficients = np.array(list(product(range(4), repeat=4)))
        table, _ = rings.gray_table(4)
        words = table[coefficients @ G % 4].reshape(256, 16)
        distance = min(int(np.count_nonzero(words[i] != words[j]))
                       for i in range(len(words)) for j in range(i))
        self.assertEqual(distance, 6)
        code.verification = Verification()  # Saved receipt carries no authority.
        replay = rings.verify_artifact(code, self.context)
        self.assertEqual(replay.status, "exact")
        self.assertEqual(replay.details["facts"], code.facts)
        for artifact in (ring,):
            self.assertEqual(rings.verify_artifact(artifact, self.context).status, "exact")

    def test_modified_generator_distance_and_proof_rejected(self):
        code = rings.gray([ring_recipe()], {}, self.context)
        for field in ("generator", "distance", "histogram", "table", "proof_method"):
            changed = deepcopy(code)
            if field == "generator":
                changed.data["ring_code"]["generator"][0][0] ^= 1
            elif field == "distance":
                changed.facts["min_distance"] = 7
            elif field == "histogram":
                changed.data["proof"]["histogram"][6] += 1
            elif field == "table":
                changed.data["gray_table"][0][0] ^= 1
            else:
                changed.data["proof"]["method"] = "arbitrary_module.some_proof"
            with self.subTest(field=field), self.assertRaises(ValueError):
                rings.verify_artifact(changed, self.context)

    def test_budget_precedes_exhaustive_allocation(self):
        context = RunContext(Path("/tmp"), {"max_codewords": 255})
        with patch.object(rings.np, "arange", side_effect=AssertionError("must not enumerate")):
            candidate = rings.gray([ring_recipe()], {}, context)
        self.assertEqual(candidate.verification.status, "unverified")
        self.assertIsNone(candidate.facts["min_distance"])
        self.assertEqual(rings.verify_artifact(candidate, context).status, "unverified")
        with self.assertRaises(ValueError):
            rings.translate_union([candidate], {"target_distance": 4}, context)
        candidate.verification = Verification(status="exact", method="forged receipt")
        candidate.facts["min_distance"] = 6
        with self.assertRaises(ValueError):
            rings.translate_union([candidate], {"target_distance": 4}, context)

    def test_huge_parent_is_implicit_and_uncertified(self):
        with patch.object(rings.np, "arange", side_effect=AssertionError("must not enumerate")):
            candidate = rings.gray([ring_recipe(23, 8)], {}, self.context)
        self.assertEqual(candidate.facts["cardinality"], 2**36)
        self.assertEqual(candidate.facts["length"], 96)
        self.assertIsNone(candidate.facts["min_distance"])
        self.assertEqual(candidate.verification.status, "unverified")
        self.assertNotIn("histogram", candidate.data["proof"])

    def test_parameter_and_resource_limits(self):
        with self.assertRaises(ValueError):
            rings.binary_factor([], {"n": 10**9}, self.context)
        with self.assertRaises(ValueError):
            rings.binary_factor([], {"n": 8}, self.context)
        pair = rings.binary_factor([], {"n": 7}, self.context)
        for q in (3, 32, True):
            with self.subTest(q=q), self.assertRaises(ValueError):
                rings.hensel([pair], {"modulus": q}, self.context)
        with self.assertRaises(ValueError):
            rings.gray_table(32)
        with self.assertRaises(ValueError):
            rings.binary_factor([], {"nn": 7}, self.context)
        with self.assertRaises(ValueError):
            rings.gray([ring_recipe()], {"materialize": True}, self.context)
        ring = ring_recipe()
        binary = rings.gray([ring], {}, self.context)
        for name, operator in rings.OPERATORS.items():
            kind = rings.INPUT_KINDS[name]
            inputs = [] if not kind else [{"polynomial_pair": pair, "ring_code": ring,
                                          "binary_code": binary}[kind[0]]]
            with self.subTest(operator=name), self.assertRaises(ValueError):
                operator(inputs, {"typo": True}, self.context)
        limited = RunContext(Path("/tmp"), {"max_dgls_words": 2**24-1})
        with patch.object(rings, "_load_parent", side_effect=AssertionError("must not start proof")):
            with self.assertRaises(ValueError):
                rings.verify_dgls([ring_recipe(23, 8)], {}, limited)

    def test_union_replays_masks_and_conservative_multiplicity(self):
        parent = rings.gray([ring_recipe(7, 8)], {}, self.context)
        union = rings.translate_union([parent], {"target_distance": 8}, self.context)
        self.assertEqual(union.data["proof"]["multiplier"], 2)
        self.assertEqual(union.facts["cardinality"], parent.facts["cardinality"]*2)
        self.assertGreaterEqual(union.facts["min_distance"], 8)
        self.assertEqual(rings.verify_artifact(union, self.context).status, "exact")
        changed = deepcopy(union)
        changed.data["proof"]["selected_repeated_masks"][1] = [0, 0, 0, 0]
        with self.assertRaises(ValueError):
            rings.verify_artifact(changed, self.context)

    def test_dgls_synthesized_two_proofs_and_implicit_union(self):
        # This dedicated compressed proof enumerates 2^24 Z4 words internally,
        # never materializes the 2^36 Gray parent or 2^37 union.
        code = rings.verify_dgls([ring_recipe(23, 8)], {}, self.context)
        self.assertEqual(code.facts, dict(length=96, cardinality=2**36, min_distance=24))
        self.assertEqual(code.data["proof"]["first"], code.data["proof"]["second"])
        self.assertEqual(code.data["proof"]["first"]["z4_words_enumerated"], 2**24)
        # Public union rechecks both exact parent proofs, even if its receipt is erased.
        code.verification = Verification()
        union = rings.translate_union([code], {"target_distance": 24}, self.context)
        self.assertEqual(union.facts, dict(length=96, cardinality=2**37, min_distance=24))
        self.assertEqual(union.data["proof"]["multiplier"], 2)
        self.assertNotIn("words", union.data)


if __name__ == "__main__":
    unittest.main()
