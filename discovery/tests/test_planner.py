"""Synthesis coverage, durable interruption recovery, and discovery-quality tests."""
import contextlib
from copy import deepcopy
import io
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from discovery.__main__ import main
from discovery.certification import verify_artifact
from discovery.core import RunContext, canonical_json
from discovery.runner import output_artifact, registry, source_fingerprints
from discovery.planner import (compile_catalogue, extension_at, extension_counts,
                               operator_catalogue, plan_search, read_ledger,
                               export_receipt, _domain, _execute, _program, _worker)


def specification(depth=2):
    return {"schema_version": 1, "max_nodes": depth, "output_kinds": ["spherical_code"],
            "catalogue": [{"op": "geometry.integer_code", "params": {"coordinates": [[1], [-1]], "norm_squared": 1}},
                          {"op": "geometry.orthogonal_union"}]}


def coset_specification(items=(2, 4)):
    return {"schema_version": 1, "max_nodes": 3, "strategy": "adaptive",
            "output_kinds": ["implicit_spherical_code"],
            "quality_targets": [{"dimension": 8, "cardinality": 240, "label": "E8 known construction"}],
            "catalogue": [
                {"op": "cosets.build", "params": {"parity_check": [[1] * 8] + [[(j >> i) & 1 for j in range(8)] for i in range(3)]}},
                {"op": "cosets.union", "domains": {"syndrome_labels": {"subsets": {"of": list(items), "include": [0]}}}},
                {"op": "geometry.sign_embed"}],
            "unavailable": [{"name": "geometry.general_compatible_union", "family": "mixed",
                             "reason": "Needed to combine implicit sign points with nonorthogonal support layers"}]}


def rows(directory):
    with sqlite3.connect(Path(directory) / "planner.sqlite3") as db:
        return [(i, json.loads(nodes), state) for i, nodes, state in db.execute("SELECT id,nodes,eval_state FROM programs WHERE id>0 ORDER BY id")]


def direct_execute(program, context, seconds, max_receipt_bytes):
    return _worker({"program": program, "workdir": str(context.workdir), "budget": context.budget,
                    "max_receipt_bytes": max_receipt_bytes})


class PlannerTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def test_all_registered_operators_visible_without_dimension_dispatch(self):
        catalogue = operator_catalogue()
        self.assertEqual(set(catalogue), set(registry()[0]))
        for op, contract in catalogue.items():
            with self.subTest(op=op):
                self.assertTrue(contract["implemented"])
                self.assertTrue(contract["output_kinds"])
                spec = specification(1)
                spec["catalogue"] = [{"op": op}]
                self.assertEqual(compile_catalogue(spec, self.context)[0].op, op)

    def test_domains_match_independent_cartesian_and_combination_enumerators(self):
        for start, stop, step in [(1, 8, 2), (8, 1, -2), (-5, 5, 3)]:
            domain = _domain({"integer_range": {"start": start, "stop": stop, "step": step}})
            self.assertEqual([domain.at(i) for i in range(domain.count)], list(range(start, stop, step)))
        # Huge finite ranges are counted, never allocated.
        domain = _domain({"integer_range": {"start": 0, "stop": 10**50}})
        self.assertEqual((domain.count, domain.at(10**40)), (10**50, 10**40))
        for reverse in (False, True):
            domain = _domain({"subsets": {"of": [2, 4, 6, 8], "include": [0], "min_size": 1, "max_size": 3,
                                          "largest_first": reverse}})
            sizes = sorted(range(1, 4), reverse=reverse)
            expected = [[0] + list(c) for n in sizes for c in itertools.combinations([2, 4, 6, 8], n)]
            self.assertEqual([domain.at(i) for i in range(domain.count)], expected)
        spec = specification(1)
        spec["catalogue"][0]["domains"] = {"weights": {"values": [[1], [2]]}, "norm_squared": {"values": [1, 2]}}
        del spec["catalogue"][0]["params"]["norm_squared"]
        entry = compile_catalogue(spec, self.context)[0]
        self.assertEqual([(entry.params_at(i)["norm_squared"], entry.params_at(i)["weights"]) for i in range(4)],
                         [(1, [1]), (1, [2]), (2, [1]), (2, [2])])

    def test_synthesized_graphs_equal_independent_brute_force_tiny_language(self):
        spec = specification(3)
        spec["catalogue"].insert(1, {"op": "groups.from_generators", "params": {"generators_real": [[[-1]]]}})
        spec["catalogue"].append({"op": "rings.extend"})  # Unreachable input type must be visible in ledger.
        entries = compile_catalogue(spec, self.context)
        expected = []
        frontier = [([], [])]
        # Independent simple product enumeration: no production ordinal decoder.
        for depth in range(1, 4):
            following = []
            for nodes, kinds in frontier:
                for entry in entries:
                    choices = [[i for i, kinds_at_i in enumerate(kinds) if need in kinds_at_i] for need in entry.inputs]
                    for binding in itertools.product(*choices):
                        child = nodes + [{"id": f"n{len(nodes)}", "op": entry.op,
                                          "params": entry.params, "inputs": [f"n{i}" for i in binding]}]
                        output_kinds = kinds + [list(entry.outputs)]
                        if depth < 3 or "spherical_code" in entry.outputs:
                            expected.append(canonical_json(child))
                        following.append((child, output_kinds))
            frontier = following
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", return_value={"state": "failed", "reason": "deliberate failure"}):
            result = plan_search(spec, self.context, Path(temporary) / "run", max_steps=10000, max_evaluations=1000)
            actual = [canonical_json(nodes) for _, nodes, _ in rows(Path(temporary) / "run")]
            self.assertEqual(actual, expected)
            self.assertTrue(result["enumeration_complete"])
            self.assertEqual(result["status"], "incomplete")
            self.assertGreater(result["pruned_non_target_leaves"], 0)
            self.assertGreater(result["pruning_by_catalogue_entry"]["2"]["type_incompatible_extensions"], 0)
            self.assertGreater(result["pruning_by_catalogue_entry"]["3"]["prefixes_without_compatible_inputs"], 0)

    def test_resume_matches_uninterrupted_program_ids_and_receipts(self):
        for strategy in ("breadth_first", "adaptive"):
            spec = specification(3)
            spec["strategy"] = strategy
            with self.subTest(strategy=strategy), tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
                full = Path(temporary) / "full"
                split = Path(temporary) / "split"
                expected = plan_search(spec, self.context, full, max_evaluations=100)
                result = plan_search(spec, self.context, split, max_steps=2)
                self.assertFalse(result["enumeration_complete"])
                while not result["enumeration_complete"]:
                    result = plan_search(spec, self.context, split, resume=True, max_steps=3)
                self.assertEqual(rows(full), rows(split))
                self.assertEqual(expected["best_exact_by_dimension"], result["best_exact_by_dimension"])
                self.assertEqual(expected["pruning_by_catalogue_entry"], result["pruning_by_catalogue_entry"])
                self.assertEqual(result, read_ledger(split))

    def test_adaptive_interleaves_roots_and_evaluates_pending_children(self):
        spec = specification(2)
        spec["strategy"] = "adaptive"
        spec["catalogue"][0]["params"].pop("coordinates")
        spec["catalogue"][0]["domains"] = {"coordinates": {"values": [
            [[1], [-1]], [[1, 0], [-1, 0], [0, 1], [0, -1]]]}}
        seen = []
        def execute(program, *args):
            seen.append(program)
            return direct_execute(program, *args)
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=execute):
            plan_search(spec, self.context, Path(temporary) / "run")
        self.assertEqual([len(p["nodes"]) for p in seen[:4]], [1, 2, 1, 2])
        # A fresh prefix reaches a child before the next root; the persistent
        # queue still gives the second constructor and its child service.
        self.assertEqual(len(seen[3]["nodes"][0]["params"]["coordinates"]), 4)

    def test_interrupt_after_committed_start_requeues_without_losing_frontier(self):
        def interrupt(event):
            if event["state"] == "started":
                raise KeyboardInterrupt()
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            directory = Path(temporary) / "run"
            with self.assertRaises(KeyboardInterrupt):
                plan_search(specification(), self.context, directory, on_progress=interrupt)
            self.assertEqual(read_ledger(directory)["evaluation_attempts"], {"running": 1})
            resumed = plan_search(specification(), self.context, directory, resume=True)
            self.assertEqual(resumed["evaluation_attempts"], {"completed": 3, "interrupted": 1})
            self.assertTrue(resumed["enumeration_complete"])
            self.assertEqual(len(rows(directory)), 3)

    def test_interrupt_during_next_expansion_replays_exactly_that_edge(self):
        calls = 0
        def interrupt(*args):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise KeyboardInterrupt()
            return extension_at(*args)
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            directory = Path(temporary) / "run"
            with patch("discovery.planner.extension_at", side_effect=interrupt), self.assertRaises(KeyboardInterrupt):
                plan_search(specification(), self.context, directory)
            resumed = plan_search(specification(), self.context, directory, resume=True)
            self.assertEqual(resumed["generated_programs"], 3)
            self.assertEqual(resumed["completed_programs"], 3)

    def test_checkpoint_identity_and_single_owner_are_enforced(self):
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            directory = Path(temporary) / "run"
            def attempt_concurrent(event):
                if event["state"] == "started":
                    with self.assertRaisesRegex(ValueError, "in use"):
                        plan_search(specification(), self.context, directory, resume=True)
            plan_search(specification(), self.context, directory, on_progress=attempt_concurrent)
            changed = specification(3)
            with self.assertRaisesRegex(ValueError, "same specification"):
                plan_search(changed, self.context, directory, resume=True)
            with self.assertRaisesRegex(ValueError, "same specification"):
                plan_search(specification(), RunContext(Path("/tmp"), {"max_pairs": 1}), directory, resume=True)
            with patch("discovery.planner.source_fingerprints", return_value={"changed": "source"}), self.assertRaisesRegex(ValueError, "source fingerprints"):
                plan_search(specification(), self.context, directory, resume=True)
            with self.assertRaises(FileExistsError):
                plan_search(specification(), self.context, directory)

    def test_resource_stops_leave_unsearched_work_and_can_resume_with_larger_limits(self):
        cases = ({"max_programs": 1}, {"max_program_bytes": 1}, {"max_evaluations": 1}, {"time_limit_seconds": 1e-9})
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            for i, limits in enumerate(cases):
                directory = Path(temporary) / str(i)
                result = plan_search(specification(), self.context, directory, **limits)
                self.assertEqual(result["status"], "incomplete")
                self.assertIsNone(result["unsearched"]["total_remaining_programs"])
                self.assertGreater(result["unsearched"]["unexpanded_prefixes"] + result["unsearched"]["generated_pending_evaluations"], 0)
                result = plan_search(specification(), self.context, directory, resume=True)
                self.assertTrue(result["enumeration_complete"])
                self.assertEqual(result["completed_programs"], 3)

    def test_real_worker_resource_failures_and_timeouts_are_not_exclusions(self):
        program = _program([{"id": "n0", **specification()["catalogue"][0]}], "test")
        self.assertEqual(_execute(program, self.context, 1e-6, 100000)["state"], "timed_out")
        self.assertEqual(_execute(program, self.context, 10, 1)["state"], "resource_limited")
        with tempfile.TemporaryDirectory() as temporary:
            result = plan_search(specification(), RunContext(Path("/tmp"), {"max_pairs": 1}), Path(temporary) / "run")
            self.assertEqual(result["program_states"], {"completed": 2, "resource_limited": 1})
            self.assertEqual(result["status"], "incomplete")
            self.assertTrue(result["enumeration_complete"])
            self.assertFalse(result["mathematical_exclusion_claim"])

    def test_retry_unresolved_keeps_attempt_history_and_recovers_with_larger_receipt_limit(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "run"
            initial = plan_search(specification(), self.context, directory, max_receipt_bytes=1)
            self.assertEqual(initial["program_states"], {"resource_limited": 3})
            resumed = plan_search(specification(), self.context, directory, resume=True, retry_unresolved=True)
            self.assertEqual(resumed["evaluation_attempts"], {"completed": 3, "resource_limited": 3})
            self.assertEqual(resumed["status"], "complete_within_declared_scope")

    def test_quality_deficit_is_visible_even_if_enumeration_completes(self):
        spec = specification(1)
        spec["quality_targets"] = [{"dimension": 28, "cardinality": 221988, "label": "dated repo baseline"}]
        for status in ("numerical", "unverified"):
            fake = {"state": "completed", "receipt": {"source_sha256": source_fingerprints()}, "outputs": {"n0": {"kind": "spherical_code",
                    "verification_status": status, "candidate_facts": {"dimension": 28, "cardinality": 6720}, "facts": {}}}}
            with self.subTest(status=status), tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", return_value=fake):
                result = plan_search(spec, self.context, Path(temporary) / "run")
                self.assertTrue(result["enumeration_complete"])
                self.assertFalse(result["target_quality_passed"])
                self.assertEqual(result["best_exact_by_dimension"], {})
                quality = result["quality_diagnostics"][0]
                self.assertEqual(quality["candidate_deficit"], 215268)
                self.assertLess(quality["candidate_fraction_of_target"], 0.031)
                self.assertIsNone(quality["best_exact_cardinality"])

    def test_conditional_group_output_retains_unverified_ambient_quality_and_limits(self):
        spec = {"schema_version": 1, "max_nodes": 2, "output_kinds": ["spherical_code", "search_result"],
                "quality_targets": [{"dimension": 2, "cardinality": 6}], "catalogue": [
                    {"op": "groups.from_generators", "params": {"generators_real": [[[-1]]]}},
                    {"op": "groups.search", "params": {"trials": 0, "phase_orders": [1], "phase_steps": 1, "max_output_points": 1}}]}
        with tempfile.TemporaryDirectory() as temporary:
            result = plan_search(spec, self.context, Path(temporary) / "run")
            self.assertEqual(result["completed_programs"], 1)
            self.assertEqual(result["best_exact_by_dimension"], {})
            quality = result["quality_diagnostics"][0]
            self.assertEqual(quality["best_observed_candidate"]["verification_status"], "unverified")
            self.assertEqual(quality["candidate_deficit"], 4)
            self.assertEqual(result["operator_search_diagnostics"][0]["diagnostics"]["output_status"], "output_limit_reached")

    def test_actual_coset_search_finds_composed_candidates_and_matches_brute_force(self):
        spec = coset_specification()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "run"
            # Real disposable Python workers; graph edges and label subsets are synthesized.
            result = plan_search(spec, self.context, directory, max_steps=1000)
            self.assertEqual(result["status"], "complete_within_declared_scope")
            self.assertEqual(result["completed_programs"], 4)
            self.assertEqual(result["best_exact_by_dimension"]["8"]["cardinality"], 48)
            self.assertFalse(result["target_quality_passed"])
            with sqlite3.connect(directory / "planner.sqlite3") as db:
                completed = db.execute("SELECT result FROM programs WHERE eval_state='completed'").fetchall()
            H = spec["catalogue"][0]["params"]["parity_check"]
            for (encoded,) in completed:
                receipt = json.loads(encoded)["receipt"]
                artifact = output_artifact(receipt)
                labels = artifact.data["binary_code"]["data"]["syndrome_labels"]
                words = []
                for word in itertools.product((0, 1), repeat=8):
                    syndrome = sum((sum(a*b for a, b in zip(row, word)) % 2) << i for i, row in enumerate(H))
                    if syndrome in labels:
                        words.append(word)
                distance = min(sum(a != b for a, b in zip(left, right)) for i, left in enumerate(words) for right in words[:i])
                binary_facts = artifact.data["binary_code"]["facts"]
                self.assertEqual((binary_facts["cardinality"], binary_facts["min_distance"]), (len(words), distance))
                self.assertEqual(verify_artifact(artifact, self.context).status, "exact")
                self.assertEqual([node["op"] for node in receipt["program"]["nodes"]],
                                 ["cosets.build", "cosets.union", "geometry.sign_embed"])

    def test_cli_checkpoint_status_export_verify_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            source = Path(temporary) / "catalogue.json"
            source.write_text(json.dumps(specification()))
            directory = Path(temporary) / "run"
            args = ["plan", str(source), "--state-dir", str(directory)]
            self.assertEqual(main(args + ["--max-steps", "1"]), 2)
            self.assertEqual(main(args + ["--resume"]), 0)
            self.assertEqual(main(["plan-status", str(directory)]), 0)
            out = Path(temporary) / "receipt.json"
            export = ["plan-export", str(directory), "--program-id", "3", "--out", str(out)]
            self.assertEqual(main(export), 0)
            self.assertEqual(main(["verify", str(out)]), 0)
            self.assertEqual(main(export), 1)
            with self.assertRaises(ValueError):
                export_receipt(directory, 9999, Path(temporary) / "bad.json")

    def test_invalid_domains_and_unimplemented_operators_fail_before_creating_state(self):
        cases = []
        for raw in ({"values": []}, {"integer_range": {"start": 0, "stop": 5, "step": 0}},
                    {"values": [1], "subsets": {"of": []}}, {"subsets": {"of": [1, 1]}},
                    {"subsets": {"of": [1], "min_size": 2}}, {"values": [float("nan")]}):
            spec = specification()
            spec["catalogue"][0]["domains"] = {"typo": raw}
            cases.append(spec)
        spec = specification();spec["catalogue"][0]["op"] = "geometry.unimplemented";cases.append(spec)
        spec = specification();spec["max_nodes"] = True;cases.append(spec)
        with tempfile.TemporaryDirectory() as temporary:
            for i, spec in enumerate(cases):
                directory = Path(temporary) / str(i)
                with self.subTest(spec=spec), self.assertRaises(ValueError):
                    plan_search(spec, self.context, directory)
                self.assertFalse(directory.exists())


if __name__ == "__main__":
    unittest.main()
