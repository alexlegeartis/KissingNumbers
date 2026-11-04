"""Finite grammar coverage, process isolation, and exact-result ranking controls."""
import contextlib
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from discovery.__main__ import main
from discovery.core import RunContext
from discovery.runner import output_artifact
from discovery.certification import verify_artifact
from discovery.scheduler import compile_grammar, program_at, search_grammar, _execute


def grammar():
    return {"schema_version": 1, "name": "small exact controls", "slots": [
        {"id": "code", "choices": [{"op": "geometry.integer_code", "params": {"norm_squared": 1},
         "parameter_grid": {"coordinates": [[[1, 0], [-1, 0]], [[1, 0], [-1, 0], [0, 1], [0, -1]]]}}]}],
        "outputs": ["code"]}


class SchedulerTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def test_deterministic_grid_and_input_dag_choices(self):
        spec = {"schema_version": 1, "slots": [
            {"id": "a", "choices": [{"op": "groups.weil", "parameter_grid": {"rank": [1, 2], "parity": [1, -1]}}]},
            {"id": "b", "choices": [{"op": "groups.from_generators", "params": {"generators_real": [[[1]]]}}]},
            {"id": "out", "choices": [{"op": "groups.search", "inputs": ["a"]}, {"op": "groups.search", "inputs": ["b"]}]}],
            "outputs": ["out"]}
        compiled, count = compile_grammar(spec, self.context)
        self.assertEqual(count, 8)
        programs = [program_at(spec, compiled, i) for i in range(count)]
        self.assertEqual([p["nodes"][-1]["inputs"] for p in programs[:2]], [["a"], ["b"]])
        self.assertEqual([p["nodes"][0]["params"] for p in programs[::2]],
                         [{"parity": 1, "rank": 1}, {"parity": 1, "rank": 2}, {"parity": -1, "rank": 1}, {"parity": -1, "rank": 2}])
        self.assertEqual(programs, [program_at(spec, compiled, i) for i in range(count)])

    def test_all_registered_operators_are_grammar_eligible(self):
        from discovery.runner import registry
        operations, kinds = registry()
        # Each operator is admitted by the grammar without dimension/family
        # dispatch. Producer type compatibility is checked by the ordinary DAG.
        for operation in operations:
            slots = []
            for i in range(len(kinds[operation])):
                slots.append({"id": f"source{i}", "choices": [{"op": "groups.weil"}]})
            slots.append({"id": "out", "choices": [{"op": operation, "inputs": [f"source{i}" for i in range(len(kinds[operation]))]}]})
            spec = {"schema_version": 1, "slots": slots, "outputs": ["out"]}
            with self.subTest(operation=operation):
                self.assertEqual(compile_grammar(spec, self.context)[1], 1)

    def test_real_workers_exact_receipts_and_dimension_ranking(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "search"
            result = search_grammar(grammar(), self.context, directory, max_programs=5, time_limit_seconds=20)
            self.assertEqual(result["status"], "complete")
            self.assertTrue(result["enumeration_complete"])
            self.assertEqual(result["attempted_programs"], 2)
            self.assertEqual(result["best_exact_by_dimension"]["1"]["cardinality"], 2)
            self.assertEqual(result["best_exact_by_dimension"]["2"]["cardinality"], 4)
            for row in result["attempts"]:
                receipt = json.loads((directory / row["receipt"]).read_text())
                self.assertEqual(verify_artifact(output_artifact(receipt), self.context).status, "exact")
            with self.assertRaises(FileExistsError):
                search_grammar(grammar(), self.context, directory)

    def test_count_cutoff_and_offset_are_explicitly_incomplete(self):
        with tempfile.TemporaryDirectory() as temporary:
            first = search_grammar(grammar(), self.context, Path(temporary) / "first", max_programs=1)
            self.assertEqual((first["status"], first["stop_reason"], first["next_index"]), ("incomplete", "max_programs_reached", 1))
            resumed = search_grammar(grammar(), self.context, Path(temporary) / "resumed", start_index=1)
            self.assertEqual(resumed["status"], "incomplete")
            self.assertEqual(resumed["stop_reason"], "grammar_exhausted")
            self.assertFalse(resumed["enumeration_complete"])

    def test_rejected_parameters_do_not_silently_exclude_or_stop_branches(self):
        spec = grammar()
        spec["slots"][0]["choices"].insert(0, {"op": "geometry.integer_code", "params": {"coordinates": [[1]], "norm_squared": 99}})
        with tempfile.TemporaryDirectory() as temporary:
            result = search_grammar(spec, self.context, Path(temporary) / "search")
            self.assertEqual(result["status"], "incomplete")
            self.assertTrue(result["enumeration_complete"])
            self.assertEqual([row["state"] for row in result["attempts"]], ["failed", "completed", "completed"])
            self.assertIn("norm", result["attempts"][0]["reason"].lower())
            self.assertEqual(result["completed_programs"], 2)

    def test_timeouts_and_receipt_size_limits(self):
        compiled, _ = compile_grammar(grammar(), self.context)
        program = program_at(grammar(), compiled, 0)
        self.assertEqual(_execute(program, self.context, 0.000001, 1000000)["state"], "timed_out")
        self.assertEqual(_execute(program, self.context, 10, 1)["state"], "failed")
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.scheduler._execute", return_value={"state": "timed_out", "reason": "timeout"}):
            result = search_grammar(grammar(), self.context, Path(temporary) / "search")
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["failed_or_timed_out_programs"], 2)
            self.assertEqual(result["best_exact_by_dimension"], {})

    def test_total_deadline_prevents_further_worker_launches(self):
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.scheduler.time.monotonic", side_effect=[0, 1, 1]), patch("discovery.scheduler._execute") as execute:
            result = search_grammar(grammar(), self.context, Path(temporary) / "search", time_limit_seconds=0.5)
            execute.assert_not_called()
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["stop_reason"], "time_limit_reached")
            self.assertEqual(result["attempted_programs"], 0)

    def test_numerical_geometry_is_not_promoted_or_exact_ranked(self):
        spec = {"schema_version": 1, "slots": [
            {"id": "representation", "choices": [{"op": "groups.from_generators", "params": {"generators_real": [[[-1]]]}}]},
            {"id": "candidate", "choices": [{"op": "groups.search", "inputs": ["representation"],
             "params": {"trials": 0, "phase_orders": [1], "phase_steps": 1}}]}], "outputs": ["candidate"]}
        with tempfile.TemporaryDirectory() as temporary:
            result = search_grammar(spec, self.context, Path(temporary) / "search")
            self.assertEqual(result["status"], "complete")
            output = result["attempts"][0]["outputs"]["candidate"]
            self.assertEqual(output["kind"], "spherical_code")
            self.assertEqual(output["verification_status"], "numerical")
            self.assertEqual(result["best_exact_by_dimension"], {})

    def test_invalid_specs_and_grammar_allocation_guard(self):
        cases = []
        for mutation in ("unsupported", "forward", "axis", "duplicate", "extra", "version"):
            spec = deepcopy(grammar())
            choice = spec["slots"][0]["choices"][0]
            if mutation == "unsupported": choice["op"] = "os.system"
            elif mutation == "forward": choice.update(op="groups.search", inputs=["code"])
            elif mutation == "axis": choice["parameter_grid"]["coordinates"] = []
            elif mutation == "duplicate": spec["slots"].append(deepcopy(spec["slots"][0]))
            elif mutation == "extra": spec["typo"] = 1
            else: spec["schema_version"] = True
            cases.append(spec)
        for spec in cases:
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                compile_grammar(spec, self.context)
        with self.assertRaises(ValueError):
            compile_grammar(grammar(), self.context, max_grammar_choices=1)

    def test_cli_search_then_existing_verify_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            spec = Path(temporary) / "grammar.json"
            spec.write_text(json.dumps(grammar()))
            directory = Path(temporary) / "search"
            self.assertEqual(main(["search", str(spec), "--out-dir", str(directory)]), 0)
            self.assertEqual(main(["verify", str(directory / "candidate-00000000.json")]), 0)
            self.assertEqual(main(["search", str(spec), "--out-dir", str(directory)]), 1)


if __name__ == "__main__":
    unittest.main()
