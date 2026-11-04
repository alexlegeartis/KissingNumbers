"""Current operator contracts, mixed input types and checkpoint integration."""
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
from discovery.core import RunContext, canonical_json
from discovery.certification import verify_artifact
from discovery.planner import (compile_catalogue, extension_at, extension_counts,
                               operator_catalogue, plan_search, read_ledger)
from discovery.planner import _domain, MAX_SUBSET_ITEMS, MAX_DOMAIN_COUNT_BITS
from discovery.runner import output_artifact, registry, source_fingerprints
from discovery.scheduler import search_grammar
from discovery.tests.test_planner import direct_execute, rows, specification


def mixed_spec():
    return {"schema_version": 1, "max_nodes": 3,
            "output_kinds": ["gram_spherical_code"], "catalogue": [
                {"op": "geometry.integer_code", "params": {"coordinates": [[1]], "norm_squared": 1}},
                {"op": "geometry.gram_code", "params": {"gram": [[1]]}},
                {"op": "gluing.union", "params": {"cross_gram": [[0]]}}]}


def quadratic_spec():
    return {"schema_version": 1, "max_nodes": 2,
            "output_kinds": ["quadratic_spherical_code"], "catalogue": [
                {"op": "geometry.integer_code", "params": {"coordinates": [[1], [-1]], "norm_squared": 1}},
                {"op": "signed_layers.coupled_lift", "params": {
                    "equator_scale": 2, "cap_rational": [], "cap_radical": [],
                    "heights": [-2, 2], "pole_amplitude": 2}}]}


class PlannerIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.context = RunContext(Path("/tmp"))

    def test_current_new_operator_output_kinds_and_cli_visibility(self):
        expected = {
            "lattice.root_system": "spherical_code", "geometry.gram_code": "gram_spherical_code",
            "geometry.section": "spherical_code", "gluing.union": "gram_spherical_code",
            "groups.real_stabilizer_sector": "implicit_spherical_code",
            "selection.isometric_exchange": "spherical_code", "selection.coupled_union": "spherical_code",
            "signed_layers.conditional_swap": "spherical_code",
            "signed_layers.coupled_lift": "quadratic_spherical_code",
            "signed_layers.free_zone_join": "spherical_code",
            "search.coset_assembly": "search_result"}
        catalogue = operator_catalogue()
        for op, kind in expected.items():
            with self.subTest(op=op):
                self.assertEqual(catalogue[op]["output_kinds"], [kind])
                self.assertTrue(catalogue[op]["planner_supported"])
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.assertEqual(main(["operators"]), 0)
        public = {row["op"]: row for row in json.loads(stream.getvalue())}
        self.assertEqual(public["gluing.union"]["inputs"],
                         [["spherical_code", "gram_spherical_code"]] * 2)
        self.assertTrue(public["gluing.union"]["planner_supported"])

    def test_union_input_intersections_count_each_producer_once(self):
        entry = compile_catalogue(mixed_spec(), self.context)[-1]
        kinds = [["spherical_code"], ["gram_spherical_code"],
                 ["spherical_code", "gram_spherical_code"], ["implicit_spherical_code"]]
        nodes = [{"id": f"n{i}"} for i in range(len(kinds))]
        expected = list(itertools.product(range(3), repeat=2))
        counts = extension_counts([entry], kinds)[0]
        self.assertEqual((counts["typed_extensions"], counts["type_incompatible_extensions"]), (9, 7))
        actual = [extension_at([entry], nodes, kinds, i)[0][-1]["inputs"] for i in range(9)]
        self.assertEqual(actual, [[f"n{i}" for i in pair] for pair in expected])
        with self.assertRaises(ValueError):
            extension_at([entry], nodes, kinds, 9)

    def test_complete_mixed_language_matches_independent_enumerator_and_exact_replay(self):
        spec = mixed_spec()
        entries = compile_catalogue(spec, self.context)
        frontier = [([], [])]
        expected = []
        for depth in range(1, 4):
            following = []
            for nodes, kinds in frontier:
                for entry in entries:
                    pools = []
                    for need in entry.inputs:
                        allowed = {need} if isinstance(need, str) else set(need)
                        pools.append([i for i, produced in enumerate(kinds) if allowed.intersection(produced)])
                    for binding in itertools.product(*pools):
                        child = nodes + [{"id": f"n{len(nodes)}", "op": entry.op,
                                          "params": entry.params, "inputs": [f"n{i}" for i in binding]}]
                        child_kinds = kinds + [list(entry.outputs)]
                        if depth < 3 or "gram_spherical_code" in entry.outputs:
                            expected.append(canonical_json(child))
                        following.append((child, child_kinds))
            frontier = following
        def evaluate(program, *args):
            try:
                return direct_execute(program, *args)
            except ValueError as error:
                return {"state": "failed", "reason": str(error)}
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=evaluate):
            directory = Path(temporary) / "run"
            result = plan_search(spec, self.context, directory, max_steps=1000)
            self.assertEqual([canonical_json(nodes) for _, nodes, _ in rows(directory)], expected)
            self.assertTrue(result["enumeration_complete"])
            # Some statically compatible gluing programs have the wrong cross-Gram
            # shape. Their explicit failures must not discard any descendants.
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["best_exact_by_dimension"]["2"]["cardinality"], 2)
            with sqlite3.connect(directory / "planner.sqlite3") as db:
                receipts = [json.loads(row[0])["receipt"] for row in
                            db.execute("SELECT result FROM programs WHERE eval_state='completed'")]
            mixed = [r for r in receipts if len(r["program"]["nodes"]) == 3 and
                     [n["op"] for n in r["program"]["nodes"][:2]] ==
                     ["geometry.integer_code", "geometry.gram_code"] and
                     r["program"]["nodes"][-1]["inputs"] == ["n0", "n1"]]
            self.assertEqual(len(mixed), 1)
            self.assertEqual(verify_artifact(output_artifact(mixed[0]), self.context).status, "exact")

    def test_quadratic_outputs_rank_in_planner_scheduler_and_compare(self):
        spec = quadratic_spec()
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            directory = Path(temporary) / "plan"
            result = plan_search(spec, self.context, directory)
            self.assertEqual(result["status"], "complete_within_declared_scope")
            self.assertEqual(result["best_exact_by_dimension"]["2"]["cardinality"], 4)
            grammar = {"schema_version": 1, "slots": [
                {"id": "n0", "choices": [spec["catalogue"][0]]},
                {"id": "n1", "choices": [{**spec["catalogue"][1], "inputs": ["n0"]}]}],
                "outputs": ["n1"]}
            scheduler_dir = Path(temporary) / "search"
            scheduled = search_grammar(grammar, self.context, scheduler_dir)
            self.assertEqual(scheduled["best_exact_by_dimension"]["2"]["cardinality"], 4)
            receipt = scheduler_dir / scheduled["attempts"][0]["receipt"]
            catalog = Path(temporary) / "catalog.json"
            catalog.write_text(json.dumps({"comparisons": [{"dimension": 2, "cardinality": 6}]}))
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                self.assertEqual(main(["compare", str(receipt), "--catalog", str(catalog)]), 0)
            self.assertEqual(json.loads(stream.getvalue())["difference"], -2)

    def test_unsupported_contracts_stay_visible_and_reject_before_state_creation(self):
        operations, inputs = registry()
        for invalid_input in (None, [()], [42]):
            new_ops, new_inputs = dict(operations), dict(inputs)
            new_ops["future.unsupported"] = lambda *args: None
            new_inputs["future.unsupported"] = invalid_input
            with self.subTest(invalid_input=invalid_input), patch("discovery.planner.registry", return_value=(new_ops, new_inputs)):
                contract = operator_catalogue()["future.unsupported"]
                self.assertTrue(contract["implemented"])
                self.assertFalse(contract["planner_supported"])
                self.assertIn("missing or invalid OUTPUT_KINDS contract", contract["unsupported_contract_reasons"])
                self.assertIn("invalid INPUT_KINDS contract", contract["unsupported_contract_reasons"])
                spec = specification(1)
                spec["catalogue"] = [{"op": "future.unsupported"}]
                with tempfile.TemporaryDirectory() as temporary:
                    directory = Path(temporary) / "run"
                    with self.assertRaisesRegex(ValueError, "unsupported by planner contracts"):
                        plan_search(spec, self.context, directory)
                    self.assertFalse(directory.exists())

    def test_contract_change_blocks_resume_but_historical_ledger_remains_readable(self):
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", side_effect=direct_execute):
            directory = Path(temporary) / "run"
            before = plan_search(mixed_spec(), self.context, directory, max_steps=1)
            changed = deepcopy(operator_catalogue())
            changed["gluing.union"]["input_kinds"] = ["spherical_code", "spherical_code"]
            with patch("discovery.planner.operator_catalogue", return_value=changed):
                self.assertEqual(read_ledger(directory), before)
                with self.assertRaisesRegex(ValueError, "source fingerprints, contracts"):
                    plan_search(mixed_spec(), self.context, directory, resume=True)

    def test_legacy_checkpoint_ledger_readable_but_incompatible_resume_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "run"
            plan_search(specification(), self.context, directory, max_steps=1)
            with sqlite3.connect(directory / "planner.sqlite3") as db:
                identity = json.loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0])
                identity["checkpoint_schema_version"] = 1
                db.execute("UPDATE metadata SET value=? WHERE key='identity'", (canonical_json(identity),))
                db.execute("DELETE FROM metadata WHERE key='expansion_clock'")
            ledger = read_ledger(directory)
            self.assertEqual(ledger["scheduling"]["policy"], "legacy prefix ordering")
            self.assertIsNone(ledger["scheduling"]["productive_expansions"])
            with self.assertRaisesRegex(ValueError, "same specification"):
                plan_search(specification(), self.context, directory, resume=True)

    def test_live_source_edits_and_mismatched_worker_sources_fail_closed(self):
        for change in ("live_source", "worker_receipt"):
            saved_sources = source_fingerprints()
            changed = [False]
            def fingerprints():
                return {"modified.py": "changed"} if changed[0] else saved_sources
            def execute(program, *args):
                result = direct_execute(program, *args)
                if change == "live_source":
                    changed[0] = True
                else:
                    result["receipt"]["source_sha256"] = {"modified.py": "changed"}
                return result
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temporary, \
                 patch("discovery.planner.source_fingerprints", side_effect=fingerprints), \
                 patch("discovery.planner._execute", side_effect=execute):
                directory = Path(temporary) / "run"
                result = plan_search(specification(), self.context, directory)
                self.assertEqual(result["stop_reason"], "source_or_contract_changed_during_invocation")
                self.assertEqual(result["status"], "incomplete")
                self.assertFalse(result["enumeration_complete"])
                self.assertEqual(result["best_exact_by_dimension"], {})
                self.assertEqual(result["evaluation_attempts"], {"failed": 1})
                self.assertIn("result discarded", result["unresolved_programs"][0]["reason"])

    def test_compiler_preflights_subset_work_and_count_bits_before_checkpoint(self):
        with patch("discovery.planner.math.comb", side_effect=AssertionError("comb reached")):
            with self.assertRaisesRegex(ValueError, "subset-items ceiling"):
                _domain({"subsets": {"of": list(range(MAX_SUBSET_ITEMS + 1))}})
        spec = specification(1)
        spec["catalogue"][0]["domains"] = {
            f"subset{i}": {"subsets": {"of": list(range(MAX_SUBSET_ITEMS))}} for i in range(4)}
        with patch("discovery.planner.math.comb", return_value=1) as comb:
            with self.assertRaisesRegex(ValueError, "subset-work ceiling"):
                compile_catalogue(spec, self.context)
            self.assertEqual(comb.call_count, 3 * (MAX_SUBSET_ITEMS + 1))
        for domains in (
                {"huge": {"integer_range": {"start": 0, "stop": 1 << MAX_DOMAIN_COUNT_BITS}}},
                {key: {"integer_range": {"start": 0, "stop": 1 << 3000}} for key in ("a", "b")}):
            spec["catalogue"][0]["domains"] = domains
            with self.subTest(domains=list(domains)), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary) / "run"
                with self.assertRaisesRegex(ValueError, "ceiling"):
                    plan_search(spec, self.context, directory)
                self.assertFalse(directory.exists())
        lazy = _domain({"integer_range": {"start": 0, "stop": 10**50}})
        self.assertEqual((lazy.count, lazy.at(10**40)), (10**50, 10**40))

    def test_registered_coset_search_generates_candidates_without_ranking_forecast(self):
        spec = {"schema_version": 1, "max_nodes": 1, "output_kinds": ["search_result"],
                "quality_targets": [{"dimension": 12, "cardinality": 841}], "catalogue": [
                    {"op": "search.coset_assembly", "params": {
                        "generator": [[1, 0, 1], [0, 1, 1]], "parity_checks": [[1, 0], [0, 1]],
                        "modulus": 8, "target_distance": 4, "required_cosets": 5, "seed": 4,
                        "offsets_per_family": 64, "information_sets": 4, "refinement_sets": 8,
                        "max_vertices": 32, "max_cliques": 16, "families": ["uniform"]}}]}
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "run"
            result = plan_search(spec, self.context, directory)
            self.assertEqual(result["status"], "complete_within_declared_scope")
            self.assertEqual(result["best_exact_by_dimension"], {})
            self.assertEqual(result["best_observed_candidates_by_dimension"], {})
            self.assertFalse(result["target_quality_passed"])
            with sqlite3.connect(directory / "planner.sqlite3") as db:
                receipt = json.loads(db.execute("SELECT result FROM programs WHERE eval_state='completed'").fetchone()[0])["receipt"]
            artifact = output_artifact(receipt)
            self.assertTrue(artifact.data["result"]["candidates"])
            self.assertEqual(artifact.data["result"]["forecast"]["total_cardinality_if_certified"], 80)
            self.assertEqual(artifact.facts, {})
            self.assertEqual(verify_artifact(artifact, self.context).status, "unverified")

    def test_quality_target_dimension_changes_ordering_only_not_enumerated_scope(self):
        spec = mixed_spec()
        spec["strategy"] = "adaptive"
        spec["quality_targets"] = [{"dimension": 96, "cardinality": 142162362560}]
        failure = {"state": "failed", "reason": "deliberate execution failure"}
        with tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute", return_value=failure):
            first = Path(temporary) / "first"
            second = Path(temporary) / "second"
            result = plan_search(spec, self.context, first, max_steps=1000)
            spec.pop("quality_targets")
            without_target = plan_search(spec, self.context, second, max_steps=1000)
            self.assertEqual({canonical_json(nodes) for _, nodes, _ in rows(first)},
                             {canonical_json(nodes) for _, nodes, _ in rows(second)})
            self.assertTrue(result["enumeration_complete"])
            self.assertFalse(result["target_quality_passed"])
            self.assertEqual(result["generated_programs"], without_target["generated_programs"])

    def test_adaptive_interleaves_constructor_choices_and_reaches_two_node_search_under_small_cap(self):
        for depth in (2, 3):
            spec = {"schema_version": 1, "max_nodes": depth, "strategy": "adaptive",
                    "output_kinds": ["spherical_code"], "catalogue": [
                        {"op": "groups.weil", "params": {"rank": 1, "parity": -1},
                         "domains": {"p": {"values": [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]}}},
                        {"op": "groups.search", "params": {"trials": 0}}]}
            with self.subTest(depth=depth), tempfile.TemporaryDirectory() as temporary, patch("discovery.planner._execute",
                    return_value={"state": "failed", "reason": "deliberate test failure"}) as execute:
                directory = Path(temporary) / "run"
                result = plan_search(spec, self.context, directory, max_programs=5)
                self.assertEqual(result["stop_reason"], "program_limit_reached")
                self.assertGreater(result["evaluation_attempts"].get("failed", 0), 0)
                self.assertGreater(result["program_states"].get("not_target", 0), 0)
                self.assertEqual(result["status"], "incomplete")
                self.assertGreater(execute.call_count, 0)
                with sqlite3.connect(directory / "planner.sqlite3") as db:
                    self.assertLess(db.execute("SELECT expansion_cursor FROM programs WHERE id=0").fetchone()[0], 10)
                    evaluated = db.execute("SELECT min(depth) FROM programs WHERE eval_state='failed'").fetchone()[0]
                    self.assertEqual(evaluated, depth)
                limits = {"max_programs": 1000, "max_steps": 10000, "max_evaluations": 1000}
                resumed = plan_search(spec, self.context, directory, resume=True, **limits)
                full_directory = Path(temporary) / "full"
                full = plan_search(spec, self.context, full_directory, **limits)
                self.assertTrue(resumed["enumeration_complete"])
                self.assertTrue(full["enumeration_complete"])
                self.assertEqual(rows(directory), rows(full_directory))
                self.assertEqual(resumed["program_states"], full["program_states"])


if __name__ == "__main__":
    unittest.main()
