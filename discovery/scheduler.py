"""Deterministic bounded enumeration of explicit construction-program grammars.

This schedules all registered operators uniformly. It exhausts only the finite
alternatives a caller declares; it is not an enumeration of all mathematics.
Each candidate executes in a disposable interpreter with a hard wall timeout.
"""
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

from .core import Artifact, RunContext, canonical_json, digest, integer, json_value
from .runner import registry, run_program, source_fingerprints
from .certification import verify_artifact


def _positive_seconds(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return float(value)


def _fields(value, allowed, required, label):
    if not isinstance(value, dict) or set(value) - set(allowed) or not set(required) <= set(value):
        raise ValueError(f"Invalid {label} fields")


def compile_grammar(grammar, context, max_grammar_choices=100000):
    """Validate and count branches without materializing their Cartesian product.

    Slots have a fixed topological order. Each choice can change the operator,
    its input edges, literal params, and finite parameter-grid axes. All grids
    are visited in sorted parameter-name order; alternatives retain list order.
    """
    integer(max_grammar_choices, "max_grammar_choices")
    _fields(grammar, {"schema_version", "name", "description", "slots", "outputs"},
            {"schema_version", "slots", "outputs"}, "grammar")
    if type(grammar["schema_version"]) is not int or grammar["schema_version"] != 1:
        raise ValueError("Only grammar schema_version 1 is supported")
    for key in ("name", "description"):
        if key in grammar and not isinstance(grammar[key], str):
            raise ValueError(f"Grammar {key} must be a string")
    slots, outputs = grammar["slots"], grammar["outputs"]
    if not isinstance(slots, list) or not slots:
        raise ValueError("Grammar needs nonempty slots")
    context.limit("max_nodes", len(slots))
    if not isinstance(outputs, list) or not outputs or any(not isinstance(x, str) for x in outputs) or len(set(outputs)) != len(outputs):
        raise ValueError("Grammar needs distinct output ids")
    operations, input_kinds = registry()
    compiled, earlier = [], set()
    total_programs, declared_choices = 1, 0
    for slot in slots:
        _fields(slot, {"id", "choices"}, {"id", "choices"}, "slot")
        name, choices = slot["id"], slot["choices"]
        if not isinstance(name, str) or re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", name) is None or name in earlier:
            raise ValueError("Slot ids must be distinct simple identifiers")
        if not isinstance(choices, list) or not choices:
            raise ValueError("Every slot needs nonempty choices")
        variants, count = [], 0
        for choice in choices:
            _fields(choice, {"op", "inputs", "params", "parameter_grid"}, {"op"}, "choice")
            operation = choice["op"]
            if not isinstance(operation, str) or operation not in operations:
                raise ValueError(f"Unsupported grammar operator {operation!r}")
            inputs = choice.get("inputs", [])
            if not isinstance(inputs, list) or any(not isinstance(x, str) or x not in earlier for x in inputs):
                raise ValueError("Choice inputs must refer to earlier slots")
            if len(inputs) != len(input_kinds[operation]):
                raise ValueError(f"Wrong input count for {operation}")
            params, grid = choice.get("params", {}), choice.get("parameter_grid", {})
            if not isinstance(params, dict) or not isinstance(grid, dict):
                raise ValueError("Choice params and parameter_grid must be objects")
            if any(not isinstance(key, str) for key in grid) or set(params) & set(grid):
                raise ValueError("Grid keys must be strings and disjoint from literal params")
            branch_count = 1
            for key in sorted(grid):
                values = grid[key]
                if not isinstance(values, list) or not values:
                    raise ValueError("Every parameter_grid axis must be a nonempty list")
                branch_count *= len(values)
                if branch_count > max_grammar_choices:
                    raise ValueError("Expanded grammar choices exceed max_grammar_choices")
            count += branch_count
            declared_choices += branch_count
            if declared_choices > max_grammar_choices:
                raise ValueError("Expanded grammar choices exceed max_grammar_choices")
            variants.append({"choice": deepcopy(choice), "count": branch_count})
        compiled.append({"id": name, "variants": variants, "count": count})
        total_programs *= count
        earlier.add(name)
    if any(output not in earlier for output in outputs):
        raise ValueError("Grammar output refers to a missing slot")
    # Reject NaN and non-JSON payloads before any worker starts.
    json_value(grammar)
    return compiled, total_programs


def _node_at(slot, branch):
    for variant in slot["variants"]:
        if branch >= variant["count"]:
            branch -= variant["count"]
            continue
        choice = variant["choice"]
        params = deepcopy(choice.get("params", {}))
        grid = choice.get("parameter_grid", {})
        selected = {}
        # Mixed radix decoding makes the final sorted axis vary fastest.
        for key in reversed(sorted(grid)):
            branch, index = divmod(branch, len(grid[key]))
            selected[key] = deepcopy(grid[key][index])
        params.update(selected)
        return {"id": slot["id"], "op": choice["op"],
                "inputs": list(choice.get("inputs", [])), "params": params}
    raise ValueError("Grammar branch out of range")


def program_at(grammar, compiled, index):
    """Decode an ordinal directly, allowing deterministic explicit start offsets."""
    integer(index, "program index", minimum=0)
    branches = []
    remaining = index
    for slot in reversed(compiled):
        remaining, branch = divmod(remaining, slot["count"])
        branches.append(branch)
    if remaining:
        raise ValueError("Program index exceeds grammar")
    branches.reverse()
    program = {"schema_version": 1, "name": grammar.get("name", "finite grammar") + f" candidate {index}",
               "nodes": [_node_at(slot, branch) for slot, branch in zip(compiled, branches)],
               "outputs": list(grammar["outputs"]),
               "provenance": {"grammar_sha256": digest(grammar), "candidate_index": index,
                              "slot_branches": branches, "enumeration": "declared list order; sorted grid axes; last slot varies fastest"}}
    return program


def _worker(payload):
    context = RunContext(Path(payload["workdir"]), payload["budget"])
    receipt = run_program(payload["program"], context)
    summaries = {}
    for name in receipt["outputs"]:
        artifact = Artifact.from_dict(receipt["artifacts"][name])
        if artifact.verification.status == "unverified":
            summaries[name] = {"kind": artifact.kind, "verification_status": "unverified",
                               "facts": {}, "candidate_facts": artifact.facts}
        else:
            report = verify_artifact(artifact, context)
            summaries[name] = {"kind": artifact.kind, "verification_status": report.status,
                               "facts": report.details.get("facts", report.details)}
    result = {"state": "completed", "receipt": receipt, "outputs": summaries}
    encoded = canonical_json(result)
    if len(encoded.encode("utf-8")) > payload["max_receipt_bytes"]:
        raise ValueError("Candidate receipt exceeds max_receipt_bytes")
    return encoded


def _execute(program, context, seconds, max_receipt_bytes):
    payload = {"program": program, "workdir": str(context.workdir.resolve()),
               "budget": context.budget, "max_receipt_bytes": max_receipt_bytes}
    try:
        process = subprocess.run([sys.executable, "-m", "discovery.scheduler", "--worker"],
                                 input=canonical_json(payload), text=True, capture_output=True,
                                 cwd=Path(__file__).resolve().parent.parent,
                                 timeout=seconds, check=False)
    except subprocess.TimeoutExpired:
        return {"state": "timed_out", "reason": "candidate wall-time limit reached"}
    if process.returncode != 0:
        return {"state": "failed", "reason": (process.stderr.strip() or f"worker exit {process.returncode}")[-4000:]}
    try:
        result = json.loads(process.stdout)
    except (ValueError, TypeError):
        return {"state": "failed", "reason": "worker returned invalid JSON"}
    if not isinstance(result, dict) or result.get("state") != "completed" or "receipt" not in result or "outputs" not in result:
        return {"state": "failed", "reason": "worker returned invalid result schema"}
    return result


def search_grammar(grammar, context, output_directory, *, max_programs=100,
                   time_limit_seconds=60, per_program_seconds=30, start_index=0,
                   max_grammar_choices=100000, max_receipt_bytes=16777216,
                   on_progress=None):
    """Enumerate a finite grammar into a fresh directory of ordinary run receipts.

    Candidate failure never removes later branches. A timeout, failure, offset,
    or enumeration cutoff prevents a claim that the entire grammar completed.
    Geometry ranking uses only independently verified exact outputs by dimension.
    """
    integer(max_programs, "max_programs")
    integer(start_index, "start_index", minimum=0)
    integer(max_receipt_bytes, "max_receipt_bytes")
    total_seconds = _positive_seconds(time_limit_seconds, "time_limit_seconds")
    candidate_seconds = _positive_seconds(per_program_seconds, "per_program_seconds")
    compiled, total = compile_grammar(grammar, context, max_grammar_choices)
    if start_index >= total:
        raise ValueError("start_index must identify a declared candidate")
    directory = Path(output_directory)
    # Exclusive directory creation gives one owner to all receipts and summary.
    directory.mkdir(parents=True, exist_ok=False)
    started, attempts, best = time.monotonic(), [], {}
    stop_reason = "grammar_exhausted"
    upper = min(total, start_index + max_programs)
    for index in range(start_index, upper):
        remaining = total_seconds - (time.monotonic() - started)
        if remaining <= 0:
            stop_reason = "time_limit_reached"
            break
        program = program_at(grammar, compiled, index)
        if on_progress:
            on_progress({"candidate_index": index, "state": "started", "program_sha256": digest(program)})
        candidate_started = time.monotonic()
        result = _execute(program, context, min(candidate_seconds, remaining), max_receipt_bytes)
        record = {"candidate_index": index, "program_sha256": digest(program),
                  "slot_branches": program["provenance"]["slot_branches"],
                  "state": result["state"], "seconds": round(time.monotonic() - candidate_started, 6)}
        if result["state"] == "completed":
            receipt, summaries = result["receipt"], result["outputs"]
            filename = f"candidate-{index:08d}.json"
            with (directory / filename).open("x") as stream:
                stream.write(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
            record.update(receipt=filename, receipt_sha256=digest(receipt), outputs=summaries)
            for name, output in summaries.items():
                facts = output["facts"]
                if output["verification_status"] != "exact" or output["kind"] not in ("spherical_code", "implicit_spherical_code", "gram_spherical_code", "quadratic_spherical_code", "weighted_quadratic_spherical_code"):
                    continue
                dimension, count = facts.get("dimension"), facts.get("cardinality")
                if type(dimension) is not int or type(count) is not int:
                    continue
                key = str(dimension)
                if key not in best or count > best[key]["cardinality"]:
                    best[key] = {"dimension": dimension, "cardinality": count, "receipt": filename,
                                 "output": name, "candidate_index": index}
        else:
            record["reason"] = result.get("reason", "worker failed")
        attempts.append(record)
        if on_progress:
            on_progress(record)
    if stop_reason != "time_limit_reached" and start_index + len(attempts) < total:
        stop_reason = "max_programs_reached" if len(attempts) == max_programs else "time_limit_reached"
    failures = sum(record["state"] != "completed" for record in attempts)
    all_enumerated = start_index == 0 and len(attempts) == total
    summary = {"schema_version": 1, "kind": "construction_grammar_search", "grammar": deepcopy(grammar),
               "grammar_sha256": digest(grammar), "source_sha256": source_fingerprints(),
               "budget": context.budget,
               "limits": {"max_programs": max_programs, "time_limit_seconds": total_seconds,
                          "per_program_seconds": candidate_seconds, "start_index": start_index,
                          "max_grammar_choices": max_grammar_choices, "max_receipt_bytes": max_receipt_bytes},
               "declared_programs": total, "attempted_programs": len(attempts),
               "completed_programs": len(attempts) - failures, "failed_or_timed_out_programs": failures,
               "enumeration_complete": all_enumerated,
               "status": "complete" if all_enumerated and not failures else "incomplete",
               "stop_reason": stop_reason, "next_index": start_index + len(attempts),
               "seconds": round(time.monotonic() - started, 6), "attempts": attempts,
               "best_exact_by_dimension": best,
               "scope": "finite declared grammar only; failures and resource bounds are not mathematical exclusions; no global optimality or SOTA claim",
               "enumeration": "slot order, choice list order, sorted parameter axes; last slot/axis varies fastest; no dimension or family eligibility filter"}
    with (directory / "search.json").open("x") as stream:
        stream.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def _worker_main():
    # Private transport entry point: the only execution is the existing registry.
    # Library operators may print diagnostic text; keep stdout JSON-only.
    import contextlib
    try:
        payload = json.load(sys.stdin)
        with contextlib.redirect_stdout(sys.stderr):
            encoded = _worker(payload)
        sys.stdout.write(encoded)
        return 0
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    if sys.argv[1:] != ["--worker"]:
        raise SystemExit("Use python -m discovery search for the public scheduler")
    raise SystemExit(_worker_main())
