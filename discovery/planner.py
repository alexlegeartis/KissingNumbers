"""Typed, finite DAG synthesis with a transactional, resumable search frontier.

The catalogue declares finite parameter domains, not graph slots or edges.
Programs are topologically labelled DAGs; equivalent/redundant programs are not
quotiented out. Only static kind incompatibility is pruned. Failed evaluations
never prune descendants. SQLite is a checkpoint, not mathematical authority.
"""
from copy import deepcopy
from dataclasses import dataclass
import fcntl
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

from .core import Artifact, RunContext, canonical_json, digest, integer, json_value
from .runner import registry, run_program, source_fingerprints
from .certification import verify_artifact

# Compiler limits apply before the invocation execution timer starts. They
# bound the finite language description, not any mathematical construction.
MAX_CATALOGUE_ENTRIES = 256
MAX_PARAMETER_AXES = 256
MAX_DOMAIN_VALUES = 4096
MAX_SUBSET_ITEMS = 512
MAX_SUBSET_COMPILE_WORK = 1000000
MAX_DOMAIN_COUNT_BITS = 4096

# Conservative contracts for current operators, including conditional returns.
# New backends may publish OUTPUT_KINDS alongside INPUT_KINDS in their module.
_OUTPUT_KINDS = {
    "groups.weil": ("unitary_representation",),
    "groups.from_generators": ("unitary_representation",),
    "groups.search": ("spherical_code", "search_result"),
    "groups.ganzhinov10": ("spherical_code",),
    "groups.ganzhinov14": ("spherical_code",),
    "rings.binary_factor": ("polynomial_pair",),
    "rings.hensel": ("polynomial_pair",),
    "rings.extend": ("ring_code",),
    "rings.gray": ("binary_code",),
    "rings.verify_dgls": ("binary_code",),
    "codes.translate_union": ("binary_code",),
    "cosets.build": ("syndrome_oracle",),
    "cosets.union": ("binary_code",),
    "geometry.integer_code": ("spherical_code",),
    "geometry.sign_embed": ("implicit_spherical_code",),
    "geometry.orthogonal_union": ("spherical_code",),
    "layers.known32": ("implicit_spherical_code",),
    "layers.signed_supports": ("implicit_spherical_code",),
    "lattice.root_system": ("spherical_code",),
    "geometry.gram_code": ("gram_spherical_code",),
    "geometry.section": ("spherical_code",),
    "gluing.union": ("gram_spherical_code",),
    "groups.real_stabilizer_sector": ("implicit_spherical_code",),
    "selection.isometric_exchange": ("spherical_code",),
    "selection.coupled_union": ("spherical_code",),
    "signed_layers.conditional_swap": ("spherical_code",),
    "signed_layers.coupled_lift": ("quadratic_spherical_code",),
    "signed_layers.free_zone_join": ("spherical_code",),
}


def operator_catalogue():
    """Describe executable contracts; never register roadmap placeholders."""
    from .runner import MODULES
    import importlib
    operations, inputs = registry()
    outputs = dict(_OUTPUT_KINDS)
    for name in MODULES:
        module = importlib.import_module(f"discovery.engines.{name}")
        outputs.update(getattr(module, "OUTPUT_KINDS", {}))
    result = {}
    for op in sorted(operations):
        kinds = outputs.get(op)
        if not isinstance(kinds, (tuple, list)) or not kinds or any(not isinstance(k, str) or not k for k in kinds):
            # The operator remains visible, but cannot be silently mistyped.
            kinds = None
        required = inputs[op]
        valid_inputs = isinstance(required, (tuple, list)) and all(
            (isinstance(item, str) and bool(item)) or
            (isinstance(item, (tuple, list)) and bool(item) and
             all(isinstance(kind, str) and bool(kind) for kind in item))
            for item in required)
        issues = []
        if kinds is None: issues.append("missing or invalid OUTPUT_KINDS contract")
        if not valid_inputs: issues.append("invalid INPUT_KINDS contract")
        result[op] = {"op": op, "input_kinds": list(required) if isinstance(required, (tuple, list)) else None,
                      "output_kinds": list(kinds) if kinds else None,
                      "family": op.split(".")[0], "implemented": True,
                      "planner_supported": not issues, "unsupported_contract_reasons": issues}

    return result


def _fields(value, allowed, required, label):
    if not isinstance(value, dict) or set(value) - set(allowed) or not set(required) <= set(value):
        raise ValueError(f"Invalid {label} fields")


def _seconds(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return float(value)


@dataclass(frozen=True)
class Domain:
    """Finite axis, decoded by ordinal without allocating integer ranges."""
    values: object = None
    start: int = 0
    step: int = 1
    count: int = 0
    subset_items: object = None
    include: tuple = ()
    sizes: tuple = ()

    def at(self, index):
        if not 0 <= index < self.count:
            raise ValueError("Parameter ordinal out of range")
        if self.subset_items is not None:
            items = self.subset_items
            for size in self.sizes:
                choices = math.comb(len(items), size)
                if index >= choices:
                    index -= choices
                    continue
                selected, first = list(self.include), 0
                for remaining in range(size, 0, -1):
                    for position in range(first, len(items) - remaining + 1):
                        choices = math.comb(len(items) - position - 1, remaining - 1)
                        if index >= choices:
                            index -= choices
                        else:
                            selected.append(items[position])
                            first = position + 1
                            break
                return deepcopy(selected)
        return deepcopy(self.values[index]) if self.values is not None else self.start + self.step * index


def _domain(raw, subset_work=None):
    _fields(raw, {"values", "integer_range", "subsets"}, set(), "parameter domain")
    if len(raw) != 1:
        raise ValueError("Domain needs exactly one of values, integer_range or subsets")
    if "values" in raw:
        values = raw["values"]
        if not isinstance(values, list) or not values:
            raise ValueError("Domain values must be a nonempty list")
        if len(values) > MAX_DOMAIN_VALUES:
            raise ValueError("Planner compiler domain-values ceiling exceeded")
        json_value(values)
        return Domain(values=deepcopy(values), count=len(values))
    if "subsets" in raw:
        subset = raw["subsets"]
        _fields(subset, {"of", "include", "min_size", "max_size", "largest_first"}, {"of"}, "subset domain")
        items, include = subset["of"], subset.get("include", [])
        if not isinstance(items, list) or not isinstance(include, list):
            raise ValueError("Subset of/include must be lists")
        if len(items) + len(include) > MAX_SUBSET_ITEMS:
            raise ValueError("Planner compiler subset-items ceiling exceeded")
        lower = integer(subset.get("min_size", 0), "min_size", minimum=0)
        upper = integer(subset.get("max_size", len(items)), "max_size", minimum=0)
        reverse = subset.get("largest_first", False)
        if upper > len(items) or lower > upper or type(reverse) is not bool:
            raise ValueError("Invalid subset size bounds/order")
        work = (len(items) + len(include) + 1) * (upper - lower + 1)
        if subset_work is None:
            subset_work = [0]
        if subset_work[0] + work > MAX_SUBSET_COMPILE_WORK:
            raise ValueError("Planner compiler subset-work ceiling exceeded")
        subset_work[0] += work
        keys = [canonical_json(v) for v in items + include]
        if len(set(keys)) != len(keys):
            raise ValueError("Subset of/include must be distinct and disjoint")
        sizes = tuple(sorted(range(lower, upper + 1), reverse=reverse))
        return Domain(subset_items=deepcopy(items), include=tuple(deepcopy(include)), sizes=sizes,
                      count=sum(math.comb(len(items), size) for size in sizes))
    interval = raw["integer_range"]
    _fields(interval, {"start", "stop", "step"}, {"start", "stop"}, "integer range")
    start, stop, step = interval["start"], interval["stop"], interval.get("step", 1)
    if any(type(v) is not int for v in (start, stop, step)) or step == 0:
        raise ValueError("Range bounds and nonzero step must be integers")
    if any(abs(v).bit_length() > MAX_DOMAIN_COUNT_BITS for v in (start, stop, step)):
        raise ValueError("Planner compiler integer-bit ceiling exceeded")
    count = max(0, (stop - start + step - (1 if step > 0 else -1)) // step)
    if not count:
        raise ValueError("Parameter range must be nonempty (stop is exclusive)")
    if count.bit_length() > MAX_DOMAIN_COUNT_BITS:
        raise ValueError("Planner compiler count-bit ceiling exceeded")
    return Domain(start=start, step=step, count=count)


@dataclass(frozen=True)
class Entry:
    op: str
    inputs: tuple
    outputs: tuple
    params: dict
    domains: dict
    count: int

    def params_at(self, ordinal):
        selected = deepcopy(self.params)
        for name in reversed(sorted(self.domains)):
            ordinal, index = divmod(ordinal, self.domains[name].count)
            selected[name] = self.domains[name].at(index)
        if ordinal:
            raise ValueError("Parameter ordinal out of range")
        return selected


def compile_catalogue(spec, context):
    """Validate a finite scope; never choose families or parameters by dimension."""
    _fields(spec, {"schema_version", "name", "description", "max_nodes", "output_kinds", "catalogue", "unavailable", "strategy", "quality_targets"},
            {"schema_version", "max_nodes", "output_kinds", "catalogue"}, "planner specification")
    if type(spec["schema_version"]) is not int or spec["schema_version"] != 1:
        raise ValueError("Only planner schema_version 1 is supported")
    for key in ("name", "description"):
        if key in spec and not isinstance(spec[key], str):
            raise ValueError(f"Planner {key} must be a string")
    if spec.get("strategy", "breadth_first") not in ("breadth_first", "adaptive"):
        raise ValueError("strategy must be breadth_first or adaptive")
    depth = integer(spec["max_nodes"], "max_nodes")
    context.limit("max_nodes", depth)
    targets = spec["output_kinds"]
    if not isinstance(targets, list) or not targets or any(not isinstance(k, str) or not k for k in targets) or len(set(targets)) != len(targets):
        raise ValueError("output_kinds must be distinct nonempty kind names")
    declared = spec["catalogue"]
    if not isinstance(declared, list) or not declared:
        raise ValueError("Planner catalogue must be nonempty")
    if len(declared) > MAX_CATALOGUE_ENTRIES:
        raise ValueError("Planner compiler catalogue-entry ceiling exceeded")
    available, entries = operator_catalogue(), []
    subset_work, axis_count = [0], 0
    for raw in declared:
        _fields(raw, {"op", "params", "domains"}, {"op"}, "catalogue entry")
        op = raw["op"]
        if not isinstance(op, str) or op not in available:
            raise ValueError(f"Operator {op!r} is not implemented/registered")
        contract = available[op]
        if not contract["planner_supported"]:
            raise ValueError(f"Operator {op} is unsupported by planner contracts: " + "; ".join(contract["unsupported_contract_reasons"]))
        params, domains = raw.get("params", {}), raw.get("domains", {})
        if not isinstance(params, dict) or not isinstance(domains, dict) or set(params) & set(domains):
            raise ValueError("Literal params and domains must be disjoint objects")
        axis_count += len(domains)
        if axis_count > MAX_PARAMETER_AXES:
            raise ValueError("Planner compiler parameter-axis ceiling exceeded")
        axes, count = {}, 1
        for key in sorted(domains):
            if not isinstance(key, str):
                raise ValueError("Parameter names must be strings")
            axes[key] = _domain(domains[key], subset_work)
            if count.bit_length() + axes[key].count.bit_length() - 1 > MAX_DOMAIN_COUNT_BITS:
                raise ValueError("Planner compiler parameter-product bit ceiling exceeded")
            count *= axes[key].count
            if count.bit_length() > MAX_DOMAIN_COUNT_BITS:
                raise ValueError("Planner compiler parameter-product bit ceiling exceeded")
        entries.append(Entry(op, tuple(contract["input_kinds"]), tuple(contract["output_kinds"]),
                             deepcopy(params), axes, count))
    unavailable = spec.get("unavailable", [])
    if not isinstance(unavailable, list):
        raise ValueError("unavailable must be a list of named capabilities with reasons")
    for item in unavailable:
        _fields(item, {"name", "family", "reason"}, {"name", "family", "reason"}, "unavailable capability")
        if any(not isinstance(v, str) or not v for v in item.values()) or item["name"] in available:
            raise ValueError("Unavailable capabilities need nonempty strings and cannot be executable operators")
    targets = spec.get("quality_targets", [])
    if not isinstance(targets, list):
        raise ValueError("quality_targets must be a list")
    dimensions = set()
    for target in targets:
        _fields(target, {"dimension", "cardinality", "label", "source"}, {"dimension", "cardinality"}, "quality target")
        dimension = integer(target["dimension"], "quality target dimension")
        integer(target["cardinality"], "quality target cardinality")
        if dimension in dimensions or any(not isinstance(target[k], str) for k in ("label", "source") if k in target):
            raise ValueError("Quality targets need distinct dimensions and string labels/sources")
        dimensions.add(dimension)
    json_value(spec)
    return entries


def _bindings(entry, kinds):
    return [[i for i, options in enumerate(kinds)
             if any(kind in options for kind in ((required,) if isinstance(required, str) else required))]
            for required in entry.inputs]


def extension_counts(entries, kinds):
    """Count typed extensions and rejected raw bindings without expanding products."""
    result = []
    for index, entry in enumerate(entries):
        pools = _bindings(entry, kinds)
        typed = math.prod(len(pool) for pool in pools)
        raw = len(kinds) ** len(entry.inputs)
        result.append({"entry": index, "op": entry.op, "parameter_variants": entry.count,
                       "typed_extensions": typed * entry.count,
                       "type_incompatible_extensions": (raw - typed) * entry.count,
                       "no_compatible_inputs": bool(entry.inputs) and not typed})
    return result


def extension_at(entries, nodes, kinds, ordinal):
    """Append one typed node; catalogue, parameter, input tuple order is stable."""
    integer(ordinal, "extension ordinal", minimum=0)
    for entry in entries:
        pools = _bindings(entry, kinds)
        bindings = math.prod(len(pool) for pool in pools)
        count = entry.count * bindings
        if ordinal >= count:
            ordinal -= count
            continue
        parameter_ordinal, binding_ordinal = divmod(ordinal, bindings)
        selected = []
        for pool in reversed(pools):
            binding_ordinal, index = divmod(binding_ordinal, len(pool))
            selected.append(pool[index])
        selected.reverse()
        node = {"id": f"n{len(nodes)}", "op": entry.op,
                "inputs": [f"n{i}" for i in selected], "params": entry.params_at(parameter_ordinal)}
        return nodes + [node], kinds + [list(entry.outputs)]
    raise ValueError("Extension ordinal out of range")


def _program(nodes, spec_digest):
    return {"schema_version": 1, "nodes": nodes, "outputs": [nodes[-1]["id"]],
            "provenance": {"planner_spec_sha256": spec_digest,
                           "search": "typed DAG synthesis; conditional kinds rechecked at runtime"}}


def _kind_distances(entries, targets):
    """Optimistic type reachability guides ordering only, never eligibility."""
    distances = {kind: 0 for kind in targets}
    changed = True
    while changed:
        changed = False
        for entry in entries:
            distance = min((distances.get(k, 1000000) for k in entry.outputs), default=1000000) + 1
            for required in entry.inputs:
                for kind in ((required,) if isinstance(required, str) else required):
                    if distance < distances.get(kind, 1000000):
                        distances[kind] = distance
                        changed = True
    return distances


def _feedback(result, targets=()):
    """Numerical/unverified counts may guide search, but never exact ranking."""
    best = 0.0
    targets = {item["dimension"]: item["cardinality"] for item in targets}
    for output in result.get("outputs", {}).values():
        facts = output.get("candidate_facts", {})
        count = facts.get("cardinality")
        if type(count) is int and count > 0:
            dimension = facts.get("dimension", output.get("candidate_ambient_dimension"))
            if targets:
                if dimension not in targets:
                    continue
                score = count / targets[dimension]
            else:
                score = math.log2(count + 1)
            best = max(best, score)
    return best


class _ResourceLimit(ValueError):
    pass


class _PlannerContext(RunContext):
    def limit(self, key, value):
        if value > self.budget[key]:
            raise _ResourceLimit(f"{key} budget exceeded: {value} > {self.budget[key]}")


def _worker(payload):
    """All facts used for exact ranking come from independent verification."""
    context = _PlannerContext(Path(payload["workdir"]), payload["budget"])
    receipt = run_program(payload["program"], context)
    outputs = {}
    for name in receipt["outputs"]:
        artifact = Artifact.from_dict(receipt["artifacts"][name])
        output = {"kind": artifact.kind, "candidate_facts": artifact.facts,
                  "verification_status": "unverified", "facts": {}}
        if artifact.verification.status != "unverified":
            report = verify_artifact(artifact, context)
            output.update(verification_status=report.status, facts=report.details.get("facts", report.details))
        complex_dimension = artifact.facts.get("complex_dimension")
        if type(complex_dimension) is int and complex_dimension > 0:
            # A diagnostic ambient bound for an unverified complex orbit result,
            # never a replacement for verified dimension/rank facts.
            output["candidate_ambient_dimension"] = 2 * complex_dimension
        diagnostic_keys = ("sampling_limit_reached", "materialization_limit_reached", "closure_failures",
                           "output_status", "limits", "solver")
        output["search_diagnostics"] = {key: artifact.data[key] for key in diagnostic_keys if key in artifact.data}
        output["search_diagnostics"].update({key: artifact.facts[key] for key in ("search_scope",) if key in artifact.facts})
        proof = artifact.data.get("proof", {})
        if isinstance(proof, dict):
            output["search_diagnostics"].update({key: proof[key] for key in ("scope", "reason") if key in proof})
        declared = operator_catalogue()[receipt["program"]["nodes"][-1]["op"]]["output_kinds"]
        if artifact.kind not in declared:
            raise ValueError("Operator output violated its planner OUTPUT_KINDS contract")
        outputs[name] = output
    result = {"state": "completed", "receipt": receipt, "outputs": outputs}
    if len(canonical_json(result).encode()) > payload["max_receipt_bytes"]:
        raise _ResourceLimit("Candidate receipt exceeds max_receipt_bytes")
    return result


def _execute(program, context, seconds, max_receipt_bytes):
    payload = {"program": program, "workdir": str(context.workdir.resolve()),
               "budget": context.budget, "max_receipt_bytes": max_receipt_bytes}
    try:
        process = subprocess.run([sys.executable, "-m", "discovery.planner", "--worker"],
                                 input=canonical_json(payload), text=True, capture_output=True,
                                 cwd=Path(__file__).resolve().parent.parent,
                                 timeout=seconds, check=False)
    except subprocess.TimeoutExpired:
        return {"state": "timed_out", "reason": "candidate wall-time limit reached"}
    if process.returncode:
        return {"state": "failed", "reason": (process.stderr.strip() or f"worker exit {process.returncode}")[-4000:]}
    try:
        result = json.loads(process.stdout)
        if not isinstance(result, dict) or result.get("state") not in ("completed", "failed", "resource_limited"):
            raise ValueError("Invalid result schema")
        if result["state"] == "completed" and ("receipt" not in result or "outputs" not in result):
            raise ValueError("Missing receipt or outputs")
        return result
    except (ValueError, TypeError):
        return {"state": "failed", "reason": "worker returned invalid JSON/result schema"}


_SCHEMA = """
CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE programs (
 id INTEGER PRIMARY KEY, nodes TEXT NOT NULL, kinds TEXT NOT NULL, depth INTEGER NOT NULL,
 eval_state TEXT NOT NULL, result TEXT, expansion_cursor INTEGER NOT NULL DEFAULT 0,
 extension_total TEXT, expansion_done INTEGER NOT NULL DEFAULT 0,
 non_target_leaf_pruned INTEGER NOT NULL DEFAULT 0,
 priority_score REAL NOT NULL DEFAULT 0, kind_distance INTEGER NOT NULL DEFAULT 1000000,
 expansion_turn INTEGER NOT NULL DEFAULT 0);
CREATE TABLE expansions (parent INTEGER NOT NULL, entry INTEGER NOT NULL, counts TEXT NOT NULL,
 PRIMARY KEY(parent, entry));
CREATE TABLE evaluations (id INTEGER PRIMARY KEY, program_id INTEGER NOT NULL,
 state TEXT NOT NULL, result TEXT, seconds REAL);
CREATE TABLE sessions (id INTEGER PRIMARY KEY, limits TEXT NOT NULL, stop_reason TEXT);
"""


class _Checkpoint:
    def __init__(self, directory, identity, resume):
        self.directory = Path(directory)
        if resume:
            if not (self.directory / "planner.sqlite3").is_file():
                raise ValueError("No planner checkpoint to resume")
        else:
            self.directory.mkdir(parents=True, exist_ok=False)
        self.lock = (self.directory / "planner.lock").open("a")
        self.db = None
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise ValueError("Planner checkpoint is already in use") from None
        try:
            self.db = sqlite3.connect(self.directory / "planner.sqlite3")
            self.db.row_factory = sqlite3.Row
            self.db.execute("PRAGMA journal_mode=WAL")
            self.db.execute("PRAGMA synchronous=FULL")
            if not resume:
                self.db.executescript(_SCHEMA)
                with self.db:
                    self.db.execute("INSERT INTO metadata VALUES ('identity', ?)", (canonical_json(identity),))
                    self.db.execute("INSERT INTO metadata VALUES ('expansion_clock', '0')")
                    self.db.execute("INSERT INTO programs(id,nodes,kinds,depth,eval_state,kind_distance) VALUES (0,'[]','[]',0,'not_target',-1)")
            else:
                saved = self.db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
                if not saved or canonical_json(json.loads(saved[0])) != canonical_json(identity):
                    raise ValueError("Resume requires the same specification, source fingerprints, contracts, budgets and workdir")
                with self.db:
                    # Results commit atomically. An interrupted attempt has no
                    # authoritative outcome and is retried, without pruning.
                    self.db.execute("UPDATE evaluations SET state='interrupted', result=? WHERE state='running'",
                                    (canonical_json({"reason": "no committed result; requeued on resume"}),))
                    self.db.execute("UPDATE programs SET eval_state='pending' WHERE eval_state='running'")
                    self.db.execute("UPDATE sessions SET stop_reason='interrupted' WHERE stop_reason IS NULL")
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.db is not None:
            self.db.close()
        self.lock.close()


def _identity(spec, context):
    return {"checkpoint_schema_version": 2, "specification": spec, "spec_sha256": digest(spec),
            "source_sha256": source_fingerprints(), "operator_catalogue": operator_catalogue(),
            "budget": context.budget, "workdir": str(context.workdir.resolve())}


def _implementation_matches(identity):
    """Detect edits across worker boundaries; this is not a hermetic snapshot."""
    return (source_fingerprints() == identity["source_sha256"] and
            canonical_json(operator_catalogue()) == canonical_json(identity["operator_catalogue"]))


def _yield_expansion(db, program_id):
    """Called in the child/cursor transaction; resume preserves scheduling turns."""
    turn = int(db.execute("SELECT value FROM metadata WHERE key='expansion_clock'").fetchone()[0]) + 1
    db.execute("UPDATE metadata SET value=? WHERE key='expansion_clock'", (str(turn),))
    db.execute("UPDATE programs SET expansion_turn=?,kind_distance=CASE WHEN id=0 THEN 1000000 ELSE kind_distance END WHERE id=?",
               (turn, program_id))


def _summary(db):
    identity = json.loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0])
    states = {r[0]: r[1] for r in db.execute("SELECT eval_state, count(*) FROM programs WHERE id>0 GROUP BY eval_state")}
    expansion_rows = db.execute("SELECT expansion_cursor, extension_total FROM programs WHERE expansion_done=0").fetchall()
    known_remaining = sum(int(row[1]) - row[0] for row in expansion_rows if row[1] is not None)
    unexplored_prefixes = sum(row[1] is None for row in expansion_rows)
    failed = sum(states.get(s, 0) for s in ("failed", "timed_out", "resource_limited"))
    pending = sum(states.get(s, 0) for s in ("pending", "running"))
    enumeration_complete = not expansion_rows and not pending
    pruning = {}
    for row in db.execute("SELECT entry, counts FROM expansions"):
        counts = json.loads(row["counts"])
        key = str(row["entry"])
        aggregate = pruning.setdefault(key, {"op": counts["op"], "expanded_prefixes": 0,
                                            "typed_extensions": 0, "type_incompatible_extensions": 0,
                                            "prefixes_without_compatible_inputs": 0})
        aggregate["expanded_prefixes"] += 1
        for field in ("typed_extensions", "type_incompatible_extensions"):
            aggregate[field] += counts[field]
        aggregate["prefixes_without_compatible_inputs"] += counts["no_compatible_inputs"]
    best, candidate_best, limited_outputs = {}, {}, []
    for row in db.execute("SELECT id, result FROM programs WHERE eval_state='completed' ORDER BY id"):
        result = json.loads(row["result"])
        for name, output in result["outputs"].items():
            candidate = output.get("candidate_facts", {})
            candidate_dimension = candidate.get("dimension", output["facts"].get("dimension", output.get("candidate_ambient_dimension")))
            candidate_count = candidate.get("cardinality")
            if output["kind"] in ("spherical_code", "implicit_spherical_code", "gram_spherical_code", "quadratic_spherical_code", "search_result") and type(candidate_dimension) is int and type(candidate_count) is int:
                key = str(candidate_dimension)
                if key not in candidate_best or candidate_count > candidate_best[key]["cardinality"]:
                    candidate_best[key] = {"dimension": candidate_dimension, "cardinality": candidate_count,
                                           "verification_status": output["verification_status"], "program_id": row["id"]}
            diagnostics = output.get("search_diagnostics", {})
            if diagnostics:
                limited_outputs.append({"program_id": row["id"], "output": name, "diagnostics": diagnostics})
            if output["verification_status"] != "exact" or output["kind"] not in ("spherical_code", "implicit_spherical_code", "gram_spherical_code", "quadratic_spherical_code"):
                continue
            facts = output["facts"]
            dimension, count = facts.get("dimension"), facts.get("cardinality")
            if type(dimension) is int and type(count) is int and count > 0:
                key = str(dimension)
                if key not in best or count > best[key]["cardinality"]:
                    best[key] = {"dimension": dimension, "cardinality": count,
                                 "program_id": row["id"], "output": name}
    spec = identity["specification"]
    quality = []
    for target in spec.get("quality_targets", []):
        key = str(target["dimension"])
        count = best.get(key, {}).get("cardinality")
        candidate = candidate_best.get(key)
        quality.append({"target": target, "best_exact_cardinality": count,
                        "best_observed_candidate": candidate,
                        "exact_fraction_of_target": count / target["cardinality"] if count is not None else None,
                        "exact_deficit": max(0, target["cardinality"] - count) if count is not None else None,
                        "candidate_fraction_of_target": candidate["cardinality"] / target["cardinality"] if candidate else None,
                        "candidate_deficit": max(0, target["cardinality"] - candidate["cardinality"]) if candidate else None,
                        "target_met_by_exact_result": count is not None and count >= target["cardinality"]})
    selected = {entry["op"] for entry in spec["catalogue"]}
    last = db.execute("SELECT id, limits, stop_reason FROM sessions ORDER BY id DESC LIMIT 1").fetchone()
    counts = db.execute("SELECT state,count(*) FROM evaluations GROUP BY state").fetchall()
    unresolved = [{"program_id": row["id"], "state": row["eval_state"],
                   "reason": json.loads(row["result"]).get("reason") if row["result"] else None}
                  for row in db.execute("SELECT id,eval_state,result FROM programs WHERE eval_state IN ('failed','timed_out','resource_limited','running','pending') ORDER BY id")]
    return {"schema_version": 1, "kind": "typed_construction_planner", "spec_sha256": identity["spec_sha256"],
            "status": "complete_within_declared_scope" if enumeration_complete and not failed else "incomplete",
            "enumeration_complete": enumeration_complete,
            "stop_reason": last["stop_reason"] if last else "not_started",
            "generated_programs": sum(states.values()), "program_states": states,
            "explored_programs": sum(states.get(s, 0) for s in ("completed", "failed", "timed_out", "resource_limited")),
            "completed_programs": states.get("completed", 0), "unresolved_evaluations": failed + pending,
            "evaluation_attempts": {r[0]: r[1] for r in counts},
            "unresolved_programs": unresolved,
            "unsearched": {"generated_pending_evaluations": pending,
                           "unexpanded_prefixes": len(expansion_rows),
                           "known_immediate_extensions": known_remaining,
                           "prefixes_with_uncounted_extensions": unexplored_prefixes,
                           "total_remaining_programs": 0 if enumeration_complete else None,
                           "descendants_of_frontier": "unknown until expanded" if expansion_rows else "none"},
            "pruning_by_catalogue_entry": pruning, "best_exact_by_dimension": best,
            "pruned_non_target_leaves": db.execute("SELECT sum(non_target_leaf_pruned) FROM programs").fetchone()[0],
            "strategy": spec.get("strategy", "breadth_first"),
            "heuristic": "adaptive: pending outputs, fresh bounded descent, then alternating oldest suspended prefix and kind/quality order; ordering only; no dimension eligibility filter",
            "scheduling": {"checkpoint_schema_version": identity["checkpoint_schema_version"],
                           "policy": "persistent alternating queue/heuristic turns" if identity["checkpoint_schema_version"] >= 2 else "legacy prefix ordering",
                           "productive_expansions": int(db.execute("SELECT value FROM metadata WHERE key='expansion_clock'").fetchone()[0]) if identity["checkpoint_schema_version"] >= 2 else None,
                           "pruned_leaf_variants": "consume steps without advancing the productive-expansion clock"},
            "best_observed_candidates_by_dimension": candidate_best,
            "operator_search_diagnostics": limited_outputs,
            "quality_diagnostics": quality,
            "target_quality_passed": all(item["target_met_by_exact_result"] for item in quality) if quality else None,
            "quality_notice": "enumeration completion is not discovery success; target deficits require better candidates or new structural operators; numerical counts do not certify lower bounds",
            "resource_limits": json.loads(last["limits"]) if last else {},
            "declared_scope": {"max_nodes": spec["max_nodes"], "output_kinds": spec["output_kinds"],
                               "parameter_domains": spec["catalogue"],
                               "registered_but_omitted": sorted(set(identity["operator_catalogue"]) - selected),
                               "unimplemented_capabilities": spec.get("unavailable", []),
                               "beyond_max_nodes": "unsearched", "outside_parameter_domains": "unsearched"},
            "scope": "finite typed DAGs only; static kind pruning is not a mathematical exclusion; no global optimum or SOTA claim",
            "mathematical_exclusion_claim": False,
            "checkpoint_notice": "saved results are receipts, not proof; export and independently verify before use"}


def plan_search(spec, context, output_directory, *, resume=False, retry_unresolved=False, max_steps=1000, max_evaluations=100,
                max_programs=10000, time_limit_seconds=60, per_program_seconds=30,
                max_program_bytes=2097152, max_receipt_bytes=16777216, on_progress=None):
    """Synthesize nodes and input edges, resuming an identical scope transactionally.

    Invocation limits may grow on resume. The catalogue, max_nodes, execution
    budget, contracts and source may not change. Failed attempts are recorded,
    and descendants remain eligible. No family dispatch or dimension pruning.
    """
    entries = compile_catalogue(spec, context)
    if retry_unresolved and not resume:
        raise ValueError("retry_unresolved requires resume")
    for name, value in (("max_steps", max_steps), ("max_evaluations", max_evaluations),
                        ("max_programs", max_programs), ("max_program_bytes", max_program_bytes),
                        ("max_receipt_bytes", max_receipt_bytes)):
        integer(value, name)
    total_seconds = _seconds(time_limit_seconds, "time_limit_seconds")
    candidate_seconds = _seconds(per_program_seconds, "per_program_seconds")
    limits = {"max_steps": max_steps, "max_evaluations": max_evaluations, "max_programs": max_programs,
              "time_limit_seconds": total_seconds, "per_program_seconds": candidate_seconds,
              "max_program_bytes": max_program_bytes, "max_receipt_bytes": max_receipt_bytes,
              "retry_unresolved": retry_unresolved}
    identity = _identity(spec, context)
    checkpoint = _Checkpoint(output_directory, identity, resume)
    db = checkpoint.db
    started, steps, evaluated = time.monotonic(), 0, 0
    spec_hash = digest(spec)
    distances = _kind_distances(entries, spec["output_kinds"])
    try:
        with db:
            session = db.execute("INSERT INTO sessions(limits) VALUES (?)", (canonical_json(limits),)).lastrowid
            if retry_unresolved:
                db.execute("UPDATE programs SET eval_state='pending' WHERE eval_state IN ('failed','timed_out','resource_limited')")
        stop = "frontier_exhausted"
        while True:
            if not _implementation_matches(identity):
                stop = "source_or_contract_changed_during_invocation"
                break
            if spec.get("strategy") == "adaptive":
                clock = int(db.execute("SELECT value FROM metadata WHERE key='expansion_clock'").fetchone()[0])
                prefix = "CASE WHEN eval_state='pending' THEN 0 ELSE 1 END, CASE WHEN expansion_turn=0 THEN 0 ELSE 1 END, "
                order = prefix + ("expansion_turn, kind_distance, priority_score DESC, id" if clock % 2 == 0 else
                                  "kind_distance, priority_score DESC, expansion_turn, id")
            else:
                order = "id"
            row = db.execute(f"SELECT * FROM programs WHERE eval_state='pending' OR expansion_done=0 ORDER BY {order} LIMIT 1").fetchone()
            if row is None:
                break
            remaining = total_seconds - (time.monotonic() - started)
            if remaining <= 0:
                stop = "time_limit_reached"
                break
            if steps >= max_steps:
                stop = "step_limit_reached"
                break
            if row["eval_state"] == "pending":
                if evaluated >= max_evaluations:
                    stop = "evaluation_limit_reached"
                    break
                program = _program(json.loads(row["nodes"]), spec_hash)
                with db:
                    db.execute("UPDATE programs SET eval_state='running' WHERE id=?", (row["id"],))
                    attempt = db.execute("INSERT INTO evaluations(program_id,state) VALUES (?,'running')", (row["id"],)).lastrowid
                if on_progress:
                    on_progress({"program_id": row["id"], "state": "started", "program_sha256": digest(program)})
                candidate_started = time.monotonic()
                result = _execute(program, context, min(candidate_seconds, remaining), max_receipt_bytes)
                implementation_changed = not _implementation_matches(identity)
                receipt_changed = (result["state"] == "completed" and
                                   result["receipt"].get("source_sha256") != identity["source_sha256"])
                if implementation_changed or receipt_changed:
                    # Never rank a result under another implementation's saved
                    # scope. Preserve the unresolved attempt and remaining work.
                    result = {"state": "failed", "reason":
                              "Source fingerprints or operator contracts changed during evaluation, "
                              "or worker receipt does not match checkpoint sources; result discarded"}
                encoded = canonical_json(result)
                with db:
                    db.execute("UPDATE programs SET eval_state=?,result=?,priority_score=? WHERE id=?",
                               (result["state"], encoded, max(row["priority_score"], _feedback(result, spec.get("quality_targets", []))), row["id"]))
                    db.execute("UPDATE evaluations SET state=?,result=?,seconds=? WHERE id=?",
                               (result["state"], canonical_json({k: v for k, v in result.items() if k != "receipt"}), time.monotonic() - candidate_started, attempt))
                steps += 1
                evaluated += 1
                if on_progress:
                    on_progress({"program_id": row["id"], "state": result["state"],
                                 **({"reason": result["reason"]} if "reason" in result else {})})
                if implementation_changed or receipt_changed:
                    stop = "source_or_contract_changed_during_invocation"
                    break
                continue
            nodes, kinds = json.loads(row["nodes"]), json.loads(row["kinds"])
            if row["extension_total"] is None:
                counts = extension_counts(entries, kinds)
                total = sum(count["typed_extensions"] for count in counts)
                with db:
                    db.execute("UPDATE programs SET extension_total=? WHERE id=?", (str(total), row["id"]))
                    db.executemany("INSERT INTO expansions VALUES (?,?,?)",
                                   [(row["id"], count["entry"], canonical_json(count)) for count in counts])
            else:
                total = int(row["extension_total"])
            if row["expansion_cursor"] == total:
                with db:
                    db.execute("UPDATE programs SET expansion_done=1 WHERE id=?", (row["id"],))
                steps += 1
                continue
            child, child_kinds = extension_at(entries, nodes, kinds, row["expansion_cursor"])
            target = bool(set(child_kinds[-1]) & set(spec["output_kinds"]))
            if len(child) == spec["max_nodes"] and not target:
                with db:
                    db.execute("UPDATE programs SET expansion_cursor=expansion_cursor+1, non_target_leaf_pruned=non_target_leaf_pruned+1 WHERE id=?", (row["id"],))
                steps += 1
                continue
            if db.execute("SELECT count(*) FROM programs WHERE id>0").fetchone()[0] >= max_programs:
                stop = "program_limit_reached"
                break
            encoded_nodes = canonical_json(child)
            if len(canonical_json(_program(child, spec_hash)).encode()) > max_program_bytes:
                stop = "program_size_limit_reached"
                break
            with db:
                db.execute("INSERT INTO programs(nodes,kinds,depth,eval_state,expansion_done,priority_score,kind_distance) VALUES (?,?,?,?,?,?,?)",
                           (encoded_nodes, canonical_json(child_kinds), len(child), "pending" if target else "not_target",
                            int(len(child) == spec["max_nodes"]), row["priority_score"],
                            min(distances.get(k, 1000000) for k in child_kinds[-1])))
                db.execute("UPDATE programs SET expansion_cursor=expansion_cursor+1 WHERE id=?", (row["id"],))
                _yield_expansion(db, row["id"])
            steps += 1
        with db:
            db.execute("UPDATE sessions SET stop_reason=? WHERE id=?", (stop, session))
        return _summary(db)
    finally:
        checkpoint.close()


def _read_db(directory):
    path = Path(directory).resolve() / "planner.sqlite3"
    if not path.is_file():
        raise ValueError("Planner checkpoint not found")
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    return db


def read_ledger(directory):
    """Read a consistent snapshot without mutating or resuming the checkpoint."""
    db = _read_db(directory)
    try:
        with db:
            db.execute("BEGIN")
            return _summary(db)
    finally:
        db.close()


def export_receipt(directory, program_id, output):
    """Export one normal runner receipt for the existing verify/compare commands."""
    integer(program_id, "program_id")
    db = _read_db(directory)
    try:
        row = db.execute("SELECT eval_state,result FROM programs WHERE id=?", (program_id,)).fetchone()
        if not row or row["eval_state"] != "completed":
            raise ValueError("Program has no completed receipt to export")
        receipt = json.loads(row["result"])["receipt"]
        with Path(output).open("x") as stream:
            stream.write(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        return receipt
    finally:
        db.close()


def _worker_main():
    import contextlib
    try:
        payload = json.load(sys.stdin)
        with contextlib.redirect_stdout(sys.stderr):
            result = _worker(payload)
    except _ResourceLimit as error:
        result = {"state": "resource_limited", "reason": str(error), "exception_type": type(error).__name__}
    except Exception as error:
        result = {"state": "failed", "reason": str(error), "exception_type": type(error).__name__}
    sys.stdout.write(canonical_json(result))
    return 0


if __name__ == "__main__":
    if sys.argv[1:] != ["--worker"]:
        raise SystemExit("Use python -m discovery plan for the public planner")
    raise SystemExit(_worker_main())
