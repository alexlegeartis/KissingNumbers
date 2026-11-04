"""Bounded exact binary syndrome oracles and unions of their kernel cosets.

Only the quotient is enumerated, never the potentially exponential kernel.
For H with r independent rows, bit i of a syndrome is row i of Hx over F2.
The Cayley BFS certifies all coset leader weights. A separate coordinate DP
certifies the minimum *nonempty* dependency among the columns of H; the empty
word of syndrome zero cannot establish the kernel's nonzero distance.

The quotient has at most 2^16 states. Both the context's max_codewords bound
and this fixed cap apply before allocating any state tables. Tables use packed
integer arrays; work is O(n 2^r), with an additional fixed transition cap.
"""
from array import array
from pathlib import Path

from discovery.core import (Artifact, RunContext, Verification, canonical_json,
                            check_params, integer)

BACKEND = "binary_syndrome_cosets"
MAX_SYNDROME_BITS = 16
MAX_LENGTH = 4096
MAX_DP_TRANSITIONS = 2**24
SYNDROME_CONVENTION = "bit i is row i of Hx modulo 2"
DISTANCE_CONVENTION = "n+1 denotes absence of a nonzero kernel word"


def _context(context):
    return context if context is not None else RunContext(Path("."))


def _matrix(H, context):
    if not isinstance(H, list) or not H:
        raise ValueError("parity_check must be a nonempty list of binary rows")
    r = len(H)
    # Check the exponent before computing 2^r or allocating quotient storage.
    if r > MAX_SYNDROME_BITS:
        raise ValueError(f"syndrome bits exceed fixed cap {MAX_SYNDROME_BITS}")
    states = 1 << r
    context.limit("max_codewords", states)
    if not isinstance(H[0], list) or not H[0]:
        raise ValueError("parity_check rows must be nonempty lists")
    n = len(H[0])
    context.limit("max_dimension", n)
    if n > MAX_LENGTH or n * states > MAX_DP_TRANSITIONS:
        raise ValueError("syndrome DP work exceeds fixed bounded backend limits")
    rows = []
    for row in H:
        if not isinstance(row, list) or len(row) != n:
            raise ValueError("parity_check must be rectangular")
        if any(type(x) is not int or x not in (0, 1) for x in row):
            raise ValueError("parity_check entries must be binary integers")
        rows.append(sum(bit << j for j, bit in enumerate(row)))
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    if len(pivots) != r:
        raise ValueError("parity_check must have full row rank over F2")
    columns = [sum(H[i][j] << i for i in range(r)) for j in range(n)]
    return n, r, states, columns


def _syndrome_weights(columns, states, n):
    """Shortest paths in the additive syndrome Cayley graph."""
    missing = n + 1
    distances = array("H", [missing]) * states
    queue = array("I", [0])
    distances[0] = 0
    head = 0
    # A shortest walk never repeats a coordinate: a repeated column cancels.
    while head < len(queue):
        source = queue[head]
        head += 1
        candidate = distances[source] + 1
        for column in columns:
            target = source ^ column
            if distances[target] == missing:
                distances[target] = candidate
                queue.append(target)
    if len(queue) != states:
        raise ValueError("full-rank syndrome graph was not connected")
    return distances


def _kernel_distance(columns, states, n):
    """Minimum nonempty zero-syndrome support, with each coordinate used once.

    Just before coordinate j is processed, best[s] is the lightest subset of
    earlier coordinates with syndrome s. A nonempty dependency whose final
    coordinate is j therefore has weight best[column_j]+1. Snapshot updates
    prevent using the new coordinate more than once. best[0]=0 is intentional;
    the candidate always includes the new coordinate and is nonempty.
    """
    missing = n + 1
    best = array("H", [missing]) * states
    best[0] = 0
    distance = missing
    for column in columns:
        distance = min(distance, best[column] + 1)
        updated = best[:]
        for source, weight in enumerate(best):
            target = source ^ column
            candidate = weight + 1
            if candidate < updated[target]:
                updated[target] = candidate
        best = updated
    return distance, best


def _oracle(H, context):
    n, r, states, columns = _matrix(H, context)
    weights = _syndrome_weights(columns, states, n)
    kernel_distance, subset_weights = _kernel_distance(columns, states, n)
    if weights != subset_weights:
        raise ValueError("Cayley BFS and independent coordinate DP disagree")
    kernel_size = 1 << (n - r)
    if (kernel_distance == n + 1) != (kernel_size == 1):
        raise ValueError("nonzero kernel distance disagrees with rank-nullity")
    facts = dict(length=n, rank=r, kernel_dimension=n-r,
                 cardinality=kernel_size,
                 min_distance=kernel_distance if kernel_size > 1 else None)
    proof = dict(method="exact syndrome BFS and nonempty coordinate dependency DP",
                 syndrome_min_weights=list(weights),
                 kernel_min_nonzero_weight=kernel_distance,
                 kernel_distance_convention=DISTANCE_CONVENTION,
                 syndrome_states=states, dp_transitions=n*states)
    data = dict(backend=BACKEND, parity_check=[row[:] for row in H],
                syndrome_convention=SYNDROME_CONVENTION, proof=proof)
    report = Verification("exact", proof["method"],
                          dict(facts=facts, implicit=True,
                               distance_scope="kernel; singleton distance is undefined"))
    return Artifact("syndrome_oracle", data, facts, report)


def _same(actual, expected, description):
    # JSON comparison also rejects bool/float counterfeits of integer facts.
    if canonical_json(actual) != canonical_json(expected):
        raise ValueError(f"{description} were modified or are unsupported")


def _rebuild_oracle(data, context):
    if not isinstance(data, dict) or set(data) != {
            "backend", "parity_check", "syndrome_convention", "proof"}:
        raise ValueError("invalid syndrome oracle schema")
    if data["backend"] != BACKEND:
        raise ValueError("unknown syndrome oracle backend")
    rebuilt = _oracle(data["parity_check"], context)
    _same(data, rebuilt.data, "syndrome oracle recipe and proof")
    return rebuilt


def build(inputs, params, context):
    """Build an exact quotient oracle from a full-rank binary parity check."""
    check_params(params, ("parity_check",))
    if inputs:
        raise ValueError("cosets.build takes no inputs")
    if "parity_check" not in params:
        raise ValueError("cosets.build requires parity_check")
    return _oracle(params["parity_check"], _context(context))


def _union(oracle, labels, context):
    if not isinstance(labels, list) or not labels:
        raise ValueError("syndrome_labels must be a nonempty list")
    states = oracle.data["proof"]["syndrome_states"]
    context.limit("max_codewords", len(labels))
    context.limit("max_pairs", len(labels)*(len(labels)-1)//2)
    for label in labels:
        integer(label, "syndrome label", minimum=0)
        if label >= states:
            raise ValueError("syndrome label is outside the quotient")
    if len(set(labels)) != len(labels):
        raise ValueError("syndrome labels must be distinct (cosets must be disjoint)")
    weights = oracle.data["proof"]["syndrome_min_weights"]
    cross = None
    for i, label in enumerate(labels):
        for previous in labels[:i]:
            distance = weights[label ^ previous]
            cross = distance if cross is None else min(cross, distance)
    kernel_distance = oracle.data["proof"]["kernel_min_nonzero_weight"]
    size = oracle.facts["cardinality"] * len(labels)
    distance = min(kernel_distance, cross) if cross is not None else kernel_distance
    facts = dict(length=oracle.facts["length"], cardinality=size,
                 min_distance=distance if size > 1 else None)
    proof = dict(method="exact disjoint syndrome coset union",
                 selected_cosets=len(labels), kernel_size=oracle.facts["cardinality"],
                 cross_min_distance=cross,
                 scope="specified cosets only; no global optimality claim")
    data = dict(backend=BACKEND, oracle=oracle.data,
                syndrome_labels=labels[:], proof=proof)
    return Artifact("binary_code", data, facts,
                    Verification("exact", proof["method"], dict(facts=facts, implicit=True)))


def union(inputs, params, context):
    """Certify exactly the union specified by distinct integer syndromes."""
    check_params(params, ("syndrome_labels",))
    if len(inputs) != 1 or inputs[0].kind != "syndrome_oracle":
        raise ValueError("cosets.union requires one syndrome_oracle")
    if "syndrome_labels" not in params:
        raise ValueError("cosets.union requires syndrome_labels")
    context = _context(context)
    oracle = _rebuild_oracle(inputs[0].data, context)
    _same(inputs[0].facts, oracle.facts, "input oracle facts")
    return _union(oracle, params["syndrome_labels"], context)


def verify_artifact(artifact, context=None):
    """Recompute rank, quotient distances, kernel distance, and union facts."""
    context = _context(context)
    try:
        if artifact.kind == "syndrome_oracle":
            rebuilt = _rebuild_oracle(artifact.data, context)
        elif artifact.kind == "binary_code":
            data = artifact.data
            if set(data) != {"backend", "oracle", "syndrome_labels", "proof"}:
                raise ValueError("invalid syndrome union schema")
            if data["backend"] != BACKEND:
                raise ValueError("unknown syndrome union backend")
            oracle = _rebuild_oracle(data["oracle"], context)
            rebuilt = _union(oracle, data["syndrome_labels"], context)
        else:
            raise ValueError("not a syndrome-coset artifact")
        _same(artifact.data, rebuilt.data, "artifact recipe and proof")
        _same(artifact.facts, rebuilt.facts, "artifact facts")
        return rebuilt.verification
    except (KeyError, TypeError, IndexError, OverflowError) as exc:
        raise ValueError("malformed syndrome-coset artifact") from exc


OPERATORS = {"cosets.build": build, "cosets.union": union}
INPUT_KINDS = {"cosets.build": (), "cosets.union": ("syndrome_oracle",)}
