"""Composable cyclic ring codes, with explicitly scoped exact certificates.

Recipes synthesize generators from binary factors, never from supplied benchmark
coordinates. Large codes stay implicit. The only specialized proof dispatch is
to the two fixed, independently implemented, published DGLS parent verifiers.
"""
from __future__ import annotations

from copy import deepcopy
import importlib.util
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp

from discovery.core import Artifact, RunContext, Verification

REPO = Path(__file__).resolve().parents[2]
_SOURCE = REPO / "verifications/improved/dim88-96-z8-golay/source"
_DEFAULT_LIMITS = dict(max_ring_length=127, max_modulus=16,
                       max_enumeration_words=2**18, max_enumeration_cells=2**24,
                       max_gray_bits=8, max_exact_words=2**24)


def _integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def _limits(context):
    result = dict(_DEFAULT_LIMITS)
    if context is not None:
        for name in result:
            if name in context.budget:
                result[name] = min(result[name], _integer(context.budget[name], name))
        if "max_codewords" in context.budget:
            result["max_enumeration_words"] = min(result["max_enumeration_words"],
                _integer(context.budget["max_codewords"], "max_codewords"))
        if "max_dgls_words" in context.budget:
            result["max_exact_words"] = min(result["max_exact_words"],
                _integer(context.budget["max_dgls_words"], "max_dgls_words"))
    return result


def _params(params, allowed):
    unknown = set(params)-set(allowed)
    if unknown:
        raise ValueError(f"unknown parameters: {sorted(unknown)}")


def _verification(status, method, facts, **details):
    return Verification(status=status, method=method,
                        details={"facts": deepcopy(facts), **facts, **details})


def _single(inputs, kind):
    if len(inputs) != 1 or inputs[0].kind != kind:
        raise ValueError(f"expected one {kind} input")
    return inputs[0]


def binary_factors(n, context=None):
    """Monic, half-degree irreducible factors of 1+x+...+x^(n-1)."""
    _integer(n, "n", 3, _limits(context)["max_ring_length"])
    if n % 2 == 0:
        raise ValueError("n must be odd (coprime square-free factorization)")
    x = sp.Symbol("x")
    f = sp.Poly(sum(x**i for i in range(n)), x, modulus=2)
    _, factors = sp.factor_list(f)
    answer = []
    for g, multiplicity in factors:
        if multiplicity != 1:
            raise ValueError("binary factors must be square-free")
        if g.degree() == (n-1)//2:
            h = f.exquo(g)
            answer.append(([int(g.nth(i)) % 2 for i in range(g.degree()+1)],
                           [int(h.nth(i)) % 2 for i in range(h.degree()+1)]))
    return sorted(answer)


def binary_factor(inputs, params, context):
    _params(params, ("n", "index"))
    if inputs:
        raise ValueError("binary_factor takes no inputs")
    n = params.get("n", 7)
    index = _integer(params.get("index", 0), "index")
    factors = binary_factors(n, context)
    if index >= len(factors):
        raise ValueError("factor index does not exist")
    g, h = factors[index]
    facts = dict(cyclic_length=n, modulus=2, degree=len(g)-1)
    return Artifact(kind="polynomial_pair", data=dict(n=n, index=index, modulus=2,
                    g=g, h=h, lift_history=[]), facts=facts,
                    verification=_verification("exact", "binary polynomial factorization", facts))


def _solve2(A, b):
    matrix = np.column_stack([A % 2, b % 2]).astype(np.int64)
    row = 0
    pivots = []
    for col in range(matrix.shape[1]-1):
        ids = np.flatnonzero(matrix[row:, col])
        if not len(ids):
            continue
        pivot = row + int(ids[0])
        matrix[[row, pivot]] = matrix[[pivot, row]]
        for i in range(len(matrix)):
            if i != row and matrix[i, col]:
                matrix[i] ^= matrix[row]
        pivots.append(col)
        row += 1
    if row != matrix.shape[1]-1 or np.any((matrix[:, :-1].sum(axis=1) == 0) & (matrix[:, -1] != 0)):
        raise ValueError("factor lift is not unique")
    result = np.zeros(matrix.shape[1]-1, dtype=np.int64)
    for i, j in enumerate(pivots):
        result[j] = matrix[i, -1]
    return result


def _lift(g, h, modulus):
    g, h = np.array(g, dtype=np.int64), np.array(h, dtype=np.int64)
    target = np.ones(len(g)+len(h)-1, dtype=np.int64)
    a, b, current = len(g)-1, len(h)-1, 2
    history = []
    while current < modulus:
        residual = target - np.convolve(g, h)
        if np.any(residual % current):
            raise ValueError("invalid input factorization")
        matrix = np.zeros((a+b, a+b), dtype=np.int64)
        for i in range(a):
            matrix[i:i+b+1, i] = h
        for i in range(b):
            matrix[i:i+a+1, a+i] = g
        correction = _solve2(matrix, (residual//current)[:-1])
        g[:-1] += current*correction[:a]
        h[:-1] += current*correction[a:]
        current *= 2
        g %= current
        h %= current
        if not np.array_equal(np.convolve(g, h) % current, target):
            raise ValueError("Hensel lift failed its product identity")
        history.append(dict(modulus=current, g=g.tolist(), h=h.tolist()))
    return g.tolist(), h.tolist(), history


def _pair(data, context):
    base = binary_factor([], {"n": data["n"], "index": data["index"]}, context)
    modulus = _integer(data["modulus"], "modulus", 2, _limits(context)["max_modulus"])
    if modulus & (modulus-1):
        raise ValueError("modulus must be a power of two")
    g, h, history = _lift(base.data["g"], base.data["h"], modulus)
    expected = dict(n=data["n"], index=data["index"], modulus=modulus, g=g, h=h, lift_history=history)
    if data != expected:
        raise ValueError("polynomial recipe or factors were modified")
    return expected


def hensel(inputs, params, context):
    _params(params, ("modulus",))
    source = _single(inputs, "polynomial_pair")
    pair = _pair(source.data, context)
    modulus = _integer(params.get("modulus", 4), "modulus", 2, _limits(context)["max_modulus"])
    if modulus & (modulus-1) or modulus < pair["modulus"]:
        raise ValueError("lift modulus must be a power of two at least the source modulus")
    base = binary_factor([], {"n": pair["n"], "index": pair["index"]}, context)
    g, h, history = _lift(base.data["g"], base.data["h"], modulus)
    data = dict(n=pair["n"], index=pair["index"], modulus=modulus, g=g, h=h, lift_history=history)
    facts = dict(cyclic_length=pair["n"], modulus=modulus, degree=len(g)-1)
    return Artifact(kind="polynomial_pair", data=data, facts=facts,
                    verification=_verification("exact", "coprime Hensel product identities", facts))


def extend(inputs, params, context):
    _params(params, ("extension",))
    pair = _pair(_single(inputs, "polynomial_pair").data, context)
    q, n, g = pair["modulus"], pair["n"], pair["g"]
    extension = params.get("extension", -1)
    if type(extension) is not int:
        raise ValueError("extension must be an integer")
    extension %= q
    rank = n-len(g)+1
    generator = np.zeros((rank, n+1), dtype=np.int64)
    for i in range(rank):
        generator[i, i:i+len(g)] = g
    generator[:, -1] = extension*generator[:, :-1].sum(axis=1) % q
    if np.any(np.tril(generator[:, :rank], -1)) or np.any(np.diag(generator[:, :rank]) % 2 != 1):
        raise ValueError("generator has no triangular unit minor")
    facts = dict(length=n+1, rank=rank, modulus=q, cardinality=q**rank, min_distance=None)
    data = dict(polynomial_pair=pair, extension=extension, generator=generator.tolist())
    return Artifact(kind="ring_code", data=data, facts=facts,
                    verification=_verification("exact", "ring structure and free cardinality only", facts,
                                               scope="No distance certificate"))


def _ring(data, context):
    pair = _pair(data["polynomial_pair"], context)
    source = Artifact(kind="polynomial_pair", data=pair)
    rebuilt = extend([source], {"extension": data["extension"]}, context)
    if rebuilt.data != data:
        raise ValueError("ring generator or recipe was modified")
    return rebuilt


def gray_table(modulus, context=None):
    q = _integer(modulus, "modulus", 4, _limits(context)["max_modulus"])
    if q & (q-1):
        raise ValueError("Gray modulus must be a power of two")
    bits = q.bit_length()-1
    size = q//2
    if size > _limits(context)["max_gray_bits"]:
        raise ValueError("Gray alphabet exceeds local table budget")
    X = list(product((0, 1), repeat=bits-1))
    table = np.array([[(a >> (bits-1) & 1) ^ (sum(((a >> j) & 1)*x[j]
                        for j in range(bits-1)) % 2) for x in X] for a in range(q)], dtype=np.uint8)
    weight = np.full(q, q//4, dtype=np.int64)
    weight[0], weight[q//2] = 0, q//2
    for a, b in product(range(q), repeat=2):
        if np.count_nonzero(table[a] != table[b]) != weight[(a-b) % q]:
            raise ValueError("Gray map fails difference isometry")
    return table, weight


def _exhaustive(ring, weight, context):
    q = ring.facts["modulus"]
    G = np.array(ring.data["generator"], dtype=np.int64)
    count = ring.facts["cardinality"]
    limits = _limits(context)
    if count > limits["max_enumeration_words"] or count*G.shape[1] > limits["max_enumeration_cells"]:
        return None
    histogram = np.zeros(G.shape[1]*q//2+1, dtype=np.int64)
    for start in range(0, count, 4096):
        labels = np.arange(start, min(start+4096, count), dtype=np.int64)
        coefficients = np.column_stack([(labels//q**j) % q for j in range(len(G))])
        words = coefficients @ G % q
        histogram += np.bincount(weight[words].sum(axis=1), minlength=len(histogram))
    if int(histogram.sum()) != count or histogram[0] != 1:
        raise ValueError("enumeration failed free-module cardinality check")
    return dict(method="exhaustive homogeneous weight", histogram=histogram.tolist(),
                words_checked=count, min_distance=int(np.flatnonzero(histogram[1:])[0]+1))


def gray(inputs, params, context):
    _params(params, ())
    ring = _ring(_single(inputs, "ring_code").data, context)
    table, weight = gray_table(ring.facts["modulus"], context)
    proof = _exhaustive(ring, weight, context)
    facts = dict(length=ring.facts["length"]*len(table[0]), cardinality=ring.facts["cardinality"],
                 min_distance=proof["min_distance"] if proof else None)
    proof = proof or dict(method="structural candidate", reason="exhaustive budget exceeded")
    return Artifact(kind="binary_code", data=dict(construction="gray", ring_code=ring.data,
                    gray_table=table.tolist(), proof=proof), facts=facts,
                    verification=_verification("exact" if facts["min_distance"] else "unverified",
                      proof["method"], facts, implicit=True))


def _load_parent(filename):
    # Fixed allowlist; artifact payloads never select executable paths.
    if filename not in ("replay.py", "verify_independent.py"):
        raise ValueError("unknown exact verifier")
    path = _SOURCE / filename
    spec = importlib.util.spec_from_file_location("_discovery_dgls_"+path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_dgls(inputs, params, context):
    _params(params, ())
    if not __debug__:
        raise ValueError("DGLS parent verifiers require Python without -O")
    ring = _ring(_single(inputs, "ring_code").data, context)
    pair = ring.data["polynomial_pair"]
    if (pair["n"], pair["modulus"], ring.data["extension"]) != (23, 8, 7):
        raise ValueError("DGLS proof requires the extended length-23 Z8 recipe")
    if _limits(context)["max_exact_words"] < 2**24:
        raise ValueError("DGLS proof requires budget for 2^24 Z4 proof words")
    table, _ = gray_table(8, context)
    raw = dict(generator=ring.data["generator"], generator_polynomial_low_first=pair["g"],
               complementary_factor_low_first=pair["h"])
    try:
        _, first = _load_parent("replay.py").verify_z8(raw)
        _, second = _load_parent("verify_independent.py").parent_proof(raw)
    except AssertionError as exc:
        raise ValueError("DGLS exact parent verification failed") from exc
    if first != second or first["grey_map"] != table.tolist():
        raise ValueError("independent DGLS proof receipts disagree")
    facts = dict(length=96, cardinality=2**36, min_distance=24)
    proof = dict(method="two independent DGLS parent proofs", first=first, second=second,
                 scope="published DGLS parent, before translate union")
    return Artifact(kind="binary_code", data=dict(construction="gray", ring_code=ring.data,
                    gray_table=table.tolist(), proof=proof), facts=facts,
                    verification=_verification("exact", proof["method"], facts, implicit=True))


def _rebuild_binary(data, context):
    construction = data.get("construction")
    if construction == "gray":
        ring = _ring(data["ring_code"], context)
        method = data["proof"]["method"]
        if method == "two independent DGLS parent proofs":
            return verify_dgls([ring], {}, context)
        if method in ("structural candidate", "exhaustive homogeneous weight"):
            return gray([ring], {}, context)
        raise ValueError("unknown binary distance proof")
    if construction == "translate_union":
        if not isinstance(data.get("parent"), dict) or data["parent"].get("construction") != "gray":
            raise ValueError("union parent must be a single Gray recipe")
        parent = _rebuild_binary(data["parent"], context)
        if parent.data != data["parent"]:
            raise ValueError("union parent certificate was modified")
        return _translate(parent, data["target_distance"], context)
    raise ValueError("unknown binary construction")


def _translate(parent, target, context):
    if parent.data.get("construction") != "gray":
        raise ValueError("repeated-mask union currently requires a Gray parent")
    distance = parent.facts["min_distance"]
    if parent.verification.status != "exact" or distance is None:
        raise ValueError("translate union requires a verified parent distance")
    target = _integer(target, "target_distance", 1, distance)
    ring = _ring(parent.data["ring_code"], context)
    table, _ = gray_table(ring.facts["modulus"], context)
    L = len(table[0])
    if L > 8:
        raise ValueError("local mask search exceeds 8-bit limit")
    representatives, seen = [], set()
    for mask_tuple in product((0, 1), repeat=L):
        mask = np.array(mask_tuple, dtype=np.uint8)
        coset = frozenset(tuple(map(int, row ^ mask)) for row in table)
        if coset not in seen:
            seen.add(coset)
            representatives.append(mask)
    selected, pair_bounds = [], []
    for mask in representatives:
        bounds = [ring.facts["length"]*min(int(np.count_nonzero(x ^ mask ^ y ^ old))
                  for x in table for y in table) for old in selected]
        if all(bound >= target for bound in bounds):
            selected.append(mask)
            pair_bounds.extend(bounds)
    cross = min(pair_bounds) if pair_bounds else None
    lower = min(distance, cross) if cross is not None else distance
    facts = dict(length=parent.facts["length"], cardinality=parent.facts["cardinality"]*len(selected),
                 min_distance=lower)
    certificate = dict(method="coordinatewise repeated-mask lower bounds", local_cosets=len(representatives),
                       selected_repeated_masks=[m.tolist() for m in selected], multiplier=len(selected),
                       cross_distance_lower_bound=cross,
                       scope="greedy restricted translates; no optimality claim")
    return Artifact(kind="binary_code", data=dict(construction="translate_union", parent=parent.data,
                    target_distance=target, proof=certificate), facts=facts,
                    verification=_verification("exact", certificate["method"], facts, implicit=True))


def translate_union(inputs, params, context):
    _params(params, ("target_distance",))
    original = _single(inputs, "binary_code")
    parent = _rebuild_binary(original.data, context)
    if parent.data != original.data or parent.facts != original.facts:
        raise ValueError("parent binary certificate or facts were modified")
    target = params.get("target_distance", parent.facts["min_distance"])
    return _translate(parent, target, context)


def verify_artifact(artifact, context=None):
    """Replay embedded recipes and proofs; never trust serialized exact status."""
    try:
        if artifact.kind == "polynomial_pair":
            data = _pair(artifact.data, context)
            rebuilt = hensel([binary_factor([], {"n": data["n"], "index": data["index"]}, context)],
                             {"modulus": data["modulus"]}, context)
        elif artifact.kind == "ring_code":
            rebuilt = _ring(artifact.data, context)
        elif artifact.kind == "binary_code":
            rebuilt = _rebuild_binary(artifact.data, context)
        else:
            raise ValueError("not a ring engine artifact")
        if rebuilt.data != artifact.data or rebuilt.facts != artifact.facts:
            raise ValueError("artifact data, facts, or proof were modified")
        return rebuilt.verification
    except (KeyError, TypeError, IndexError, OverflowError) as exc:
        raise ValueError("malformed ring artifact") from exc


OPERATORS = {"rings.binary_factor": binary_factor, "rings.hensel": hensel,
             "rings.extend": extend, "rings.gray": gray,
             "rings.verify_dgls": verify_dgls, "codes.translate_union": translate_union}
INPUT_KINDS = {"rings.binary_factor": (), "rings.hensel": ("polynomial_pair",),
               "rings.extend": ("polynomial_pair",), "rings.gray": ("ring_code",),
               "rings.verify_dgls": ("ring_code",), "codes.translate_union": ("binary_code",)}
