"""Bounded finite unitary orbit discovery and an independently checked control.

All floating point orbit closure and optimization results remain numerical.
Projective gauge is used only to identify lines: phase blocks are always built
from closed *vector* orbits under the original generators.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import isfinite, lcm
from hashlib import sha256

import numpy as np

from discovery.core import Artifact, RunContext, Verification


def _integer(params, name, default, minimum=0, maximum=100000):
    value = params.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def _limit(params, context, name, default, maximum=100000):
    value = _integer(params, name, default, 1, maximum)
    budget = getattr(context, "budget", {}) or {}
    if name in budget:
        value = min(value, _integer(budget, name, value, 1, maximum))
    return value


def _known_params(params, allowed):
    if not isinstance(params, dict):
        raise ValueError("group operator parameters must be an object")
    unknown = set(params) - set(allowed)
    if unknown:
        raise ValueError("unknown group operator parameters: " + ", ".join(sorted(unknown)))


def canonical_line(vector):
    vector = np.asarray(vector, dtype=complex)
    norm = np.linalg.norm(vector)
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("seed must be a finite nonzero vector")
    vector = vector / norm
    pivot = np.flatnonzero(abs(vector) > 1e-9)[0]
    return vector * np.conj(vector[pivot]) / abs(vector[pivot])


def _key(vector):
    return tuple(np.round(np.r_[vector.real, vector.imag], 7))


def orbit(seed, generators, limit=500, projective=False):
    """Return (partial vectors, status); never pass a truncated orbit as closed."""
    if limit < 1:
        raise ValueError("orbit limit must be positive")
    normalize = canonical_line if projective else lambda v: v / np.linalg.norm(v)
    initial = normalize(np.asarray(seed, dtype=complex))
    vectors, seen = [initial], {_key(initial)}
    for vector in vectors:
        for generator in generators:
            candidate = normalize(generator @ vector)
            key = _key(candidate)
            if key not in seen:
                if len(vectors) == limit:
                    return np.asarray(vectors), "limit_reached"
                seen.add(key)
                vectors.append(candidate)
    return np.asarray(vectors), "closed_numerically"


def _generator_data(real, imag, context=None):
    """Check JSON shape and allocation ceilings before creating dense arrays."""
    if not isinstance(real, list) or not real or not isinstance(real[0], list) or not real[0]:
        raise ValueError("generators_real must be a nonempty list of square matrices")
    dimension, count = len(real[0]), len(real)
    budget = getattr(context, "budget", {}) or {}
    if dimension > budget.get("max_dimension", 256):
        raise ValueError("representation dimension exceeds budget")
    # Count both scalar components even when imaginary zeros are implicit.
    if 2 * count * dimension * dimension > budget.get("max_materialized_points", 20000):
        raise ValueError("generator scalar entries exceed materialization budget")
    if imag is not None and (not isinstance(imag, list) or len(imag) != count):
        raise ValueError("generators_imag must match generators_real")
    for component in ([real] if imag is None else [real, imag]):
        for matrix in component:
            if not isinstance(matrix, list) or len(matrix) != dimension:
                raise ValueError("generators must be equal square matrices")
            for row in matrix:
                if not isinstance(row, list) or len(row) != dimension:
                    raise ValueError("generators must be equal square matrices")
                for value in row:
                    try:
                        finite = type(value) in (int, float) and isfinite(value)
                    except OverflowError:
                        finite = False
                    if not finite:
                        raise ValueError("generator entries must be finite real numbers")
    real_array = np.asarray(real, dtype=float)
    imag_array = np.zeros_like(real_array) if imag is None else np.asarray(imag, dtype=float)
    return real_array + 1j * imag_array


def _generators(artifact, context=None):
    if artifact.kind != "unitary_representation":
        raise ValueError("groups.search requires a unitary_representation")
    data = artifact.data
    generators = _generator_data(data["generators_real"], data.get("generators_imag"), context)
    if not np.isfinite(generators).all():
        raise ValueError("generators must be finite")
    if any(np.max(abs(g.conj().T @ g - np.eye(g.shape[0]))) > 1e-9 for g in generators):
        raise ValueError("generators must be unitary to numerical tolerance")
    return generators


def from_generators(inputs, params, context):
    """Accept arbitrary bounded JSON unitary generators without a finiteness claim."""
    if inputs:
        raise ValueError("groups.from_generators takes no inputs")
    _known_params(params, {"generators_real", "generators_imag", "label"})
    if "generators_real" not in params:
        raise ValueError("groups.from_generators requires generators_real")
    label = params.get("label", "supplied unitary generators")
    if not isinstance(label, str) or not label or len(label) > 200:
        raise ValueError("label must be a nonempty string of at most 200 characters")
    generators = _generator_data(params["generators_real"], params.get("generators_imag"), context)
    artifact = Artifact(kind="unitary_representation", data={"generators_real": generators.real.tolist(),
                        "generators_imag": generators.imag.tolist(), "family": "supplied", "label": label},
                        facts={"complex_dimension": generators.shape[1], "generator_count": len(generators)})
    _generators(artifact, context)
    artifact.verification = verify_artifact(artifact, context)
    return artifact


def weil(inputs, params, context):
    """Fourier/quadratic-phase generators on an inversion parity subspace."""
    if inputs:
        raise ValueError("groups.weil takes no inputs")
    _known_params(params, {"p", "rank", "parity", "max_ambient_dimension"})
    p = _integer(params, "p", 3, 3, 31)
    rank = _integer(params, "rank", 2, 1, 5)
    parity = params.get("parity", 1)
    if parity not in (-1, 1) or isinstance(parity, bool):
        raise ValueError("parity must be -1 or 1")
    if any(p % k == 0 for k in range(2, int(p**0.5) + 1)):
        raise ValueError("p must be an odd prime")
    ambient = p**rank
    if ambient > _limit(params, context, "max_ambient_dimension", 243, 2000):
        raise ValueError("Weil ambient dimension exceeds materialization budget")
    budget = getattr(context, "budget", {}) or {}
    if (ambient + parity) // 2 > budget.get("max_dimension", 256):
        raise ValueError("Weil representation dimension exceeds budget")
    dimension = (ambient + parity) // 2
    generator_count = 2 * rank + rank * (rank - 1) // 2
    if 2 * generator_count * dimension * dimension > budget.get("max_materialized_points", 20000):
        raise ValueError("Weil generator scalar entries exceed materialization budget")
    points = list(product(range(p), repeat=rank))
    index = {point: i for i, point in enumerate(points)}
    basis, seen = [], set()
    for point in points:
        if point in seen:
            continue
        opposite = tuple(-x % p for x in point)
        seen.update((point, opposite))
        vector = np.zeros(ambient)
        if point == opposite:
            if parity < 0:
                continue
            vector[index[point]] = 1
        else:
            vector[index[point]], vector[index[opposite]] = 1 / np.sqrt(2), parity / np.sqrt(2)
        basis.append(vector)
    basis = np.asarray(basis).T
    root = np.exp(2j * np.pi / p)
    generators = []
    for axis in range(rank):
        fourier = np.zeros((ambient, ambient), complex)
        for point in points:
            for value in range(p):
                target = list(point)
                target[axis] = value
                fourier[index[tuple(target)], index[point]] = root**(value * point[axis]) / np.sqrt(p)
        generators.append(basis.T @ (fourier * (-1j if p % 4 == 3 else 1)) @ basis)
        diagonal = np.asarray([root**(x[axis] ** 2) for x in points])
        generators.append(basis.T @ (diagonal[:, None] * basis))
    for a in range(rank):
        for b in range(a):
            diagonal = np.asarray([root**(x[a] * x[b]) for x in points])
            generators.append(basis.T @ (diagonal[:, None] * basis))
    generators = np.asarray(generators)
    artifact = Artifact(kind="unitary_representation", data={
        "generators_real": generators.real.tolist(), "generators_imag": generators.imag.tolist(),
        "family": "weil", "p": p, "rank": rank, "parity": parity,
    }, facts={"complex_dimension": len(basis.T), "generator_count": len(generators)})
    artifact.verification = verify_artifact(artifact, context)
    return artifact


def max_cosine(left, right=None):
    same = right is None
    if same:
        right = left
    best = -1.0
    for start in range(0, len(left), 128):
        gram = (left[start:start + 128] @ right.conj().T).real
        if same:
            rows = np.arange(len(gram))
            gram[rows, start + rows] = -np.inf
        if gram.size:
            best = max(best, float(gram.max()))
    return best


def phase_block(vector_orbit, order, offset=0.0):
    """Return full block and one representative per constituent group orbit."""
    phases = np.exp(1j * (offset + 2 * np.pi * np.arange(order) / order))
    full = np.concatenate([phase * vector_orbit for phase in phases])
    _, indices = np.unique(np.round(np.c_[full.real, full.imag], 7), axis=0, return_index=True)
    return full[np.sort(indices)], phases[:, None] * vector_orbit[0]


def _block_fingerprint(full):
    """Order-independent numerical identity at the same seven-digit screening precision."""
    quantized = np.rint(np.c_[full.real, full.imag] * 10**7).astype(np.int64)
    order = np.lexsort(quantized[:, ::-1].T)
    return (quantized.shape, sha256(quantized[order].astype("<i8", copy=False).tobytes()).hexdigest())


def _bounded_phase_block(vector_orbit, order, limit):
    """Deduplicate while constructing; a raw m*N estimate can greatly overcount."""
    if len(vector_orbit) > limit:
        return None
    if len(vector_orbit) * order <= limit:
        return phase_block(vector_orbit, order)
    phases = np.exp(2j * np.pi * np.arange(order) / order)
    rows, seen = [], set()
    for phase in phases:
        for vector in vector_orbit * phase:
            key = _key(vector)
            if key in seen:
                continue
            if len(rows) == limit:
                return None
            seen.add(key)
            rows.append(vector)
    return np.asarray(rows), phases[:, None] * vector_orbit[0]


def compressed_max_cosine(representatives, full, exclude_identical=False):
    gram = (representatives @ full.conj().T).real
    if exclude_identical:
        # Only self copies are removed. Distinct near coincident vectors must
        # remain conflicts; the same seven-decimal keys define deduplication.
        full_keys = [_key(v) for v in full]
        for i, rep in enumerate(representatives):
            key = _key(rep)
            for j, other in enumerate(full_keys):
                if key == other:
                    gram[i, j] = -np.inf
    return max(-1.0, float(gram.max())) if gram.size else -1.0


def _sample_orbits(generators, trials, linecap, seed, max_candidates, max_orbits):
    rng = np.random.default_rng(seed)
    dimension = generators.shape[1]
    candidates = [("coordinate:" + str(i), vector) for i, vector in enumerate(np.eye(dimension, dtype=complex))][:max_candidates]
    words = []
    candidate_limit = dimension > max_candidates
    for trial in range(trials):
        word = rng.integers(len(generators), size=int(rng.integers(3, 18))).tolist()
        matrix = np.eye(dimension, dtype=complex)
        for generator in word:
            matrix = generators[generator] @ matrix
        values, vectors = np.linalg.eig(matrix)
        for j, value in enumerate(values):
            if np.count_nonzero(abs(values - value) < 1e-6) == 1:
                candidates.append((f"cyclic:{trial}:{j}", vectors[:, j]))
        if trial and trial % 3 == 0 and dimension > 1:
            other = words[int(rng.integers(len(words)))]
            other_values = np.linalg.eigvals(other)
            for value in np.unique(np.round(values, 8)):
                for second in np.unique(np.round(other_values, 8)):
                    _, singular, vh = np.linalg.svd(np.vstack((matrix - value * np.eye(dimension), other - second * np.eye(dimension))))
                    if singular[-1] < 1e-6 and singular[-2] > 1e-5:
                        candidates.append((f"joint:{trial}", vh[-1].conj()))
        words.append(matrix)
        if len(candidates) >= max_candidates:
            candidates = candidates[:max_candidates]
            candidate_limit = True
            break
    seen, orbits, diagnostics = set(), [], []
    for label, vector in candidates:
        if _key(canonical_line(vector)) in seen:
            continue
        lines, status = orbit(vector, generators, linecap, True)
        diagnostics.append({"seed": label, "lines": len(lines), "status": status})
        if status == "closed_numerically":
            seen.update(_key(v) for v in lines)
            orbits.append(lines)
            if len(orbits) == max_orbits:
                return orbits, diagnostics, True
    return orbits, diagnostics, candidate_limit


def search(inputs, params, context):
    if len(inputs) != 1:
        raise ValueError("groups.search takes one representation")
    _known_params(params, {"trials", "seed", "linecap", "vectorcap", "max_materialized_points",
                          "max_candidates", "max_orbits", "max_blocks", "max_output_points",
                          "phase_steps", "phase_orders", "solver_seconds"})
    generators = _generators(inputs[0], context)
    trials = _integer(params, "trials", 40, 0, 10000)
    budget = getattr(context, "budget", {}) or {}
    if generators.shape[1] > budget.get("max_dimension", 256):
        raise ValueError("representation dimension exceeds budget")
    trials = min(trials, _integer(budget, "max_trials", trials, 0, 10000))
    seed = _integer(params, "seed", 20261104, 0, 2**63 - 1)
    linecap = _limit(params, context, "linecap", 500)
    vectorcap = _limit(params, context, "vectorcap", 6 * linecap)
    if "max_orbit_points" in budget:
        cap = _integer(budget, "max_orbit_points", vectorcap, 1, 1000000)
        linecap, vectorcap = min(linecap, cap), min(vectorcap, cap)
    max_points = _limit(params, context, "max_materialized_points", 100000, 1000000)
    max_candidates = _limit(params, context, "max_candidates", 2000)
    max_orbits = _limit(params, context, "max_orbits", 100)
    max_blocks = _limit(params, context, "max_blocks", 1000)
    max_output = _limit(params, context, "max_output_points", 10000, 1000000)
    phase_steps = _integer(params, "phase_steps", 12, 1, 120)
    orders = params.get("phase_orders", sorted({1, 2, 3, 4, 6, 2 * inputs[0].data.get("p", 3)}))
    if not isinstance(orders, list) or not orders or any(isinstance(m, bool) or not isinstance(m, int) or not 1 <= m <= 60 for m in orders):
        raise ValueError("phase_orders must contain integers in [1, 60]")
    orbits, diagnostics, sampling_limit = _sample_orbits(generators, trials, linecap, seed, max_candidates, max_orbits)
    blocks, representatives, inventory = [], [], []
    materialized, materialization_limit = 0, False
    block_fingerprints, duplicates_skipped = set(), 0
    for i, lines in enumerate(orbits):
        vectors, status = orbit(lines[0], generators, vectorcap)
        entry = {"orbit": i, "lines": len(lines), "vectors": len(vectors), "vector_status": status, "accepted_phase_orders": [],
                 "admitted_blocks": 0, "duplicate_blocks_skipped": 0, "capacity_blocks_skipped": 0}
        inventory.append(entry)
        if status != "closed_numerically":
            continue
        if compressed_max_cosine(vectors[:1], vectors, True) > 0.5 + 1e-8:
            continue
        prepared = []
        for order in sorted(set(orders)):
            bounded = _bounded_phase_block(vectors, order, max_points)
            if bounded is None:
                materialization_limit = True
                continue
            full, reps = bounded
            if compressed_max_cosine(reps, full, True) > 0.5 + 1e-8:
                continue
            entry["accepted_phase_orders"].append(order)
            prepared.append((order, full, reps))
        # Give each phase order its first candidate before filling one order's
        # entire grid. Antipodal vector orbits often make orders 1/2 or 3/6
        # identical, and invariant phase shifts also repeat whole blocks.
        for step in range(phase_steps):
            for order, full, reps in prepared:
                phase = np.exp(2j * np.pi * step / (2 * phase_steps))
                candidate = full * phase
                fingerprint = _block_fingerprint(candidate)
                if fingerprint in block_fingerprints:
                    duplicates_skipped += 1
                    entry["duplicate_blocks_skipped"] += 1
                    continue
                if len(blocks) >= max_blocks or materialized + len(full) > max_points:
                    materialization_limit = True
                    entry["capacity_blocks_skipped"] += 1
                    continue
                block_fingerprints.add(fingerprint)
                blocks.append((f"o{i}:m{order}:phase{step}/{2 * phase_steps}", candidate))
                representatives.append(reps * phase)
                materialized += len(full)
                entry["admitted_blocks"] += 1
    selection_params = dict(params)
    seconds = params.get("solver_seconds", 10)
    if isinstance(seconds, bool) or not isinstance(seconds, (int, float)) or not np.isfinite(seconds) or not 0 < seconds <= 3600:
        raise ValueError("solver_seconds must be in (0, 3600]")
    if "solver_seconds" in budget:
        selection_params["solver_seconds"] = min(seconds, budget["solver_seconds"])
    chosen, solver = _select_blocks(blocks, representatives, selection_params)
    count = sum(len(blocks[i][1]) for i in chosen)
    facts = {"cardinality": count, "complex_dimension": generators.shape[1], "seed": seed,
             "search_scope": "sampled subgroup eigenlines and finite phase grid only",
             "global_optimality_claim": False}
    data = {"inventory": inventory, "seed_diagnostics": diagnostics, "solver": solver,
            "selected_blocks": [blocks[i][0] for i in chosen],
            "limits": {"trials": trials, "linecap": linecap, "vectorcap": vectorcap,
                       "max_candidates": max_candidates, "max_orbits": max_orbits,
                       "max_blocks": max_blocks, "max_materialized_points": max_points,
                       "max_output_points": max_output},
            "sampling_limit_reached": sampling_limit, "materialization_limit_reached": materialization_limit,
            "block_coverage": {"unique_blocks_admitted": len(blocks), "admitted_point_rows": materialized,
                               "duplicate_blocks_skipped": duplicates_skipped,
                               "identity_precision_decimal_digits": 7,
                               "schedule": "phase step then phase order within each closed orbit; orbit discovery order unchanged"},
            "closure_failures": sum(item["status"] != "closed_numerically" for item in diagnostics),
            "cross_checks": "one representative per constituent vector orbit against full other block"}
    within_pairs = count * (count - 1) // 2 <= budget.get("max_pairs", 10000000)
    if count and count <= max_output and within_pairs:
        full = np.concatenate([blocks[i][1] for i in chosen])
        data["real_coordinates"] = np.c_[full.real, full.imag].tolist()
        artifact = Artifact(kind="spherical_code", data=data, facts=facts)
        artifact.verification = verify_artifact(artifact, context)
        if artifact.verification.status != "numerical":
            raise ValueError("selected group blocks failed full numerical verification")
        artifact.facts["dimension"] = artifact.verification.details["dimension"]
    else:
        data["output_status"] = "empty" if count == 0 else ("output_limit_reached" if count > max_output else "verification_pair_limit_reached")
        artifact = Artifact(kind="search_result", data=data, facts=facts,
                            verification=Verification(status="unverified", method="bounded unitary orbit screening", details={"reason": data["output_status"]}))
    return artifact


def _select_blocks(blocks, representatives, params):
    if not blocks:
        return [], {"status": "no_candidates", "candidate_count": 0}
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    edges = [(i, j) for i in range(len(blocks)) for j in range(i)
             if compressed_max_cosine(representatives[i], blocks[j][1]) > 0.5 + 1e-8]
    matrix = lil_matrix((len(edges), len(blocks)))
    for row, (i, j) in enumerate(edges):
        matrix[row, i] = matrix[row, j] = 1
    seconds = params.get("solver_seconds", 10)
    if isinstance(seconds, bool) or not isinstance(seconds, (int, float)) or not np.isfinite(seconds) or not 0 < seconds <= 3600:
        raise ValueError("solver_seconds must be in (0, 3600]")
    result = milp(-np.asarray([len(block) for _, block in blocks], float),
                  integrality=np.ones(len(blocks)), bounds=Bounds(0, 1),
                  constraints=LinearConstraint(matrix.tocsr(), -np.inf, 1),
                  options={"time_limit": float(seconds)})
    chosen = [] if result.x is None else np.flatnonzero(result.x > 0.5).tolist()
    chosen_set = set(chosen)
    if any(i in chosen_set and j in chosen_set for i, j in edges):
        raise ValueError("solver returned an infeasible incumbent")
    bound = getattr(result, "mip_dual_bound", None)
    return chosen, {"status": int(result.status), "message": str(result.message),
                    "candidate_count": len(blocks), "conflicts": len(edges),
                    "finite_candidate_upper_bound": None if bound is None or not np.isfinite(bound) else float(-bound),
                    "bound_scope": "materialized candidate conflict graph; numerical screening"}


def _build10():
    """Generate the 510 point control using rational Q(sqrt(-3)) arithmetic."""
    F = Fraction
    zero, one = (F(0), F(0)), (F(1), F(0))
    roots = [one, (F(-1, 2), F(1, 2)), (F(-1, 2), F(-1, 2))]
    points = list(product(range(3), repeat=2))
    index = {point: i for i, point in enumerate(points)}

    def add(x, y):
        return x[0] + y[0], x[1] + y[1]

    def mul(x, y):
        return x[0] * y[0] - 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]

    def transform(vector, generator):
        if generator in (0, 2):
            axis = generator // 2
            answer = []
            for target in points:
                value = zero
                for coordinate in range(3):
                    source = list(target)
                    source[axis] = coordinate
                    value = add(value, mul(roots[target[axis] * coordinate % 3], vector[index[tuple(source)]]))
                answer.append(mul((F(0), F(-1, 3)), value))
            return tuple(answer)
        axis = (generator - 1) // 2
        return tuple(mul(roots[(point[axis] ** 2 if generator < 4 else point[0] * point[1]) % 3], vector[i])
                     for i, point in enumerate(points))

    def exact_orbit(seed):
        vectors, seen = [tuple(seed)], {tuple(seed)}
        for vector in vectors:
            for generator in range(5):
                candidate = transform(vector, generator)
                if candidate not in seen:
                    if len(vectors) == 1000:
                        raise ValueError("exact control orbit exceeded construction ceiling")
                    seen.add(candidate)
                    vectors.append(candidate)
        return vectors

    first, second = [zero] * 9, [zero] * 9
    first[0] = one
    for i in (3, 6):
        second[i] = (F(1, 2), F(0))
    for i in (4, 8):
        second[i] = (F(1, 4), F(1, 4))
    a, b = exact_orbit(first), exact_orbit(second)
    vectors = [tuple(mul(root, x) for x in vector) for root in roots for vector in a] + b
    denominator = lcm(*(value.denominator for vector in vectors for pair in vector for value in pair))
    coordinates = [[int(pair[0] * denominator) for pair in vector] + [int(pair[1] * denominator) for pair in vector] for vector in vectors]
    return coordinates, [1] * 9 + [3] * 9, denominator**2, [len(a), len(b)]


def ganzhinov10(inputs, params, context):
    if inputs or params:
        raise ValueError("groups.ganzhinov10 takes no inputs or parameters")
    coordinates, weights, norm, sizes = _build10()
    artifact = Artifact(kind="spherical_code", data={"coordinates": coordinates, "weights": weights,
                        "norm_squared": norm, "construction": "Q(sqrt(-3)) Fourier/phase orbits; 3*80 + 270",
                        "orbit_sizes": sizes}, facts={"cardinality": 510, "dimension": 10,
                        "role": "known construction control", "global_optimality_claim": False})
    artifact.verification = verify_artifact(artifact, context)
    if artifact.verification.status != "exact":
        raise ValueError("exact group control failed verification")
    return artifact


def _build14():
    """Replay the E7 orthogonal-frame recipe without reading a coordinate file."""
    roots = []
    for axis in range(7):
        vector = [0] * 7
        vector[axis] = 2
        roots.append(tuple(vector))
    for line in range(7):
        support = sorted(set(range(7)) - {line, (line + 1) % 7, (line + 3) % 7})
        for signs in product((-1, 1), repeat=3):
            vector = [0] * 7
            for axis, sign in zip(support, (1, *signs)):
                vector[axis] = sign
            roots.append(tuple(vector))
    adjacency = [{j for j, second in enumerate(roots) if sum(x * y for x, y in zip(first, second)) == 0}
                 for first in roots]
    frames = []

    def extend(frame, candidates):
        if len(frame) == 7:
            frames.append(tuple(frame))
            return
        remaining = set(candidates)
        while len(frame) + len(remaining) >= 7:
            vertex = min(remaining)
            remaining.remove(vertex)
            extend(frame + [vertex], remaining & adjacency[vertex])

    extend([], set(range(len(roots))))
    containing = [[frame for frame in frames if vertex in frame] for vertex in range(len(roots))]

    def partition(remaining):
        if not remaining:
            return []
        pivot = min(remaining, key=lambda vertex: sum(set(frame) <= remaining for frame in containing[vertex]))
        for frame in containing[pivot]:
            if set(frame) <= remaining:
                rest = partition(remaining - set(frame))
                if rest is not None:
                    return [frame] + rest
        return None

    selected = partition(set(range(len(roots))))
    if selected is None or len(selected) != 9:
        raise ValueError("E7 control frame partition failed")
    coordinates = []
    for frame in selected:
        for first in frame:
            for second in frame:
                if first == second:
                    continue
                for left, right in product((-1, 1), repeat=2):
                    coordinates.append([left * x for x in roots[first]] + [right * x for x in roots[second]])
    for root in roots:
        for left, right in product((-1, 1), repeat=2):
            coordinates.append([left * x for x in root] + [right * x for x in root])
    for side in (0, 7):
        for first, second in combinations(range(7), 2):
            for left, right in product((-1, 1), repeat=2):
                vector = [0] * 14
                vector[side + first], vector[side + second] = 2 * left, 2 * right
                coordinates.append(vector)
    return coordinates, len(frames)


def ganzhinov14(inputs, params, context):
    if inputs or params:
        raise ValueError("groups.ganzhinov14 takes no inputs or parameters")
    coordinates, frame_count = _build14()
    artifact = Artifact(kind="spherical_code", data={"coordinates": coordinates, "weights": [1] * 14,
                        "norm_squared": 8, "construction": "E7 orthogonal-frame construction",
                        "frame_count": 9, "available_orthogonal_frames": frame_count,
                        "components": [1512, 252, 168],
                        "provenance": "known construction recipe replay; not autonomous U(3,3) subgroup discovery"},
                        facts={"cardinality": 1932, "dimension": 14, "role": "known construction control",
                               "global_optimality_claim": False})
    # The shared checker recomputes distinctness, all norms, all pairs and exact
    # rank from the freshly generated integers, with proven overflow protection.
    from discovery.certification import verify_integer_spherical
    artifact.verification = verify_integer_spherical(artifact.data, context)
    if artifact.verification.details["cardinality"] != 1932 or artifact.verification.details["dimension"] != 14:
        raise ValueError("E7 control cardinality or rank mismatch")
    return artifact


def verify_artifact(artifact, context=None):
    """Recheck coordinates rather than trusting stored facts or status."""
    try:
        if artifact.kind == "unitary_representation":
            generators = _generators(artifact, context)
            if "complex_dimension" in artifact.facts and artifact.facts["complex_dimension"] != generators.shape[1]:
                raise ValueError("claimed complex dimension mismatch")
            if "generator_count" in artifact.facts and artifact.facts["generator_count"] != len(generators):
                raise ValueError("claimed generator count mismatch")
            return Verification(status="numerical", method="unitarity residual", details={
                "max_residual": float(max(np.max(abs(g.conj().T @ g - np.eye(g.shape[0]))) for g in generators)),
                "complex_dimension": generators.shape[1], "generator_count": len(generators),
                "finite_group_proved": False})
        if artifact.kind != "spherical_code":
            raise ValueError("unsupported artifact kind")
        data = artifact.data
        rows = data.get("coordinates", data.get("real_coordinates", []))
        budget = getattr(context, "budget", {}) or {}
        if len(rows) > budget.get("max_materialized_points", 100000):
            raise ValueError("verification materialization budget exceeded")
        if len(rows) * (len(rows) - 1) // 2 > budget.get("max_pairs", 10000000):
            raise ValueError("verification pair budget exceeded")
        if "coordinates" in data:
            coordinates = data["coordinates"]
            weights, norm = data["weights"], data["norm_squared"]
            if not coordinates or not weights or isinstance(norm, bool) or not isinstance(norm, int) or norm <= 0:
                raise ValueError("invalid exact code schema")
            if any(isinstance(w, bool) or not isinstance(w, int) or w <= 0 for w in weights):
                raise ValueError("weights must be positive integers")
            if any(len(row) != len(weights) or any(isinstance(x, bool) or not isinstance(x, int) for x in row) for row in coordinates):
                raise ValueError("coordinates must be integer rows matching weights")
            if len(weights) > budget.get("max_dimension", 256):
                raise ValueError("exact coordinate dimension exceeds budget")
            if len(set(map(tuple, coordinates))) != len(coordinates):
                raise ValueError("duplicate coordinates")
            if any(sum(w * x * x for w, x in zip(weights, row)) != norm for row in coordinates):
                raise ValueError("incorrect exact norm")
            maximum = None
            for i, row in enumerate(coordinates):
                for other in coordinates[:i]:
                    inner = sum(w * x * y for w, x, y in zip(weights, row, other))
                    if inner * 2 > norm:
                        raise ValueError("pair exceeds kissing threshold")
                    maximum = inner if maximum is None else max(maximum, inner)
            import sympy
            rank = int(sympy.Matrix(coordinates).rank())
            if "cardinality" in artifact.facts and artifact.facts["cardinality"] != len(coordinates):
                raise ValueError("claimed cardinality mismatch")
            if "dimension" in artifact.facts and artifact.facts["dimension"] != rank:
                raise ValueError("claimed dimension mismatch")
            for key, value in {"ambient_dimension": len(weights), "norm_squared": norm, "max_inner_product": maximum}.items():
                if key in artifact.facts and artifact.facts[key] != value:
                    raise ValueError("claimed " + key + " mismatch")
            return Verification(status="exact", method="integer weighted Gram and rational rank", details={
                "cardinality": len(coordinates), "dimension": rank, "rank_exact": rank, "ambient_dimension": len(weights),
                "norm_squared": norm, "max_inner_product": maximum,
                "pairs_checked": len(coordinates) * (len(coordinates) - 1) // 2})
        # Validate the JSON shape and width before converting to an ndarray.
        # The materialization guard above bounds rows, but not column allocation.
        import math
        raw = data["real_coordinates"]
        if not isinstance(raw, list) or not raw or not isinstance(raw[0], list) or not raw[0]:
            raise ValueError("invalid numerical coordinate rows")
        width = len(raw[0])
        if width > budget.get("max_dimension", 256):
            raise ValueError("numerical coordinate dimension exceeds budget")
        for row in raw:
            if not isinstance(row, list) or len(row) != width:
                raise ValueError("numerical coordinates must be rectangular")
            if any(type(x) not in (int, float) or not math.isfinite(x) for x in row):
                raise ValueError("numerical coordinates must be finite real scalars")
        real = np.asarray(raw, dtype=float)
        if len(np.unique(np.round(real, 7), axis=0)) != len(real):
            raise ValueError("duplicate numerical coordinates")
        if np.max(abs(np.sum(real**2, axis=1) - 1)) > 1e-8:
            raise ValueError("incorrect numerical norm")
        maximum = max_cosine(real.astype(complex))
        if maximum > 0.5 + 1e-8:
            raise ValueError("numerical pair exceeds kissing threshold")
        rank = int(np.linalg.matrix_rank(real, tol=1e-7))
        if "cardinality" in artifact.facts and artifact.facts["cardinality"] != len(real):
            raise ValueError("claimed cardinality mismatch")
        if "dimension" in artifact.facts and artifact.facts["dimension"] != rank:
            raise ValueError("claimed dimension mismatch")
        for key, value in {"ambient_dimension": real.shape[1], "norm_squared": 1, "max_inner_product": maximum}.items():
            if key in artifact.facts and artifact.facts[key] != value:
                raise ValueError("claimed " + key + " mismatch")
        return Verification(status="numerical", method="full floating point Gram check", details={
            "cardinality": len(real), "dimension": rank, "rank_numerical": rank,
            "max_cosine": maximum, "tolerance": 1e-8, "certificate": False})
    except (ValueError, KeyError, TypeError, IndexError, OverflowError) as error:
        return Verification(status="unverified", method="group artifact verification failed", details={"reason": str(error)})


OPERATORS = {"groups.weil": weil, "groups.search": search, "groups.ganzhinov10": ganzhinov10,
             "groups.ganzhinov14": ganzhinov14, "groups.from_generators": from_generators}
INPUT_KINDS = {"groups.weil": (), "groups.search": ("unitary_representation",), "groups.ganzhinov10": (),
               "groups.ganzhinov14": (), "groups.from_generators": ()}
