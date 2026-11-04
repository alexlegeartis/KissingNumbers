"""Exact point selection from two explicitly realized isometric code copies.

The optimizer works on points, including accessory points.  It never assumes
that a selected subset is a union of group orbits or product fibers.  Geometry
of the returned ordinary integer-coordinate artifact is checked independently.
"""
from collections import deque
from fractions import Fraction
from math import gcd, isqrt, lcm

import numpy as np

from ..core import Artifact, check_params, digest, integer
from ..exact_gram import bounded_fraction, check_fraction_size


def _checked(value, context):
    return check_fraction_size(Fraction(value), context)


def rational_isometry(dimension, weights, params, context):
    """Return (integer row-map numerator, denominator) with A W A^t = W.

    Exactly one of ``matrix`` and ``normals`` is required.  Matrix entries and
    Householder normal entries use rational JSON scalars.  Rows act on the
    right: ``image = coordinates @ matrix``.  A weighted Householder normal v
    acts by x -> x - 2*(x W v)/(v W v)*v.
    """
    if ("matrix" in params) == ("normals" in params):
        raise ValueError("Provide exactly one of matrix or normals")
    context.limit("max_gram_entries", dimension * dimension)
    if "matrix" in params:
        raw = params["matrix"]
        if not isinstance(raw, list) or len(raw) != dimension or any(
                not isinstance(row, list) or len(row) != dimension for row in raw):
            raise ValueError("Isometry matrix must match the source coordinate width")
        matrix = [[bounded_fraction(x, context, "isometry entry") for x in row] for row in raw]
    else:
        normals = params["normals"]
        if not isinstance(normals, list):
            raise ValueError("normals must be a list of rational normal vectors")
        context.limit("max_nodes", len(normals))
        context.limit("max_rank_work", (len(normals) + dimension) * dimension ** 2 * 3)
        matrix = [[Fraction(i == j) for j in range(dimension)] for i in range(dimension)]
        for raw in normals:
            if not isinstance(raw, list) or len(raw) != dimension:
                raise ValueError("Householder normal must match the source coordinate width")
            normal = [bounded_fraction(x, context, "Householder normal") for x in raw]
            if not any(normal):
                raise ValueError("Householder normal must be nonzero")
            square = sum((_checked(w * x * x, context) for w, x in zip(weights, normal)), Fraction(0))
            _checked(square, context)
            for row in matrix:
                dot = _checked(sum((_checked(x * w * v, context) for x, w, v in zip(row, weights, normal)), Fraction(0)), context)
                factor = _checked(2 * dot / square, context)
                row[:] = [_checked(x - _checked(factor * v, context), context) for x, v in zip(row, normal)]
    context.limit("max_rank_work", dimension ** 3 * 3)
    for i in range(dimension):
        for j in range(i, dimension):
            dot = _checked(sum((_checked(matrix[i][k] * weights[k] * matrix[j][k], context)
                                for k in range(dimension)), Fraction(0)), context)
            if dot != (weights[i] if i == j else 0):
                raise ValueError("Rational map does not preserve the exact weighted metric")
    denominator = 1
    for row in matrix:
        for value in row:
            denominator = lcm(denominator, value.denominator)
            context.limit("max_rational_bits", denominator.bit_length())
    numerator = [[int(value * denominator) for value in row] for row in matrix]
    for row in numerator:
        for value in row:
            context.limit("max_rational_bits", value.bit_length())
    return numerator, denominator


def transform_rows(rows, numerator, denominator, context, image_rows=None):
    """Common-scale integer coordinates; never round a rational image."""
    dimension = len(numerator)
    image_rows = rows if image_rows is None else image_rows
    # Bound every scalar product before NumPy can choose a fixed-width dtype.
    largest = max(abs(x) for row in rows + image_rows for x in row)
    map_bound = max(sum(abs(numerator[k][j]) for k in range(dimension)) for j in range(dimension))
    safe = largest * map_bound < 2 ** 62 and largest * denominator < 2 ** 62
    dtype = np.int64 if safe else object
    source = np.asarray(rows, dtype=dtype)
    image = np.asarray(image_rows, dtype=dtype) @ np.asarray(numerator, dtype=dtype)
    left = source * denominator
    left_rows = [[int(x) for x in row] for row in left]
    right_rows = [[int(x) for x in row] for row in image]
    for row in left_rows + right_rows:
        for value in row:
            context.limit("max_rational_bits", value.bit_length())
    return left_rows, right_rows


def conflict_graph(left, right, weights, norm, context):
    """All cross conflicts using bounded row blocks and exact integer dots."""
    context.limit("max_pairs", len(left) * len(right))
    if not right:
        return [[] for _ in left], "int64 proven overflow-safe"
    context.limit("max_gram_entries", min(64, len(left)) * len(right))
    bound = max(abs(x) for row in left for x in row) * max(abs(x) for row in right for x in row) * sum(weights)
    safe = bound < 2 ** 62 and max(weights) < 2 ** 62
    dtype = np.int64 if safe else object
    second = np.asarray(right, dtype=dtype)
    weight_array = np.asarray(weights, dtype=dtype)
    graph = []
    edges = 0
    for start in range(0, len(left), 64):
        first = np.asarray(left[start:start + 64], dtype=dtype)
        dots = (first * weight_array) @ second.T
        block = [np.flatnonzero(row > norm // 2).tolist() for row in dots]
        edges += sum(map(len, block))
        context.limit("max_gram_entries", edges)
        graph.extend(block)
    return graph, "int64 proven overflow-safe" if safe else "Python integers"


def maximum_independent_bipartite(graph, right_count):
    """Matching plus a checked equal-size vertex cover certifies optimality.

    The returned left/right indices are the complement of that cover.  The
    solver's result is not trusted: matching and cover properties are checked
    against every edge of the supplied exact graph.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching
    lengths = [len(row) for row in graph]
    indptr = np.r_[0, np.cumsum(lengths)]
    indices = np.asarray([v for row in graph for v in row], dtype=np.int64)
    sparse = csr_matrix((np.ones(len(indices), dtype=np.int8), indices, indptr), shape=(len(graph), right_count))
    matching = maximum_bipartite_matching(sparse, perm_type="column").tolist()
    reverse = [-1] * right_count
    for u, v in enumerate(matching):
        if v < 0:
            continue
        if v >= right_count or v not in graph[u] or reverse[v] >= 0:
            raise ValueError("Invalid bipartite matching result")
        reverse[v] = u
    reachable_left = {u for u, v in enumerate(matching) if v < 0}
    reachable_right = set()
    queue = deque(sorted(reachable_left))
    while queue:
        u = queue.popleft()
        for v in graph[u]:
            if v == matching[u] or v in reachable_right:
                continue
            reachable_right.add(v)
            if reverse[v] < 0:
                raise ValueError("Matching has an augmenting path")
            other = reverse[v]
            if other not in reachable_left:
                reachable_left.add(other)
                queue.append(other)
    cover_left = sorted(set(range(len(graph))) - reachable_left)
    cover_right = sorted(reachable_right)
    cover_left_set, cover_right_set = set(cover_left), set(cover_right)
    if any(u not in cover_left_set and v not in cover_right_set for u, row in enumerate(graph) for v in row):
        raise ValueError("Bipartite cover misses a conflict edge")
    pairs = [[u, v] for u, v in enumerate(matching) if v >= 0]
    if len(cover_left) + len(cover_right) != len(pairs):
        raise ValueError("Matching and cover cardinalities disagree")
    return {"left_indices": sorted(reachable_left),
            "right_indices": sorted(set(range(right_count)) - reachable_right),
            "matching": pairs, "cover_left": cover_left, "cover_right": cover_right,
            "optimal_in_supplied_graph": True}


def select_coupled_rows(rows, image_rows, weights, norm, params, context):
    """Search primitive; caller must first verify both codes' within geometry.

    Sources must already have the same norm and diagonal weighted metric.
    """
    count, image_count, dimension = len(rows), len(image_rows), len(rows[0])
    total = count + image_count
    context.limit("max_materialized_points", total)
    context.limit("max_dimension", dimension)
    context.limit("max_pairs", total * (total - 1) // 2)
    numerator, denominator = rational_isometry(dimension, weights, params, context)
    left, full_right = transform_rows(rows, numerator, denominator, context, image_rows)
    left_set = set(map(tuple, left))
    right_indices = [i for i, row in enumerate(full_right) if tuple(row) not in left_set]
    right = [full_right[i] for i in right_indices]
    scaled_norm = norm * denominator ** 2
    context.limit("max_rational_bits", scaled_norm.bit_length())
    graph, backend = conflict_graph(left, right, weights, scaled_norm, context)
    certificate = maximum_independent_bipartite(graph, len(right))
    chosen = [left[i] for i in certificate["left_indices"]] + [right[i] for i in certificate["right_indices"]]
    # Common gcd reduction changes neither normalized geometry nor matching.
    common = 0
    for row in chosen:
        for value in row:
            common = gcd(common, value)
    if common > 1:
        if scaled_norm % (common * common):
            raise ValueError("Nonintegral norm after common coordinate reduction")
        chosen = [[x // common for x in row] for row in chosen]
        scaled_norm //= common * common
    receipt = {**certificate, "source_cardinality": count, "image_cardinality": image_count,
               "image_new_points": len(right), "overlap_points": image_count - len(right),
               "cross_pairs_checked": count * len(right),
               "conflict_edges": sum(map(len, graph)), "graph_sha256": digest(graph),
               "integer_backend": backend, "matrix_numerator": numerator,
               "matrix_denominator": denominator, "coordinate_gcd_reduction": common,
               "right_source_indices": right_indices,
               "scope": "maximum cardinality subset of this two-copy union; no global optimality or novelty claim"}
    return {"coordinates": chosen, "weights": weights[:], "norm_squared": scaled_norm}, receipt


def select_isometric_rows(rows, weights, norm, params, context):
    """Specialize coupled selection to a code and its rational isometric image."""
    return select_coupled_rows(rows, rows, weights, norm, params, context)


def isometric_exchange(inputs, params, context):
    """Typed exact operator returning a plain independently replayable code."""
    from ..certification import verify_artifact, verify_integer_spherical
    check_params(params, ("matrix", "normals"))
    if len(inputs) != 1 or inputs[0].kind != "spherical_code" or "coordinates" not in inputs[0].data:
        raise ValueError("selection.isometric_exchange needs one explicit weighted integer code")
    source = inputs[0]
    rows = source.data.get("coordinates")
    if not isinstance(rows, list) or not rows or not isinstance(rows[0], list) or not rows[0]:
        raise ValueError("Nonempty explicit integer source rows required")
    count, dimension = len(rows), len(rows[0])
    # Preflight the entire possible union before verification or matrix work.
    context.limit("max_materialized_points", 2 * count)
    context.limit("max_dimension", dimension)
    context.limit("max_pairs", 2 * count * (2 * count - 1) // 2)
    weights = source.data.get("weights", [1] * dimension)
    if not isinstance(weights, list) or len(weights) != dimension:
        raise ValueError("Wrong source diagonal quadratic form")
    for weight in weights:
        integer(weight, "quadratic-form weight")
        context.limit("max_rational_bits", weight.bit_length())
    # Reject a malformed/nonisometric map before spending source pair work.
    rational_isometry(dimension, weights, params, context)
    if verify_artifact(source, context).status != "exact":
        raise ValueError("An independently exact source is required")
    data, receipt = select_isometric_rows(rows, weights, source.data["norm_squared"], params, context)
    receipt["source_data_sha256"] = digest(source.data)
    report = verify_integer_spherical(data, context)
    report.details["isometric_exchange"] = receipt
    result = Artifact("spherical_code", data, report.details["facts"], report)
    return result


def coupled_union(inputs, params, context):
    """Maximum point subset of two explicit valid codes in one exact metric.

    Codes may have different cardinalities and coordinate alphabets. Rational
    common scaling aligns norms when their ratio is a rational square. A
    supplied map places the second code relative to the first; omitted maps
    use the given coordinate frames. Implicit/algebraic-source adapters are
    deliberately outside this explicit-coordinate implementation.
    """
    from ..certification import verify_artifact, verify_integer_spherical
    check_params(params, ("matrix", "normals"))
    if len(inputs) != 2 or any(source.kind != "spherical_code" or "coordinates" not in source.data for source in inputs):
        raise ValueError("selection.coupled_union needs two explicit integer codes")
    raw_rows = [source.data.get("coordinates") for source in inputs]
    if any(not isinstance(rows, list) or not rows or not isinstance(rows[0], list) or not rows[0] for rows in raw_rows):
        raise ValueError("Nonempty explicit integer source rows required")
    dimension = len(raw_rows[0][0])
    if len(raw_rows[1][0]) != dimension:
        raise ValueError("Coupled sources require equal ambient coordinate widths")
    total = sum(map(len, raw_rows))
    context.limit("max_materialized_points", total)
    context.limit("max_dimension", dimension)
    context.limit("max_pairs", total * (total - 1) // 2)
    weights = [source.data.get("weights", [1] * dimension) for source in inputs]
    if weights[0] != weights[1]:
        raise ValueError("Coupled sources require the same explicit diagonal metric")
    for metric in weights:
        if not isinstance(metric, list) or len(metric) != dimension:
            raise ValueError("Wrong source diagonal quadratic form")
        for value in metric:
            integer(value, "quadratic-form weight")
            context.limit("max_rational_bits", value.bit_length())
    params = params or {"normals": []}
    rational_isometry(dimension, weights[0], params, context)
    norms = [integer(source.data.get("norm_squared"), "norm_squared") for source in inputs]
    for value in norms:
        context.limit("max_rational_bits", value.bit_length())
    common = gcd(*norms)
    squares = [norms[1] // common, norms[0] // common]
    scales = [isqrt(value) for value in squares]
    if any(scale * scale != value for scale, value in zip(scales, squares)):
        raise ValueError("Source norm ratio must be a rational square for integer-coordinate common scaling")
    aligned_norm = norms[0] * scales[0] ** 2
    context.limit("max_rational_bits", aligned_norm.bit_length())
    for source in inputs:
        if verify_artifact(source, context).status != "exact":
            raise ValueError("Independently exact sources required")
    scaled_rows = [[[value * scale for value in row] for row in rows] for rows, scale in zip(raw_rows, scales)]
    data, receipt = select_coupled_rows(scaled_rows[0], scaled_rows[1], weights[0], aligned_norm, params, context)
    receipt["source_data_sha256"] = [digest(source.data) for source in inputs]
    receipt["source_norm_alignment_scales"] = scales
    report = verify_integer_spherical(data, context)
    report.details["coupled_union"] = receipt
    return Artifact("spherical_code", data, report.details["facts"], report)


OPERATORS = {"selection.isometric_exchange": isometric_exchange, "selection.coupled_union": coupled_union}
INPUT_KINDS = {"selection.isometric_exchange": ("spherical_code",),
               "selection.coupled_union": ("spherical_code", "spherical_code")}
