"""Finite exact compatibility of explicit coordinate and rational Gram blocks.

Cross inner products are supplied mathematical data, never guessed from two
unrelated frames. The entire unscaled source block matrix must be PSD before
scaling, adding a shared orthogonal head space, or identifying overlaps.
"""
from fractions import Fraction
from pathlib import Path

from .core import Artifact, RunContext, Verification, canonical_json, digest
from .exact_gram import (bounded_fraction, check_fraction_size, positive_semidefinite_rank,
                         rational_json, verify_rational_gram)


SUPPORTED_KINDS = ("spherical_code", "gram_spherical_code")


def _context(context):
    return context or RunContext(Path("."))


def _source_count(source, context):
    if not isinstance(source, Artifact) or source.kind not in SUPPORTED_KINDS:
        raise ValueError("Unsupported compatibility backend: exact integer coordinates or rational Gram required; implicit/interval backends are not implemented")
    key = "coordinates" if source.kind == "spherical_code" else "gram"
    rows = source.data.get(key)
    if not isinstance(rows, list) or not rows:
        raise ValueError("Compatibility source needs nonempty explicit coordinates or Gram")
    count = len(rows)
    context.limit("max_materialized_points", count)
    context.limit("max_gram_entries", count * count)
    context.limit("max_pairs", count * (count - 1) // 2)
    return count


def _replay_depth(sources, context):
    """Bound nested gluing evidence before invoking recursive source replay."""
    pending = [{"data": source.data} for source in sources]
    visited = 0
    while pending:
        source = pending.pop()
        visited += 1
        context.limit("max_nodes", visited)
        if not isinstance(source, dict) or not isinstance(source.get("data"), dict):
            raise ValueError("Malformed embedded compatibility source")
        recipe = source["data"].get("gluing")
        if recipe is not None:
            children = recipe.get("sources") if isinstance(recipe, dict) else None
            if not isinstance(children, list) or len(children) != 2:
                raise ValueError("Malformed nested gluing sources")
            pending.extend(children)


def _source_gram(source, context):
    from .certification import verify_artifact
    count = _source_count(source, context)
    if source.kind == "spherical_code":
        rows = source.data["coordinates"]
        if not isinstance(rows[0], list) or not rows[0]:
            raise ValueError("Invalid coordinate source shape")
        width = len(rows[0])
        # Three scalar operations per weighted product term, before Gram allocation.
        context.limit("max_rank_work", 3 * width * count * (count + 1) // 2)
    report = verify_artifact(source, context)
    if report.status != "exact":
        raise ValueError("Compatibility requires independently exact sources")
    if source.kind == "gram_spherical_code":
        matrix = [[bounded_fraction(value, context, "source Gram entry") for value in row]
                  for row in source.data["gram"]]
    else:
        rows = source.data["coordinates"]
        weights = source.data.get("weights", [1] * len(rows[0]))
        matrix = [[Fraction(0) for _ in range(count)] for _ in range(count)]
        for i, first in enumerate(rows):
            for j in range(i + 1):
                value = check_fraction_size(Fraction(sum(w * x * y for w, x, y in zip(weights, first, rows[j]))), context)
                matrix[i][j] = matrix[j][i] = value
    return matrix, report


def _maps(scales, heads, context):
    if scales is None:
        scales = [1, 1]
    if heads is None:
        heads = [[], []]
    if not isinstance(scales, list) or len(scales) != 2:
        raise ValueError("scales must contain two exact rational values")
    values = [bounded_fraction(value, context, "gluing scale") for value in scales]
    if any(value == 0 for value in values):
        raise ValueError("Unsupported noninjective zero-scale map; gluing requires nonzero rational scales")
    if not isinstance(heads, list) or len(heads) != 2 or any(not isinstance(row, list) for row in heads) or len(heads[0]) != len(heads[1]):
        raise ValueError("heads must contain two rational vectors of equal width")
    context.limit("max_dimension", len(heads[0]))
    context.limit("max_materialized_points", 2 * len(heads[0]))
    vectors = [[bounded_fraction(value, context, "head coordinate") for value in row] for row in heads]
    return values, vectors


def _dot(first, second, context):
    value = Fraction(0)
    for x, y in zip(first, second):
        value = check_fraction_size(value + check_fraction_size(x * y, context), context)
    return value


def _transformed(matrix, scale, head, context):
    square = check_fraction_size(scale * scale, context)
    height = _dot(head, head, context)
    return [[check_fraction_size(square * value + height, context) for value in row] for row in matrix]


def within(source, context=None, *, scale=1, head=None):
    """Independently certify the mapped within-source geometry (nonzero scale)."""
    context = _context(context)
    _source_count(source, context)
    _replay_depth([source], context)
    vector = [] if head is None else head
    values, heads = _maps([scale, 1], [vector, vector], context)
    matrix, source_report = _source_gram(source, context)
    transformed = _transformed(matrix, values[0], heads[0], context)
    report = verify_rational_gram({"gram": [[rational_json(value) for value in row] for row in transformed]}, context)
    report.method = "exact scaled source Gram with rational orthogonal head and within-block checks"
    report.details["source_verification"] = source_report.to_dict()
    report.details["map"] = {"scale": rational_json(values[0]), "head": [rational_json(x) for x in heads[0]], "injective": True}
    return report


def _assemble(first, second, cross_gram, scales, heads, context):
    count_a, count_b = _source_count(first, context), _source_count(second, context)
    total = count_a + count_b
    context.limit("max_materialized_points", total)
    context.limit("max_gram_entries", total * total)
    context.limit("max_pairs", total * (total - 1) // 2)
    _replay_depth([first, second], context)
    values, vectors = _maps(scales, heads, context)
    if not isinstance(cross_gram, list) or len(cross_gram) != count_a or any(not isinstance(row, list) or len(row) != count_b for row in cross_gram):
        raise ValueError("cross_gram must have one row per first-source point and one column per second-source point")
    cross = [[bounded_fraction(value, context, "cross Gram entry") for value in row] for row in cross_gram]
    gram_a, report_a = _source_gram(first, context)
    gram_b, report_b = _source_gram(second, context)
    base = [row[:] + cross[i][:] for i, row in enumerate(gram_a)]
    base.extend([[cross[i][j] for i in range(count_a)] + row[:] for j, row in enumerate(gram_b)])
    # This check occurs before heads can conceal an inconsistent relative frame,
    # and before duplicate removal can conceal inconsistent repeated rows.
    base_rank, base_work = positive_semidefinite_rank(base, context)
    mapped_a = _transformed(gram_a, values[0], vectors[0], context)
    mapped_b = _transformed(gram_b, values[1], vectors[1], context)
    if mapped_a[0][0] != mapped_b[0][0]:
        raise ValueError("Scaled/headed source blocks do not have the same norm")
    norm = mapped_a[0][0]
    coefficient = check_fraction_size(values[0] * values[1], context)
    head_dot = _dot(vectors[0], vectors[1], context)
    mapped_cross = [[check_fraction_size(coefficient * value + head_dot, context) for value in row] for row in cross]
    full = [row[:] + mapped_cross[i][:] for i, row in enumerate(mapped_a)]
    full.extend([[mapped_cross[i][j] for i in range(count_a)] + row[:] for j, row in enumerate(mapped_b)])
    full_rank, full_work = positive_semidefinite_rank(full, context)

    # In a PSD common-norm Gram, g_ij=norm iff ||x_i-x_j||^2=0.
    # The PSD check above prevents a fake overlap from hiding inconsistent rows.
    parents = list(range(total))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    for i in range(total):
        for j in range(i):
            if full[i][j] == norm:
                a, b = find(i), find(j)
                parents[max(a, b)] = min(a, b)
    classes = {}
    for index in range(total):
        classes.setdefault(find(index), []).append(index)
    representatives = list(classes)
    unique = [[full[i][j] for j in representatives] for i in representatives]
    plain = {"gram": [[rational_json(value) for value in row] for row in unique], "norm_squared": rational_json(norm)}
    report = verify_rational_gram(plain, context)
    if report.details["facts"]["rank_exact"] != full_rank:
        raise ValueError("Deduplication unexpectedly changed Gram span")
    labels = lambda index: [0, index] if index < count_a else [1, index - count_a]
    evidence = {"sources": [first.to_dict(), second.to_dict()],
                "source_sha256": [digest(first.to_dict()), digest(second.to_dict())],
                "scales": [rational_json(value) for value in values],
                "heads": [[rational_json(value) for value in row] for row in vectors],
                "cross_gram": [[rational_json(value) for value in row] for row in cross],
                "classes": [[labels(index) for index in members] for members in classes.values()]}
    data = {**plain, "gluing": evidence}
    facts = dict(report.details["facts"])
    facts.update(input_cardinalities=[count_a, count_b], duplicates_removed=total - len(unique), global_optimality_claim=False)
    report.method = "exact heterogeneous Gram gluing: source replay, joint PSD, injective maps, exact overlap quotient and final Gram checks"
    report.details.update(facts=facts, **facts, data_sha256=digest(data),
                          source_verifications=[report_a.to_dict(), report_b.to_dict()],
                          joint_base_rank=base_rank, joint_base_rank_scalar_work=base_work,
                          mapped_rank_before_dedup=full_rank, mapped_rank_scalar_work=full_work,
                          input_ambient_dimensions=[report_a.details.get("facts", report_a.details).get("ambient_dimension"), report_b.details.get("facts", report_b.details).get("ambient_dimension")],
                          head_width=len(vectors[0]),
                          cross_frame_scope="provided unscaled cross inner products validated by the full joint PSD matrix; no relative frame inferred",
                          source_metadata_scope="source metadata retained; historical and novelty claims are not established by this geometric replay")
    return data, facts, report


def between(first, second, cross_gram, context=None, *, scales=None, heads=None):
    """Certify finite exact joint compatibility, including overlaps and rank."""
    return _assemble(first, second, cross_gram, scales, heads, _context(context))[2]


def certified_union(first, second, cross_gram, context=None, *, scales=None, heads=None):
    """Return an independently replayable Gram union of the two exact sources."""
    data, facts, report = _assemble(first, second, cross_gram, scales, heads, _context(context))
    return Artifact("gram_spherical_code", data, facts, report)


def verify_glued_artifact(artifact, context=None):
    """Reconstruct everything from the stored source/maps/cross data, not status."""
    from .certification import compare_facts
    context = _context(context)
    if artifact.kind != "gram_spherical_code" or set(artifact.data) != {"gram", "norm_squared", "gluing"}:
        raise ValueError("Invalid glued Gram artifact schema")
    recipe = artifact.data["gluing"]
    if not isinstance(recipe, dict) or set(recipe) != {"sources", "source_sha256", "scales", "heads", "cross_gram", "classes"}:
        raise ValueError("Invalid gluing recipe fields")
    sources = recipe["sources"]
    if not isinstance(sources, list) or len(sources) != 2:
        raise ValueError("Gluing needs two embedded sources")
    # Count nested descriptors before Artifact.from_dict recursively normalizes.
    pending, nodes = list(sources), 0
    while pending:
        source = pending.pop()
        nodes += 1
        context.limit("max_nodes", nodes)
        if not isinstance(source, dict) or not isinstance(source.get("data"), dict):
            raise ValueError("Malformed stored gluing source")
        nested = source["data"].get("gluing")
        if nested is not None:
            if not isinstance(nested, dict) or not isinstance(nested.get("sources"), list) or len(nested["sources"]) != 2:
                raise ValueError("Malformed nested gluing recipe")
            pending.extend(nested["sources"])
    first, second = (Artifact.from_dict(source) for source in sources)
    data, facts, report = _assemble(first, second, recipe["cross_gram"], recipe["scales"], recipe["heads"], context)
    if canonical_json(data) != canonical_json(artifact.data):
        raise ValueError("Stored glued Gram, maps, source hashes or overlap classes differ from independent replay")
    compare_facts(artifact.facts, facts)
    return report
