"""Fresh classical root systems with exact integer weighted coordinates.

These are construction controls, not a catalogue of supplied record points.
A and D are rank-parametrized; E exists only at its mathematical ranks 6,7,8.
All norm, distinctness, pair and span claims are replayed by the shared checker.
"""
from itertools import combinations, product
from math import lcm
from pathlib import Path

from ..core import Artifact, RunContext, check_params, integer
from ..certification import verify_integer_spherical


def _type_d(rank, amplitude=1):
    for first, second in combinations(range(rank), 2):
        for left, right in product((-1, 1), repeat=2):
            vector = [0] * rank
            vector[first], vector[second] = amplitude * left, amplitude * right
            yield vector


def _type_e8():
    # Twice the conventional E8 roots: squared norm 8 in Z^8.
    yield from _type_d(8, amplitude=2)
    for signs in product((-1, 1), repeat=8):
        if signs.count(-1) % 2 == 0:
            yield list(signs)


def _type_a(rank):
    """Orthogonal Helmert coordinates on sum(x_i)=0, with no square roots.

    b_k=(1,...,1,-k,0,...) has norm k(k+1), and the b_k are
    orthogonal. For x in the zero-sum subspace, sum_k (x.b_k)^2 /
    [k(k+1)] = ||x||^2. Multiplying the metric by the common L keeps
    all coordinates and weights integral and preserves normalized angles.
    """
    denominators = [k * (k + 1) for k in range(1, rank + 1)]
    common = lcm(*denominators)
    coordinates = []
    for first in range(rank + 1):
        for second in range(rank + 1):
            if first == second:
                continue
            vector = [0] * (rank + 1)
            vector[first], vector[second] = 1, -1
            partial = vector[0]
            row = []
            for k in range(1, rank + 1):
                row.append(partial - k * vector[k])
                partial += vector[k]
            coordinates.append(row)
    return coordinates, [common // value for value in denominators], 2 * common


def root_system(inputs, params, context):
    """Generate A_rank, D_rank or E_rank; no dimension-dependent dispatch."""
    check_params(params, ("family", "rank"))
    if inputs:
        raise ValueError("lattice.root_system takes no input artifacts")
    family = params.get("family")
    if not isinstance(family, str) or family not in ("A", "D", "E"):
        raise ValueError("family must be A, D or E")
    rank = integer(params.get("rank"), "root-system rank")
    context = context or RunContext(Path("."))
    context.limit("max_dimension", rank)
    if family == "D" and rank < 2:
        raise ValueError("D roots require rank at least 2")
    if family == "E" and rank not in (6, 7, 8):
        raise ValueError("E roots exist in ranks 6, 7 and 8")
    cardinality = rank * (rank + 1) if family == "A" else (
        2 * rank * (rank - 1) if family == "D" else {6: 72, 7: 126, 8: 240}[rank])
    # These precede all root-list, coordinate-matrix, or pair allocation.
    context.limit("max_materialized_points", cardinality)
    context.limit("max_pairs", cardinality * (cardinality - 1) // 2)
    if family == "A":
        coordinates, weights, norm = _type_a(rank)
        construction = "A roots e_i-e_j; exact orthogonal Helmert coordinates with integral diagonal metric"
    elif family == "D":
        coordinates, weights, norm = list(_type_d(rank)), [1] * rank, 2
        construction = "D roots: all four signs of e_i+e_j for distinct coordinate indices"
    else:
        coordinates = []
        for vector in _type_e8():
            if rank <= 7 and vector[6] + vector[7] != 0:
                continue
            if rank == 6 and vector[5] + vector[7] != 0:
                continue
            coordinates.append(vector[:rank])
        weights = {6: [1] * 5 + [3], 7: [1] * 6 + [2], 8: [1] * 8}[rank]
        norm = 8
        construction = ("E8 integer roots: signed 2e_i+2e_j and even all-sign vectors; "
                        "E7/E6 exact rational hyperplane sections with isometric weighted coordinate projection")
    if len(coordinates) != cardinality:
        raise ValueError("Root generator cardinality disagrees with its formula")
    data = {"coordinates": coordinates, "weights": weights, "norm_squared": norm,
            "construction": construction, "root_system": {"family": family, "rank": rank},
            "provenance": "classical root-system formula replay; no imported record coordinates or novelty claim"}
    report = verify_integer_spherical(data, context)
    if report.status != "exact" or report.details.get("facts", report.details).get("dimension") != rank:
        raise ValueError("Root system did not independently verify its intended exact rank")
    facts = dict(report.details.get("facts", report.details))
    facts.update(root_system=f"{family}{rank}", classification="classical_construction_control",
                 global_optimality_claim=False)
    return Artifact("spherical_code", data, facts, report)


OPERATORS = {"lattice.root_system": root_system}
INPUT_KINDS = {"lattice.root_system": ()}
