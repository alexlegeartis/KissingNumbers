"""Solve a finite Barnes–Wall tail-pair refill exported by bw_refill.py."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ortools.sat.python import cp_model

from source import COMPLEMENT, ROOT


def add_conflicts(model, variables, bodies, tails):
    conflicts = 0
    b = bodies.astype(np.float32)
    t = tails.astype(np.float32)
    for start in range(0, len(b), 128):
        left = slice(start, start + 128)
        gram = b[left] @ b[start:].T + np.abs(t[left] @ t[start:].T)
        first, second = np.where(np.triu(gram > 64, 1))
        for i, j in zip(first, second):
            model.add_at_most_one(variables[start + int(i)], variables[start + int(j)])
        conflicts += len(first)
    return conflicts


def assemble(data, chosen):
    shell = data["shell"]
    parents = data["parents"][chosen]
    keep = ~np.isin(np.arange(len(shell)), np.union1d(data["deleted"], parents))
    equator = np.pad(shell[keep], ((0, 0), (0, 3)))
    insertions = np.zeros((int(chosen.sum()), 27), dtype=np.int16)
    insertions[:, COMPLEMENT] = data["bodies"][chosen]
    insertions[:, 24:] = data["tails"][data["colors"][chosen]]
    reflected = insertions.copy()
    reflected[:, 24:] *= -1
    return np.r_[equator, data["fixed"], insertions, reflected]


def solve(data, target, seconds, workers, output):
    parents = data["parents"]
    colors = data["colors"]
    old = data["incumbent"]
    model = cp_model.CpModel()
    choices = [model.new_bool_var(f"edge_{i}") for i in range(len(parents))]
    incidence = {}
    for i, row in enumerate(parents):
        for parent in row:
            incidence.setdefault(int(parent), []).append(i)
    twice = []
    for parent, indices in incidence.items():
        used = model.new_bool_var(f"used_{parent}")
        degree_two = model.new_bool_var(f"degree_two_{parent}")
        variables = [choices[i] for i in indices]
        model.add_max_equality(used, variables)
        model.add(sum(variables) == used + degree_two)
        # Two orbits of the same tail line at a shared parent conflict.
        for color in np.unique(colors[indices]):
            same_color = [choices[i] for i in indices if colors[i] == color]
            if len(same_color) > 1:
                model.add_at_most_one(same_color)
        degree = int(old[indices].sum())
        model.add_hint(used, degree > 0)
        model.add_hint(degree_two, degree == 2)
        twice.append(degree_two)

    for variable, present in zip(choices, old):
        model.add_hint(variable, bool(present))
    base = len(data["shell"]) - len(data["deleted"]) + len(data["fixed"])
    score = sum(twice)
    model.add(score + base >= target)
    model.maximize(score)
    conflicts = add_conflicts(model, choices, data["bodies"], data["tails"][colors])
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.symmetry_level = 0
    outcome = solver.solve(model)

    result = {
        "tail_pairs": len(choices), "parent_endpoints": len(twice),
        "pair_conflicts": conflicts, "fixed_points_before_refill": base,
        "target": target, "status": solver.status_name(outcome),
    }
    if outcome in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        chosen = np.array([solver.boolean_value(x) for x in choices])
        points = assemble(data, chosen)
        assert len(points) == base + round(solver.objective_value)
        sys.path.insert(0, str(ROOT / "research" / "26_27_native_projection_obstructions_20261118"))
        from verify_exact import verify

        result["verified_geometry"] = verify(points)
        np.savez_compressed(output, points=points, squared_norm=128, inner_product_limit=64)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pool", type=Path)
    parser.add_argument("--target", type=int, default=210109)
    parser.add_argument("--seconds", type=int, default=240)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with np.load(args.pool, allow_pickle=False) as data:
        print(json.dumps(solve(data, args.target, args.seconds, args.workers, args.out), indent=2))
