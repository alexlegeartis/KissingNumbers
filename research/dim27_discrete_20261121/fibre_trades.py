"""Search signed trades in a single native fibre, allowing exterior point deletions."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import numpy as np
from numba import njit
from ortools.sat.python import cp_model

from bw_refill import TAILS
from fibres import native_fibres
from midpoint_switches import source_layers, switch_groups
from source import COMPLEMENT, OCTAD, checked_archive, leech_shell, load_points


def initialize():
    global shell, mids, heads, labels, groups, plane, other_line, roots_by_fibre
    global owner, deleted, incident, body_matrix, head_matrix
    shell = 2 * leech_shell().astype(np.int16)
    mids, labels, all_heads = source_layers()
    heads = all_heads[all_heads[np.arange(len(all_heads)), 24 + np.argmax(all_heads[:, 24:] != 0, axis=1)] > 0]
    groups, _, _ = switch_groups(mids, labels, all_heads)
    plane = np.argmin(np.abs(TAILS), axis=1)
    other_line = np.array([np.flatnonzero((plane == p) & (np.arange(6) != i))[0] for i, p in enumerate(plane)])
    frames, _ = native_fibres(shell)
    number = {w.tobytes(): i for i, w in enumerate(frames.reshape(-1, 16))}
    fibre = np.full(len(shell), -1, dtype=np.int32)
    for i in np.flatnonzero((shell[:, COMPLEMENT] ** 2).sum(1) == 64):
        fibre[i] = number[shell[i, COMPLEMENT].tobytes()]
    roots_by_fibre = [np.flatnonzero(fibre == f) for f in range(len(number))]
    assert all(len(r) == 16 for r in roots_by_fibre)
    with checked_archive("d27_certificate.npz") as data:
        raw = data["added_integer"][data["midpoint_rows"]]
        pairs = np.lexsort(shell.T[::-1])[data["midpoint_endpoints"]]
    parents = {v.tobytes(): p for v, p in zip(raw, pairs)}
    pairs = np.array([parents[v.tobytes()] for v in mids])
    owner = np.where(fibre[pairs[:, 0]] == fibre[pairs[:, 1]], fibre[pairs[:, 0]], -1)
    deleted = np.zeros(len(shell), dtype=bool)
    deleted[pairs] = True
    incident = [[] for _ in shell]
    for i, (a, b) in enumerate(pairs):
        incident[a].append(i)
        incident[b].append(i)
    body_matrix = mids[:, :24].astype(np.float32).T
    head_matrix = heads[:, :24].astype(np.float32).T


@njit
def price_planes(products, owner, fibre, groups, labels, plane):
    """Minimum signed-midpoint cuts for one candidate, before head costs."""
    result = np.zeros((len(products), 3), dtype=np.int32)
    for i, row in enumerate(products):
        # For a fixed plane, cost is independent of the new line up to switching every group in that plane.
        same = np.zeros((32, 6), dtype=np.int32)
        perp = np.zeros((32, 6), dtype=np.int32)
        cross = np.zeros(3, dtype=np.int32)
        for m, score in enumerate(row):
            if owner[m] == fibre or score <= 32:
                continue
            c = labels[m]
            p = plane[c]
            same[groups[m], c] += 1 + int(score > 96)
            perp[groups[m], c] += 2 * int(score > 64)
            cross[p] += int(score > 48) + int(score > 80)
        for p in range(3):
            colors = np.where(plane == p)[0]
            a, b = colors[0], colors[1]
            total = cross.sum() - cross[p]
            for j in range(32):
                total += min(same[j, a] + perp[j, b], same[j, b] + perp[j, a])
            result[i, p] = total
    return result


def add_candidate(model, phase, paid, rows, variables, point, required, conditional, budget):
    minimum = len(required) + sum(min(map(len, pair)) for pair in conditional.values())
    if minimum > budget:
        return
    selected = model.new_bool_var("")
    rows.append(point)
    variables.append(selected)

    def cut(row):
        if row not in paid:
            paid[row] = model.new_bool_var("")
        return paid[row]

    for row in required:
        model.add_implication(selected, cut(row))
    for group, pair in conditional.items():
        for state, blockers in enumerate(pair):
            wrong_phase = phase[group] if state == 0 else ~phase[group]
            if minimum - min(map(len, pair)) + len(blockers) > budget:
                model.add_implication(selected, wrong_phase)
            else:
                for row in blockers:
                    model.add_bool_or([~selected, wrong_phase, cut(row)])


def search(fibre):
    roots = roots_by_fibre[fibre]
    mcount = len(mids)
    incumbent = int(np.sum(~deleted[roots]) + 2 * np.sum(owner == fibre))
    # Residuals have norm 64 and pairwise nonpositive products in dimension 11.
    budget = 22 - incumbent - 1
    a, b = np.where(np.triu(shell[roots] @ shell[roots].T == 64, 1))
    bodies = (shell[roots[a]] + shell[roots[b]]) // 2
    assert len(bodies) == 112
    products = bodies.astype(np.float32) @ body_matrix
    hproducts = bodies.astype(np.float32) @ head_matrix
    prices = price_planes(products, owner, fibre, groups, labels, plane)
    model = cp_model.CpModel()
    phase = [model.new_bool_var("") for _ in range(32)]
    paid, rows, variables = {}, [], []
    for i, body in enumerate(bodies):
        if np.min(prices[i]) > budget:
            continue
        near = np.flatnonzero((products[i] > 32) & (owner != fibre))
        scores, old_lines = products[i, near], labels[near]
        for c, tail in enumerate(TAILS):
            hd = tail @ heads[:, 24:].T
            hpos, hneg = hproducts[i] + hd > 64, hproducts[i] - hd > 64
            if prices[i, plane[c]] + hpos.sum() + hneg.sum() > budget:
                continue
            original = tail @ TAILS[old_lines].T
            switched = tail @ TAILS[other_line[old_lines]].T
            for sign in (1, -1):
                pos = hpos if sign == 1 else hneg
                neg = hneg if sign == 1 else hpos
                required = list(2 * mcount + np.flatnonzero(pos))
                required.extend(2 * mcount + len(heads) + np.flatnonzero(neg))
                conditional = {}
                for old_sign, offset in ((1, 0), (-1, mcount)):
                    before = scores + sign * old_sign * original > 64
                    after = scores + sign * old_sign * switched > 64
                    required.extend(offset + near[before & after])
                    for k in np.flatnonzero(before != after):
                        pair = conditional.setdefault(groups[near[k]], ([], []))
                        pair[int(after[k])].append(offset + near[k])
                add_candidate(model, phase, paid, rows, variables, np.r_[body[OCTAD], sign * tail],
                              required, conditional, budget)
    for r in roots:
        blockers = [m for m in incident[r] if owner[m] != fibre]
        add_candidate(model, phase, paid, rows, variables, np.r_[shell[r, OCTAD], [0, 0, 0]],
                      blockers + [m + mcount for m in blockers], {}, budget)
    if len(rows) <= incumbent:
        return "count"
    v = np.array(rows, dtype=np.int16)
    assert np.all((v * v).sum(1) == 64)
    for i, j in np.argwhere(np.triu(v @ v.T > 0, 1)):
        model.add_at_most_one(variables[i], variables[j])
    # A nonpositive code with m points and d dimensions has at least 2m-3d antipodal pairs.
    active_pairs = []
    lookup = {tuple(row): i for i, row in enumerate(v)}
    for i, row in enumerate(v):
        j = lookup.get(tuple(-row), -1)
        if i < j:
            used = model.new_bool_var("")
            model.add_implication(used, variables[i])
            model.add_implication(used, variables[j])
            active_pairs.append(used)
    selected = sum(variables)
    model.add(selected <= 22)
    model.add(sum(active_pairs) >= 2 * selected - 33)
    model.add(selected - sum(paid.values()) > incumbent)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(model)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        source = load_points()
        # Reconstruct the entire witness, including the chosen switches outside this fibre.
        remaining = source[(source[:, 24:] ** 2).sum(1) != 32]
        roots_here = {r.tobytes() for r in shell[roots]}
        remaining = [r for r in remaining if np.any(r[24:]) or r[:24].tobytes() not in roots_here]
        switched = mids.copy()
        mask = np.array([solver.boolean_value(x) for x in phase])
        switched[:, 24:] = TAILS[np.where(mask[groups], other_line[labels], labels)]
        minus = switched.copy()
        minus[:, 24:] *= -1
        hminus = heads.copy()
        hminus[:, 24:] *= -1
        mutable = np.vstack([switched, minus, heads, hminus])
        drop = {mutable[i].tobytes() for i, x in paid.items() if solver.boolean_value(x)}
        kept = [r for r in remaining if r.tobytes() not in drop]
        kept.extend(r for i, r in enumerate(mutable[:2*mcount]) if owner[i % mcount] != fibre and r.tobytes() not in drop)
        local = v[[solver.boolean_value(x) for x in variables]]
        added = np.zeros((len(local), 27), dtype=np.int16)
        added[:, COMPLEMENT] = shell[roots[0], COMPLEMENT]
        added[:, OCTAD] = local[:, :8]
        added[:, 24:] = local[:, 8:]
        points = np.vstack([kept, added])
        assert len(points) > len(source)
        np.savez_compressed(output / f"fibre_{fibre}.npz", points=points)
    return solver.status_name(status)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seconds", type=float, default=20)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fibre", type=int)
    args = parser.parse_args()
    seconds, output = args.seconds, args.out
    output.mkdir(exist_ok=True, parents=True)
    initialize()
    ids = range(4320) if args.fibre is None else [args.fibre]
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=get_context("fork")) as pool:
        print(dict(Counter(pool.map(search, ids, chunksize=4))), flush=True)
