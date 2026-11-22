"""Exact one- and two-class exchanges of signed native fibre midpoints."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from source import COMPLEMENT, OCTAD, RESIDUAL, leech_shell, load_points


def native_fibres(shell):
    norm = np.sum(shell[:, OCTAD] ** 2, axis=1)
    fibres = defaultdict(list)
    for row in shell[norm == 64]:
        fibres[tuple(row[COMPLEMENT])].append(tuple(row[OCTAD]))
    classes = defaultdict(list)
    for projection, frame in fibres.items():
        classes[tuple(sorted(frame))].append(projection)
    projections = np.array(list(classes.values()), dtype=np.int16)
    frames = np.array(list(classes.keys()), dtype=np.int16)
    assert projections.shape == (135, 32, 16)
    assert frames.shape == (135, 16, 8)
    return projections, frames


def identify_source_fibres(points, projections):
    index = {tuple(w): i for i, group in enumerate(projections) for w in group}
    tail = np.sum(points[:, 24:] ** 2, axis=1)
    body = np.sum(points[:, OCTAD] ** 2, axis=1)
    comp = np.sum(points[:, COMPLEMENT] ** 2, axis=1)
    native = (comp == 64) & (((tail == 0) & (body == 64)) | ((tail == 32) & (body == 32)))
    location = np.full(len(points), -1, dtype=np.int16)
    states = defaultdict(set)
    for i in np.flatnonzero(native):
        w = tuple(points[i, COMPLEMENT])
        location[i] = index[w]
        states[w].add(tuple(points[i, RESIDUAL]))
    exceptional = {
        c for c, group in enumerate(projections)
        if len({frozenset(states[tuple(w)]) for w in group}) != 1
    }
    incumbent = [states[tuple(group[0])] for group in projections]
    return location, exceptional, incumbent


def signed_candidates(frame):
    """Every native contact midpoint with a signed D3 tail, and all roots."""
    roots = np.pad(frame, ((0, 0), (0, 3)))
    i, j = np.where(np.triu(frame @ frame.T == 0, 1))
    midpoints = (frame[i] + frame[j]) // 2
    tails = np.array([t for t in product((-4, 0, 4), repeat=3) if sum(v != 0 for v in t) == 2])
    return np.r_[roots, np.c_[np.repeat(midpoints, len(tails), axis=0),
                             np.tile(tails, (len(midpoints), 1))]].astype(np.int16)


def fixed_projection_maxima(projections, points):
    """Collect max w · b for each actual exterior residual and each class."""
    residuals, index = np.unique(points[:, RESIDUAL], axis=0, return_inverse=True)
    coordinates = points[:, COMPLEMENT].astype(np.float32)
    maxima = np.full((len(projections), len(residuals)), -128, dtype=np.int16)
    for c, group in enumerate(projections):
        dot = (group.astype(np.float32) @ coordinates.T).max(axis=0)
        np.maximum.at(maxima[c], index, dot.astype(np.int16))
    return residuals, maxima


def prepare():
    points = load_points()
    projections, frames = native_fibres(2 * leech_shell().astype(np.int16))
    locations, exceptional, source = identify_source_fibres(points, projections)
    assert len(exceptional) == 2
    free_classes = np.array([c for c in range(len(frames)) if c not in exceptional])
    fixed = points[~np.isin(locations, free_classes)]
    normals, max_dot = fixed_projection_maxima(projections, fixed)

    options, incumbent = [], []
    for c, frame in enumerate(frames):
        if c in exceptional:
            options.append(np.empty((0, 11), dtype=np.int16))
            incumbent.append(np.empty(0, dtype=bool))
            continue
        rows = signed_candidates(frame)
        allowed = (rows.astype(np.float32) @ normals.astype(np.float32).T <= 64 - max_dot[c]).all(1)
        rows = rows[allowed]
        in_source = np.array([tuple(row) in source[c] for row in rows])
        assert in_source.sum() == len(source[c])
        options.append(rows)
        incumbent.append(in_source)

    w = projections.reshape(-1, 16).astype(np.float32)
    max_between = (w @ w.T).reshape(135, 32, 135, 32).max(axis=(1, 3)).astype(np.int16)
    assert len(fixed) + 32 * sum(x.sum() for x in incumbent) == len(points)
    return free_classes, options, incumbent, max_between


def blocker_sets(options, incumbent, max_between):
    old = np.concatenate([rows[mask] for rows, mask in zip(options, incumbent)])
    owners = np.repeat(np.arange(len(options)), [mask.sum() for mask in incumbent])
    old = old.astype(np.float32)
    blockers = []
    for c, rows in enumerate(options):
        bad = rows.astype(np.float32) @ old.T > 64 - max_between[c, owners]
        bad[:, owners == c] = False
        blockers.append([frozenset(owners[row].tolist()) for row in bad])
    return blockers


def has_independent_set(bad, target):
    """Bitset branch and bound. Each disjoint conflict costs one vertex."""
    neighbors = [sum(1 << j for j in np.flatnonzero(row).tolist()) for row in bad]

    def search(rest, needed):
        while needed > 0:
            if rest.bit_count() < needed:
                return False
            scan = rest
            branch, maximum_degree = 0, 0
            while scan:
                bit = scan & -scan
                scan ^= bit
                index = bit.bit_length() - 1
                count = (neighbors[index] & rest).bit_count()
                if count > maximum_degree:
                    branch, maximum_degree = index, count
            if maximum_degree == 0:
                return True
            rest ^= 1 << branch
            if search(rest & ~neighbors[branch], needed - 1):
                return True
        return True

    return search((1 << len(bad)) - 1, target)


def audit(width):
    classes, options, incumbent, max_between = prepare()
    outside = blocker_sets(options, incumbent, max_between)
    report = Counter()
    for scope in combinations(classes.tolist(), width):
        free = frozenset(scope)
        candidates = [options[c][[x <= free for x in outside[c]]] for c in scope]
        owners = np.repeat(scope, [len(rows) for rows in candidates])
        rows = np.concatenate(candidates)
        source_size = sum(int(incumbent[c].sum()) for c in scope)
        assert len(rows) >= source_size
        if len(rows) == source_size:
            report["no_additional_candidates"] += 1
            continue
        mask = rows @ rows.T > 64 - max_between[owners[:, None], owners[None, :]]
        np.fill_diagonal(mask, False)
        if has_independent_set(mask, source_size + 1):
            raise RuntimeError(f"An improving search window was found: {scope}")
        report["excluded_by_branching"] += 1
    return {
        "native_classes": 135,
        "exceptional_classes_frozen": 2,
        "classes_per_window": width,
        "signed_native_candidates": sum(map(len, options)),
        "windows": sum(report.values()),
        **dict(report),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, choices=(1, 2), default=2)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(args.width), indent=2) + "\n"
    print(result, end="")
    if args.json:
        args.json.write_text(result)


if __name__ == "__main__":
    main()
