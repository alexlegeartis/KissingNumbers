"""Exclude omitting up to six incumbent E8 head lines with a pure midpoint refill."""
from collections import Counter, defaultdict
from itertools import combinations, product
import json

import numpy as np

from e8_input import midpoint_options, positive
from source import COMPLEMENT, OCTAD, RESIDUAL, leech_shell, load_points


def head_blockers(options, heads):
    cost = np.abs(options[:, :8] @ heads[:, :8].T) + np.abs(options[:, 8:] @ heads[:, 8:].T)
    return (cost > 64).astype(np.int64) @ (1 << np.arange(len(heads), dtype=np.int64))


def midpoint_conflicts(points):
    return np.abs(points[:, :8] @ points[:, :8].T) + np.abs(points[:, 8:] @ points[:, 8:].T) > 64


def root_requirements(options, parents, blocked_by, max_omissions):
    required = defaultdict(int)
    for root in range(120):
        incident = np.flatnonzero((parents == root).any(1))
        compatible = ~midpoint_conflicts(options[incident])
        first, second = np.where(np.triu(compatible, 1))
        signatures = {int(blocked_by[i] | blocked_by[j]) for i, j in zip(incident[first], incident[second])}
        minimal = []
        for mask in sorted(signatures, key=lambda x: (x.bit_count(), x)):
            if mask.bit_count() <= max_omissions and not any(mask & smaller == smaller for smaller in minimal):
                minimal.append(mask)
                required[mask] |= 1 << root
    return required


def root_pair_bound(choices, parents):
    conflicts = midpoint_conflicts(choices)
    np.fill_diagonal(conflicts, False)  # Two root options may use the same midpoint.
    states_by_root = []
    for root in range(120):
        ids = np.flatnonzero((parents == root).any(1))
        first, second = np.where(np.triu(~conflicts[np.ix_(ids, ids)], 1))
        if len(first):
            states_by_root.append(np.c_[ids[first], ids[second]])

    states = np.concatenate(states_by_root)
    offsets = np.cumsum([0, *map(len, states_by_root[:-1])])
    bad = np.zeros((len(states), len(states)), dtype=bool)
    for i, j in product(range(2), repeat=2):
        bad |= conflicts[states[:, i, None], states[None, :, j]]
    root_bad = np.logical_and.reduceat(bad, offsets, axis=0)
    root_bad = np.logical_and.reduceat(root_bad, offsets, axis=1)
    np.fill_diagonal(root_bad, False)

    matched = np.zeros(len(root_bad), dtype=bool)
    matching_size = 0
    for u in np.argsort(root_bad.sum(1), kind="stable"):
        if not matched[u]:
            v = np.flatnonzero(root_bad[u] & ~matched)
            if len(v):
                matched[[u, v[0]]] = True
                matching_size += 1
    # Each disjoint incompatible pair forces at least one missing degree-two root.
    return len(states_by_root) - matching_size


def audit():
    shell = 2 * leech_shell().astype(np.int16)
    source = load_points()
    outside = np.any(source[:, COMPLEMENT], axis=1)
    exterior = np.unique(source[outside][:, RESIDUAL], axis=0)
    kernel = source[~outside][:, RESIDUAL]
    roots = shell[~np.any(shell[:, COMPLEMENT], axis=1)][:, OCTAD]
    points, parents = midpoint_options(roots, exterior)

    heads = kernel[np.sum(kernel[:, 8:] ** 2, axis=1) == 64]
    heads = heads[positive(heads[:, :8]) & positive(heads[:, 8:])]
    assert len(heads) == 24
    blocked_by = head_blockers(points, heads)
    requirements = root_requirements(points, parents, blocked_by, 6)
    counts = []
    for size in range(7):
        report = Counter()
        maximum = 0
        for removed in combinations(range(len(heads)), size):
            mask = sum(1 << index for index in removed)
            current = mask
            reachable = requirements.get(current, 0)
            while current:
                current = (current - 1) & mask
                reachable |= requirements.get(current, 0)
            first_bound = reachable.bit_count()
            maximum = max(maximum, first_bound)
            target = 65 + 2 * size
            if first_bound < target:
                report["excluded_by_incidence"] += 1
            else:
                keep = (blocked_by & ~mask) == 0
                bound = root_pair_bound(points[keep], parents[keep])
                if bound < target:
                    report["excluded_by_matching"] += 1
                else:
                    report["remaining"] += 1
        counts.append({
            "omitted_head_lines": size, "degree_two_target": 65 + 2 * size,
            "maximum_first_bound": maximum, **dict(report),
        })
    return counts


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
