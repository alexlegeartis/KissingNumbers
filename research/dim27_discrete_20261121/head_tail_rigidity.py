"""Test length-eight tails in three head-body families after deleting the poles."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, product
import json

import numpy as np

from bw_tail_rigidity import body_profiles, spherical_vertices, tail_blocks, vertex_recipes
from source import COMPLEMENT, OCTAD, leech_shell, load_points


def orthogonal_half_sums(roots):
    """Enumerate each signed body once, even when it has many representations."""
    roots_f = roots.astype(np.float32)
    candidates = []
    for start in range(0, len(roots), 128):
        a, b = np.where(roots_f[start:start + 128] @ roots_f[start:].T == 0)
        allowed = a < b
        h = (roots[a[allowed] + start] + roots[b[allowed] + start]) // 2
        if len(h):
            candidates.append(np.unique(h, axis=0))
    return np.unique(np.concatenate(candidates), axis=0)


def sparse_bodies():
    """All signed four-coordinate bodies meeting both octad blocks."""
    points = []
    octad = set(OCTAD.tolist())
    signs = np.array(list(product((-4, 4), repeat=4)), dtype=np.int16)
    for support in combinations(range(24), 4):
        if 0 < len(octad.intersection(support)) < 4:
            rows = np.zeros((len(signs), 24), dtype=np.int16)
            rows[:, support] = signs
            points.append(rows)
    return np.concatenate(points)


def initialize(fixed):
    global BLOCKS
    BLOCKS = fixed


def sample(rows):
    return body_profiles(rows.astype(np.float32), BLOCKS)


def audit(family, workers):
    shell = 2 * leech_shell().astype(np.int16)
    points = load_points()
    heights = np.sum(points[:, 24:] ** 2, axis=1)
    fixed = points[(heights > 0) & (heights < 128)]
    if family == "bw":
        coordinates, zero = COMPLEMENT, OCTAD
        bodies = orthogonal_half_sums(shell[~np.any(shell[:, zero], axis=1)][:, coordinates])
    elif family == "e8":
        coordinates, zero = OCTAD, COMPLEMENT
        bodies = orthogonal_half_sums(shell[~np.any(shell[:, zero], axis=1)][:, coordinates])
    else:
        coordinates = np.arange(24)
        bodies = sparse_bodies()
    assert np.all(np.sum(bodies ** 2, axis=1) == 64)
    normals, blocks = tail_blocks(fixed, coordinates)
    assert sum(np.count_nonzero(row) == 1 for row in normals) == 3
    signed = np.r_[normals, -normals]
    recipes = vertex_recipes(signed)

    counts = Counter()
    batches = (bodies[i:i + 512] for i in range(0, len(bodies), 512))
    if workers == 1:
        initialize(blocks)
        for report in map(sample, batches):
            counts.update(report)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=initialize, initargs=(blocks,)) as pool:
            for report in pool.map(sample, batches):
                counts.update(report)

    axial_only = 0
    for constraint, number in counts.items():
        tails = spherical_vertices(signed, recipes, constraint, 64)
        if len(tails):
            assert np.all(np.sort(np.abs(tails), axis=1) == [0, 0, 8])
            axial_only += number
    return {"body_family": family, "body_count": len(bodies), "profiles": len(counts),
            "bodies_with_axial_tails": axial_only, "bodies_with_off_axis_tails": 0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=("bw", "e8", "mixed_sparse"))
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    print(json.dumps(audit(args.family, args.workers), indent=2))
