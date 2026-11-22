"""Exhaust arbitrary real length-sqrt(32) tails of the pure BW midpoint refill."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json

import numpy as np

from bw_refill import contact_blocks, tilted_section
from source import COMPLEMENT, keys, leech_shell


def initialize(shell, frozen):
    global ROOTS, NORMALS, BODY_BLOCKS
    ROOTS = shell[:, COMPLEMENT]
    NORMALS, BODY_BLOCKS = tail_blocks(frozen, COMPLEMENT)


def tail_blocks(frozen, coordinates):
    palette = np.unique(frozen[:, 24:], axis=0)
    palette = [t for t in palette if any(t) and next(v for v in t if v) > 0]
    normals = np.array(palette, dtype=np.int64)
    blocks = []
    for tail in normals:
        plus = frozen[(frozen[:, 24:] == tail).all(1)][:, coordinates]
        minus = frozen[(frozen[:, 24:] == -tail).all(1)][:, coordinates]
        assert np.array_equal(np.sort(keys(plus)), np.sort(keys(minus)))
        blocks.append(np.unique(plus, axis=0).astype(np.float32).T)
    return normals, blocks


def profiles(edges):
    bodies = ((ROOTS[edges[:, 0]] + ROOTS[edges[:, 1]]) // 2).astype(np.float32)
    return body_profiles(bodies, BODY_BLOCKS)


def body_profiles(bodies, blocks):
    rhs = np.empty((len(bodies), len(blocks)), dtype=np.int16)
    remaining = np.arange(len(bodies))
    for i, fixed in enumerate(blocks):
        if not len(remaining):
            break
        rim = 64 - (bodies[remaining] @ fixed).max(axis=1).astype(np.int16)
        rhs[remaining, i] = rim
        remaining = remaining[rim >= 0]
    return Counter(tuple(row) for row in rhs[remaining].tolist())


def vertex_recipes(normals):
    """Return integer numerators and positive determinants for every invertible triple."""
    recipes = []
    for indices in combinations(range(len(normals)), 3):
        a, b, c = normals[list(indices)]
        adjugate = np.stack([np.cross(b, c), np.cross(c, a), np.cross(a, b)], axis=1)
        determinant = int(a @ adjugate[:, 0])
        if determinant:
            recipes.append((indices, adjugate * (1 if determinant > 0 else -1), abs(determinant)))
    return tuple(np.array(values) for values in zip(*recipes))


def spherical_vertices(normal, recipes, bound, squared_norm=32):
    choices, numerators, denominators = recipes
    rhs = np.tile(bound, 2).astype(np.int64)
    points = (numerators @ rhs[choices, None])[:, :, 0]
    length = (points * points).sum(axis=1)
    possible = length >= squared_norm * denominators ** 2
    p, d = points[possible], denominators[possible]
    keep = (p @ normal.T <= d[:, None] * rhs[None]).all(1)
    p, d = p[keep], d[keep]
    # If all vertices have norm at most the target radius, strict convexity of squared
    # norm implies that every point on the sphere is one of these vertices.
    if np.any((p * p).sum(axis=1) > squared_norm * d ** 2):
        raise RuntimeError(f"The polytope reaches beyond the sphere for profile {bound}")
    assert np.all(p % d[:, None] == 0)
    return np.unique(p // d[:, None], axis=0)


def audit(extra_class, workers):
    shell = 2 * leech_shell().astype(np.int16)
    frozen, _, _ = tilted_section(shell, extra_class)
    initialize(shell, frozen)
    signed = np.r_[NORMALS, -NORMALS]
    # The three pairs of coordinate-axis normals make all polytopes bounded.
    assert sum(np.count_nonzero(row) == 1 for row in NORMALS) == 3
    recipes = vertex_recipes(signed)
    distinct = Counter()
    contacts = 0

    def blocks():
        nonlocal contacts
        for edges in contact_blocks(shell):
            contacts += len(edges)
            yield edges

    if workers == 1:
        for counts in map(profiles, blocks()):
            distinct.update(counts)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=initialize,
                                 initargs=(shell, frozen)) as pool:
            for counts in pool.map(profiles, blocks()):
                distinct.update(counts)

    legal_bodies = legal_profiles = signed_points = 0
    lines = set()
    for profile, multiplicity in distinct.items():
        vertices = spherical_vertices(signed, recipes, profile)
        if len(vertices):
            assert np.all(np.sort(np.abs(vertices), axis=1) == [0, 4, 4])
            legal_profiles += 1
            legal_bodies += multiplicity
            signed_points += multiplicity * len(vertices)
            lines.update(tuple(map(int, t)) for t in vertices if next(v for v in t if v) > 0)

    assert contacts == 604800
    assert signed_points % 2 == 0
    return {
        "extra_fibre_class": extra_class, "contact_pairs": contacts,
        "constraint_profiles": len(distinct), "profiles_reaching_the_sphere": legal_profiles,
        "bodies_reaching_the_sphere": legal_bodies, "admissible_tail_pairs": signed_points // 2,
        "admissible_tail_lines": sorted(lines), "additional_real_tails": 0,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extra-class", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    print(json.dumps(audit(args.extra_class, args.workers), indent=2))
