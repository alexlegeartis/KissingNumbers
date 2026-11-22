"""Test real two-for-three rotations of opposite residuals in every native fibre."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
from multiprocessing import get_context

from numba import njit
import numpy as np

from fibres import native_fibres
from source import COMPLEMENT, RESIDUAL, leech_shell, load_points


def vertex_recipes(normals):
    triples, adjugates, determinants = [], [], []
    for triple in combinations(range(len(normals)), 3):
        a, b, c = normals[list(triple)]
        adjugate = np.stack([np.cross(b, c), np.cross(c, a), np.cross(a, b)], axis=1)
        determinant = a @ adjugate[:, 0]
        if determinant:
            triples.append(triple)
            adjugates.append(np.sign(determinant) * adjugate)
            determinants.append(abs(determinant))
    return np.array(triples), np.array(adjugates), np.array(determinants)


@njit
def support_limits(slacks, slopes, tail_types, type_count):
    size = (type_count, 1 + np.max(slacks))
    largest = np.full(size, -999, np.int32)
    smallest = np.full(size, 999, np.int32)
    removed = 0
    for i, slack in enumerate(slacks):
        slope = slopes[i]
        if slack < abs(slope):
            removed += 1
            continue
        # ||tail||^2 <= 128. Cauchy bounds its entire support for |s| <= 1.
        if slack * slack >= slope * slope + 8192:
            continue
        tail = tail_types[i]
        largest[tail, slack] = max(largest[tail, slack], slope)
        smallest[tail, slack] = min(smallest[tail, slack], slope)
    return largest, smallest, removed


@njit
def reaches_radius(caps, required_squared, normals, triples, adjugates, determinants):
    for i, triple in enumerate(triples):
        numer = np.zeros(3, np.int64)
        for j in range(3):
            numer += caps[triple[j]] * adjugates[i, :, j]
        det = determinants[i]
        if np.sum(numer * numer) < required_squared * det * det:
            continue
        for j, normal in enumerate(normals):
            if np.sum(normal * numer) > caps[j] * det:
                break
        else:
            return True
    return False


@njit
def excludes_belt(slopes, tail_indices, slacks, factors, labels, normals,
                  triples, adjugates, determinants, orthogonal):
    denominator = 256
    intervals = [(0, 182)]  # Deliberately contains [0, 1/sqrt(2)].
    tests = 0
    while intervals:
        left, right = intervals.pop()
        tests += 1
        # Vertices below live in U = 8 * denominator * u.
        caps = 64 * denominator * np.sum(np.abs(normals), axis=1)
        for axis in orthogonal:
            caps[axis] = 0
        for i, slack in enumerate(slacks):
            tail = tail_indices[i]
            slope = slopes[tail, slack]
            cap = factors[tail] * (denominator * slack - min(slope * left, slope * right))
            label = labels[tail]
            caps[label] = min(caps[label], cap)
        radius = 4096 * (denominator ** 2 - right ** 2)
        if not reaches_radius(caps, radius, normals, triples, adjugates, determinants):
            continue
        if right - left == 1:
            return False, tests
        middle = (left + right) // 2
        intervals.extend(((left, middle), (middle, right)))
    return True, tests


def initialize():
    global projections, fibres, fixed_body, fixed_residual, tail_types
    global normals, labels, factors, recipes
    source = load_points().astype(np.int16)
    projections, _ = native_fibres(2 * leech_shell().astype(np.int16))
    projections = projections.reshape(-1, len(COMPLEMENT))
    locations = {tuple(w): i for i, w in enumerate(projections)}
    fibres = [set() for _ in projections]
    for p in source[np.sum(source[:, COMPLEMENT] ** 2, axis=1) == 64]:
        fibres[locations[tuple(p[COMPLEMENT])]].add(tuple(p[RESIDUAL]))

    fixed = source[np.any(source[:, 24:] != 0, axis=1)]
    assert np.max(np.sum(fixed[:, 24:] ** 2, axis=1)) == 128
    tails, tail_types = np.unique(fixed[:, 24:], axis=0, return_inverse=True)
    divisors = np.gcd.reduce(np.abs(tails), axis=1)
    assert np.all(8 % divisors == 0)
    normals, labels = np.unique(tails // divisors[:, None], axis=0, return_inverse=True)
    normals = normals.astype(np.int64)
    factors = 8 // divisors
    recipes = vertex_recipes(normals)
    # All products and partial sums are exactly representable in float32.
    fixed_body = fixed[:, COMPLEMENT].astype(np.float32).T
    fixed_residual = fixed[:, RESIDUAL].astype(np.float32).T


def check_fibre(index):
    incumbent = fibres[index]
    pairs = [row for row in incumbent if row[np.argmax(np.array(row) != 0)] > 0
             and tuple(-np.array(row)) in incumbent]
    assert pairs
    slacks = (64 - projections[index].astype(np.float32) @ fixed_body).astype(np.int32)
    slopes = (np.array(pairs, dtype=np.float32) @ fixed_residual).astype(np.int32)
    excluded = tests = max_tests = 0
    for row, coefficients in zip(pairs, slopes):
        tail = np.array(row[-3:])
        largest, smallest, removed = support_limits(slacks, coefficients, tail_types, len(labels))
        assert removed == (2 if np.any(tail) else 0)
        tail_indices, used_slacks = np.where(largest != -999)
        orthogonal = np.empty(0, dtype=np.int64)
        if np.any(tail):
            primitive = tail // np.gcd.reduce(np.abs(tail))
            orthogonal = np.flatnonzero(np.all(normals == primitive, axis=1)
                                        | np.all(normals == -primitive, axis=1))
            assert len(orthogonal) == 2
        for extrema in (largest, -smallest):
            ok, count = excludes_belt(extrema, tail_indices, used_slacks, factors, labels,
                                      normals, *recipes, orthogonal)
            excluded += ok
            tests += count
            max_tests = max(max_tests, count)
    return len(pairs), excluded, tests, max_tests


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    initialize()
    first = check_fibre(0)
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        results = np.array([first, *pool.map(check_fibre, range(1, len(fibres)), chunksize=4)])
    pairs, excluded, tests = map(int, results[:, :3].sum(axis=0))
    print(json.dumps(dict(fibres=len(fibres), opposite_pairs=pairs, excluded_signs=excluded,
                          unresolved_signs=2*pairs-excluded, slabs=tests,
                          most_slabs_per_sign=int(results[:, 3].max())), indent=2))
