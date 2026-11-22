"""Count and optimize whole projected shells in the published rank-45 P48m section."""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import gcd, lcm, prod

import numpy as np
import sympy as sp


GRAM = np.array([[6, 2, 2], [2, 6, -2], [2, -2, 6]])
DET = 128


def integer_row(entries):
    values = list(map(Fraction, entries))
    denom = lcm(*(v.denominator for v in values))
    values = [int(v * denom) for v in values]
    divisor = gcd(*values) or 1
    return [v // divisor for v in values]


def unique_counts(rows, totals):
    # Independence modulo a prime certifies independence over the rationals.
    echelon = {}
    selected = []
    for i, row in enumerate(rows):
        residue = np.array([int(x) % 1_000_003 for x in row], dtype=np.int64)
        for col in sorted(echelon):
            residue = (residue - residue[col] * echelon[col]) % 1_000_003
        pivots = np.flatnonzero(residue)
        if len(pivots):
            col = pivots[0]
            echelon[col] = residue * pow(int(residue[col]), -1, 1_000_003) % 1_000_003
            selected.append(i)
        if len(selected) == len(row):
            break
    assert len(selected) == len(rows[0])
    matrix = sp.Matrix([list(rows[i]) for i in selected]).to_DM()
    rhs = sp.Matrix([totals[i] for i in selected]).to_DM()
    numerator, denominator = matrix.solve_den(rhs)
    counts = [Fraction(int(v), int(denominator)) for v in numerator.to_Matrix()]
    assert all(v.denominator == 1 and v > 0 for v in counts)
    counts = np.array(list(map(int, counts)), dtype=object)
    assert np.array_equal(np.array(rows, dtype=object) @ counts, totals)
    return counts


@lru_cache(None)
def pairings(indices):
    if not indices:
        return 1
    return sum(int(GRAM[indices[0], indices[i]]) * pairings(indices[1:i] + indices[i+1:])
               for i in range(1, len(indices)))


def spherical_moments(degrees, norm, population):
    totals = [Fraction(population * norm**(len(e)//2) * pairings(e),
                       prod(range(48, 48+len(e), 2))) for e in degrees]
    assert all(v.denominator == 1 for v in totals)
    return np.array(list(map(int, totals)), dtype=object)


def harmonic12(inner, norm, direction_norm):
    # Zonal harmonic in R^48, with leading coefficient one.
    inner = inner.astype(object)
    value = inner**12
    coefficient = Fraction(1)
    for j in range(1, 7):
        coefficient *= Fraction(-(14-2*j)*(13-2*j), 2*j*(70-2*j))
        value += coefficient * (norm*direction_norm)**j * inner**(12-2*j)
    return value


def largest_configuration(families, label_index, distance):
    families = sorted(families, key=lambda item: -item[2])
    adjacent = []
    for i, (a, r, _) in enumerate(families):
        neighbors = set()
        for j in range(i+1, len(families)):
            b, s, _ = families[j]
            difference = tuple((x-y) % DET for x, y in zip(a, b))
            twice_dot = r+s-distance[label_index[difference]]
            if twice_dot <= 0 or twice_dot**2 <= r*s:
                neighbors.add(j)
        adjacent.append(neighbors)
    best_weight, best_set = 0, []

    def visit(choices, used, weight):
        nonlocal best_weight, best_set
        if weight > best_weight:
            best_weight, best_set = weight, used
        remaining = sum(families[i][2] for i in choices)
        for k, i in enumerate(choices):
            if weight+remaining <= best_weight:
                break
            w = families[i][2]
            visit([j for j in choices[k+1:] if j in adjacent[i]], used+[i], weight+w)
            remaining -= w

    visit(list(range(len(families))), [], 0)
    return best_weight, [families[i] for i in best_set]


def main():
    adj = np.array(sp.Matrix(GRAM).adjugate(), dtype=np.int64)
    assert np.array_equal(GRAM @ adj, DET*np.eye(3, dtype=np.int64))
    # Cauchy bounds each individual inner product by sqrt(8*6).
    dual = np.array(list(product(range(-6, 7), repeat=3)))
    normal = ((dual @ adj)*dual).sum(axis=1)
    dual, normal = dual[normal <= 8*DET], normal[normal <= 8*DET]
    labels, ids = np.unique((dual @ adj) % DET, axis=0, return_inverse=True)
    lookup = {tuple(label): i for i, label in enumerate(labels)}
    assert len(labels) == DET
    negative = [lookup[tuple(-v % DET)] for v in labels]
    _, orbit = np.unique(np.minimum(negative, np.arange(DET)), return_inverse=True)
    shortest = np.array([normal[ids == i].min() for i in range(DET)])
    nearest = normal == shortest[ids]
    groups = [nearest & (orbit[ids] == j) for j in range(max(orbit)+1)]
    assert len(groups) == 68 and max(shortest) == 5*DET//2

    degrees = [e for d in range(0, 11, 2) for e in combinations_with_replacement(range(3), d)]
    powers = np.array([dual[:, e].prod(axis=1) for e in degrees])
    rows = np.array([powers[:, sel].sum(axis=1) for sel in groups]).T.tolist()
    # Published origin kissing number counts three translates at squared radius 23/4.
    cluster = np.array([[1, 0, 1], [0, -1, 1], [0, 0, 1]])
    cluster_orbits = [orbit[lookup[tuple(x)]] for x in (cluster @ adj) % DET]
    calibration = np.bincount(cluster_orbits, minlength=len(groups)).tolist()
    pure6 = (ids == 0) & (normal == 6*DET)
    rhs6 = spherical_moments(degrees, 6, 52_416_000) - powers @ pure6
    six = unique_counts(rows+[calibration], [*rhs6, 7_379_838])[orbit]

    q = sp.Symbol('q')
    e4 = 1 + sum(240*sp.divisor_sigma(i, 3)*q**i for i in range(1, 5))
    delta = q-24*q**2+252*q**3-1472*q**4
    theta8 = int(sp.series(e4**6-1440*e4**3*delta+125280*delta**2, q, 0, 5).coeff(q, 4))
    known8 = ((ids == 0) & (normal == 8*DET)).astype(int)
    known8 = known8 + six[ids]*(normal == shortest[ids]+2*DET)
    rhs8 = spherical_moments(degrees, 8, theta8) - powers @ known8
    direction = np.array([1, 2, 4])
    norm_y = int(direction @ GRAM @ direction)
    h6 = harmonic12(dual @ direction, 6, norm_y)
    h8 = harmonic12(dual @ direction, 8, norm_y)
    # The weighted theta series has weight 36 and begins at q^3, hence equals a multiple of Delta^3.
    total = -72*(h6 @ (six[ids]*nearest + pure6)) - h8 @ known8
    harmonic = integer_row([*(sum(h8[sel]) for sel in groups), total])
    eight = unique_counts(rows+[harmonic[:-1]], [*rhs8, harmonic[-1]])[orbit]

    calibration_label = lookup[tuple((np.array([-3, -3, -1]) @ adj) % DET)]
    assert (six[calibration_label], eight[calibration_label]) == (49, 941_720)
    families = [(tuple(map(int, label)), int(norm*DET-shortest[i]), int(weights[i]))
                for norm, weights in [(6, six), (8, eight)] for i, label in enumerate(labels)
                if norm*DET-shortest[i] <= 6*DET]
    maximum, chosen = largest_configuration(families, lookup, 6*DET-shortest)
    assert maximum == 7_392_584
    print(f"{len(labels)} residue classes, {len(families)} complete families")
    print("Calibrating body counts at normal (-3, -3, -1):", six[calibration_label], eight[calibration_label])
    for label, radius, count in chosen:
        print(label, "squared radius", Fraction(radius, DET), "population", count)
    print("Maximum under the zero-head class-distance test:", maximum)


if __name__ == '__main__':
    main()
