"""Verify projected shell sizes as functions of four counts in a D4 section."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import lcm, prod
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
K = np.array([[6,-3,-3,-3], [-3,6,0,0], [-3,0,6,0], [-3,0,0,6]], dtype=np.int64)
DET = 324
PRIME = 1_000_003


@lru_cache(None)
def pairings(indices):
    if not indices:
        return 1
    return sum(int(K[indices[0], indices[i]]) * pairings(indices[1:i] + indices[i+1:])
               for i in range(1, len(indices)))


def shell_moment(norm, size, indices):
    half = len(indices)//2
    result = Fraction(size * norm**half * pairings(indices), prod(range(48, 48+2*half, 2)))
    assert result.denominator == 1
    return result.numerator


def add_pivot(row, pivots):
    v = np.array([int(x) % PRIME for x in row], dtype=np.int64)
    for j in sorted(pivots):
        v = (v - v[j]*pivots[j]) % PRIME
    nz = np.flatnonzero(v)
    if not len(nz):
        return False
    j = nz[0]
    pivots[j] = (v * pow(int(v[j]), -1, PRIME)) % PRIME
    return True


def zonal(inner, norm, direction_norm):
    coefficients = [Fraction(1)]
    for j in range(1, 7):
        coefficients.append(coefficients[-1] * Fraction(-(14-2*j)*(13-2*j), 2*j*(70-2*j)))
    scale = lcm(*(v.denominator for v in coefficients))
    x = inner.astype(object)
    return sum(int(scale*c) * (norm*direction_norm)**j * x**(12-2*j)
               for j, c in enumerate(coefficients))


def shell_sizes():
    q = sp.Symbol("q")
    e4 = 1 + sum(240*sp.divisor_sigma(i, 3)*q**i for i in range(1, 5))
    delta = q - 24*q**2 + 252*q**3 - 1472*q**4
    theta = sp.series(e4**6-1440*e4**3*delta+125280*delta**2, q, 0, 5).removeO()
    assert [theta.coeff(q,i) for i in (0,1,2,3)] == [1,0,0,52_416_000]
    return int(theta.coeff(q,3)), int(theta.coeff(q,4))


def verify_families():
    with np.load(HERE/"families.npz", allow_pickle=False) as data:
        assert all(data[k].dtype.kind in "iu" for k in data)
        template = data["affine_counts"].astype(object)
        selected, radii, calibration = data["labels"], data["radii"], data["calibrations"]
    adj = np.array(sp.Matrix(K).adjugate(), dtype=np.int64)
    assert np.array_equal(K @ adj, DET*np.eye(4, dtype=np.int64))
    # Cauchy's inequality bounds each inner product by sqrt(6*8).
    dual = np.array(list(product(range(-6,7), repeat=4)), dtype=np.int64)
    normal = ((dual @ adj)*dual).sum(axis=1)
    dual, normal = dual[normal <= 8*DET], normal[normal <= 8*DET]
    labels, ids = np.unique((dual @ adj) % DET, axis=0, return_inverse=True)
    lookup = {tuple(x): i for i, x in enumerate(labels)}
    assert len(labels) == DET and not np.any(labels[0])
    opposite = [lookup[tuple(-x % DET)] for x in labels]
    _, orbit = np.unique(np.minimum(opposite, np.arange(DET)), return_inverse=True)
    size = max(orbit)+1
    assert size == 164 and template.shape == (2*size, 5)
    six, eight = template[:size], template[size:]
    shortest = np.array([min(normal[ids == i]) for i in range(DET)])
    assert max(shortest) == 3*DET
    nearest = normal == shortest[ids]
    next_nearest = normal == shortest[ids]+2*DET
    pure6 = (ids == 0) & (normal == 6*DET)
    pure8 = (ids == 0) & (normal == 8*DET)
    indices = orbit[ids]
    first = [np.flatnonzero(nearest & (indices == j)) for j in range(size)]
    second = [np.flatnonzero(next_nearest & (indices == j)) for j in range(size)]

    def group(values, slices):
        return np.array([values[s].sum() for s in slices], dtype=object)

    moments = [e for d in range(0, 11, 2) for e in combinations_with_replacement(range(4), d)]
    size6, size8 = shell_sizes()
    pivots = {}
    for e in moments:
        values = dual[:, e].prod(axis=1)
        a, b = group(values, first), group(values, second)
        add_pivot(a, pivots)
        rhs6 = shell_moment(6, size6, e) - int(values[pure6].sum())
        rhs8 = shell_moment(8, size8, e) - int(values[pure8].sum())
        assert (a @ six).tolist() == [rhs6,0,0,0,0]
        assert (a @ eight + b @ six).tolist() == [rhs8,0,0,0,0]
    assert len(pivots) == 160

    # Fix four entries at norm six. The ordinary equations then determine every entry.
    calibrated = pivots.copy()
    for j, v in enumerate(calibration):
        k = lookup[tuple(v @ adj % DET)]
        assert v @ adj @ v == shortest[k]
        unit = np.eye(size, dtype=np.int64)[orbit[k]]
        assert add_pivot(unit, calibrated)
        assert six[orbit[k]].tolist() == [0, *(int(i==j) for i in range(4))]
    assert len(calibrated) == size

    # A degree-12 harmonic theta series is a multiple of Delta^3. Its
    # coefficient at q^4 is -72 times its coefficient at q^3.
    harmonic = pivots.copy()
    for direction in product(range(1,5), repeat=4):
        t = np.array(direction)
        inner = dual @ t
        norm_y = int(t @ K @ t)
        h8 = zonal(inner, 8, norm_y)
        row8 = group(h8, first)
        if not add_pivot(row8, harmonic):
            continue
        h6 = zonal(inner, 6, norm_y)
        row6 = 72*group(h6, first) + group(h8, second)
        rhs = -72*h6[pure6].sum() - h8[pure8].sum()
        assert (row6 @ six + row8 @ eight).tolist() == [rhs,0,0,0,0]
        if len(harmonic) == size:
            break
    assert len(harmonic) == size

    # All distances and squared radii are scaled by det K in this check.
    family_forms = []
    total = np.zeros(5, dtype=object)
    assert len(selected) == len(radii) == 17
    assert len({(tuple(c), int(r)) for c,r in zip(selected,radii)}) == len(radii)
    for i, (c, r) in enumerate(zip(selected, radii)):
        k = lookup[tuple(c)]
        assert 0 < r <= 6*DET
        norm = int(r + shortest[k])
        assert norm in (6*DET, 8*DET)
        family_forms.append((six if norm == 6*DET else eight)[orbit[k]])
        total += family_forms[-1]
        for other, s in zip(selected[:i], radii[:i]):
            difference = lookup[tuple((c-other) % DET)]
            doubled_product = r+s-(6*DET-shortest[difference])
            assert doubled_product <= 0 or doubled_product**2 <= r*s
    assert total.tolist() == [4_940_961, -70, -70, -70, 0]
    return calibration, np.array(family_forms, dtype=np.int64)


if __name__ == "__main__":
    _, forms = verify_families()
    print("Population as a function of the four calibrations:", forms.sum(axis=0))
