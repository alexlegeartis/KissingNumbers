"""Independent polynomial-coefficient replay of all signed-tail profiles.

No imports from the certificate generator. Coefficients of
(1+t(X+X^-1))^a (1+2t)^(k-a) enumerate products against a fixed
weight-a tail. Direct small-orbit enumeration checks the coefficient model.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, isqrt, prod
from pathlib import Path
import json

import numpy as np

HERE = Path(__file__).resolve().parent


def coefficients(k, left, right):
    states = {(0, 0): 1}
    for coordinate in range(k):
        result = Counter(states)
        for (used, dot), count in states.items():
            if used == right:
                continue
            if coordinate < left:
                result[used+1, dot+1] += count
                result[used+1, dot-1] += count
            else:
                result[used+1, dot] += 2*count
        states = result
    return {dot: count for (used, dot), count in states.items() if used == right}


def direct(k, a, b):
    result = Counter()
    for support in combinations(range(k), b):
        for signs in product((-1, 1), repeat=b):
            result[sum(s for i, s in zip(support, signs) if i < a)] += 1
    return dict(result)


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    report = json.loads((HERE/'signed_norm10_certificate.json').read_text())
    shell = {int(j): count for j, count in report['norm10_absolute_product_distribution'].items()}
    for moment in range(6):
        expected = Q(15281788354560 * prod(range(1, 2*moment, 2)) * 100**moment,
                     prod(range(72, 72+2*moment, 2)))
        assert sum(count*j**(2*moment) for j, count in shell.items()) == expected
    # Cross-check polynomial enumeration including empty intersections,
    # full overlap, signs, and the zero-product coefficient.
    for k in range(1, 6):
        for a in range(1, k+1):
            for b in range(1, k+1):
                assert coefficients(k, a, b) == direct(k, a, b)
    for dim, row in report['dimensions'].items():
        k = int(dim)-72
        weights = row['tail_support_weights']
        rates = row['priority_rates']
        sizes = [2**w*comb(k, w) for w in weights]
        assert sizes == row['orbit_populations']
        assert weights[0] == 1 and row['tail_rank'] == k
        scale = Q(row['head_scale_squared'])
        p, d = scale.numerator, scale.denominator
        residual = 8*d-10*p
        matrix = []
        for a in weights:
            entries = []
            for b, size in zip(weights, sizes):
                hist = coefficients(k, a, b)
                assert sum(hist.values()) == size
                entry = size
                for j in (4, 5, 6):
                    # dot integral: dot>sqrt(C) iff dot>floor(sqrt(C)).
                    threshold = isqrt(((4*d-j*p)**2*a*b)//residual**2)
                    count = sum(n for dot, n in hist.items() if dot > threshold)
                    entry += shell[j]//2*count
                entries.append(entry)
            matrix.append(entries)
        assert matrix == row['closed_neighborhood_matrix']
        closed = [sum(c*r for c, r in zip(entries, rates)) for entries in matrix]
        alpha = (15281788354560//2)*sum((Q(n*r, degree)
                    for n, r, degree in zip(sizes, rates, closed)), Q(0))
        assert alpha == Q(int(row['priority_expectation_numerator']),
                          int(row['priority_expectation_denominator']))
        pairs = -(-alpha.numerator//alpha.denominator)
        poles = []
        for i, j in combinations(range(k), 2):
            for s, t in product((-2, 2), repeat=2):
                v = np.zeros(k, dtype=np.int64)
                v[i], v[j] = s, t
                poles.append(v)
        poles = np.array(poles)
        assert len(poles) == 2*k*(k-1) == row['pure_tail_points']
        dots = poles@poles.T
        assert np.all(np.diag(dots) == 8)
        np.fill_diagonal(dots, 0)
        assert np.max(dots) == 4
        assert np.array_equal(poles.T@poles, 16*(k-1)*np.eye(k, dtype=np.int64))
        assert row['bound'] == 6218175600+2*pairs+len(poles)
        print('PASS independent polynomial profiles, poles and rational bound:', dim, row['bound'], flush=True)


if __name__ == '__main__':
    main()
