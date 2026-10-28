"""Exact norm-ten existence bounds from signed constant-weight tail orbits.

Every orbit consists of all vectors with w nonzero coordinates in {+1,-1}
in R^k. Hyperoctahedral symmetry makes each neighborhood profile constant
within an orbit. Integer rates and exact combinatorial profiles suffice;
the large tail sets and selected norm-ten heads need not be materialized.
"""
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json
import sys

from norm10_bound_certificate import theta_coefficients, solve_moments

HERE = Path(__file__).resolve().parent
SCALE = Q(8, 13)
SPECS = {
    10: ([1, 4], [19871466, 10000000], 997443),
    11: ([1, 4], [21748661, 10000000], 1416419),
    12: ([1, 4], [23745110, 10000000], 1924770),
    13: ([1, 4, 9], [26122963, 10000000, 3247], 2523388),
    15: ([1, 4, 9], [31936463, 10000000, 4877], 4009295),
}


def population(k, weight):
    return 2**weight * comb(k, weight)


def histogram(k, left_weight, right_weight):
    """Products against a fixed positive left tail, counting every right tail."""
    result = Counter()
    for overlap in range(max(0, left_weight+right_weight-k),
                         min(left_weight, right_weight)+1):
        placements = comb(left_weight, overlap) * comb(k-left_weight, right_weight-overlap)
        outside_signs = 2**(right_weight-overlap)
        for positive in range(overlap+1):
            result[2*positive-overlap] += placements * outside_signs * comb(overlap, positive)
    assert sum(result.values()) == population(k, right_weight)
    assert all(result[v] == result[-v] for v in result)
    return result


def certificate():
    theta, coefficients = theta_coefficients()
    n8, n10 = theta[4:6]
    assert (n8, n10) == (6218175600, 15281788354560)
    small = solve_moments(n8, 80, {})
    shell = solve_moments(n10, 100, {6: small[4], 10: 2})
    assert sum(shell.values()) == n10
    assert all(shell[j] % 2 == 0 for j in (4, 5, 6))
    p, d = SCALE.numerator, SCALE.denominator
    residual = 8*d-10*p
    assert Q(3, 5) < SCALE <= Q(16, 25)
    assert 3*SCALE+(8-10*SCALE) <= 4
    assert 10*SCALE-(8-10*SCALE) > 4
    rows = {}
    for k, (weights, rates, expected_pairs) in SPECS.items():
        assert weights[0] == 1 and all(1 <= w <= k for w in weights)
        assert len(weights) == len(rates) and all(r > 0 for r in rates)
        # The weight-one orbit spans R^k. All other tails and pure D_k
        # roots use the same k coordinates. Their norms are w and 8.
        sizes = [population(k, w) for w in weights]
        matrix = []
        profiles = []
        for a in weights:
            row = []
            raw = []
            for b, size in zip(weights, sizes):
                hist = histogram(k, a, b)
                counts = []
                for j in (4, 5, 6):
                    assert 4*d-j*p > 0
                    # Strict inequality: boundary products are compatible.
                    counts.append(sum(count for dot, count in hist.items()
                                      if dot > 0 and residual**2*dot**2 > (4*d-j*p)**2*a*b))
                raw.append(counts)
                row.append(size + sum(shell[j]//2*count for j, count in zip((4, 5, 6), counts)))
            matrix.append(row)
            profiles.append(raw)
        for i in range(len(weights)):
            for j in range(len(weights)):
                assert sizes[i]*matrix[i][j] == sizes[j]*matrix[j][i]
        closed = [sum(a*r for a, r in zip(row, rates)) for row in matrix]
        alpha = (n10//2) * sum((Q(size*rate, degree)
                               for size, rate, degree in zip(sizes, rates, closed)), Q(0))
        pairs = -(-alpha.numerator//alpha.denominator)
        assert pairs == expected_pairs
        poles = 2*k*(k-1)  # 2(+/-e_i +/-e_j), i<j, norm eight
        rows[str(72+k)] = {
            'bound': n8+2*pairs+poles,
            'equator': n8,
            'tail_rank': k,
            'head_scale_squared': str(SCALE),
            'tail_support_weights': weights,
            'tail_squared_norms': weights,
            'priority_rates': rates,
            'orbit_populations': sizes,
            'tail_vectors': sum(sizes),
            'positive_conflict_populations_4_5_6': profiles,
            'closed_neighborhood_matrix': matrix,
            'weighted_closed_degrees': closed,
            'priority_expectation_numerator': str(alpha.numerator),
            'priority_expectation_denominator': str(alpha.denominator),
            'guaranteed_antipodal_pairs': pairs,
            'pure_tail_points': poles,
        }
    return {
        'status': 'PASS',
        'certificate_type': 'Exact signed-orbit weighted-priority existence proof; no selected head list',
        'theta_coefficients': theta,
        'modular_form_coefficients': coefficients,
        'norm10_absolute_product_distribution': shell,
        'dimensions': rows,
        'literature_inputs': [
            'Gamma72 even unimodular, minimum 8 (Nebe)',
            'All shells are spherical 11-designs (Venkov)',
            'Classical modular-form description of even unimodular theta series',
        ],
    }


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    result = certificate()
    path = HERE/'signed_norm10_certificate.json'
    if '--check' in sys.argv:
        assert json.loads(path.read_text()) == json.loads(json.dumps(result)), 'stale certificate'
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
