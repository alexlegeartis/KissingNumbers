#!/usr/bin/env python3
"""Exact verification of a binary-code compression from length 62 to dimension 61.

For a distance-16 binary code, split a word as (i,j,s), where s has 60 signs.
Map it to (s/8, (-1)^i/4) if j=0, or
(sqrt(17)*s/32, (-1)^i/16) if j=1. All codewords survive.

The [62,26,16] code follows the construction in Grassl's code tables:
https://www.codetables.de/BKLC/BKLC.php?q=2&n=62&k=26
The companion binary_62_code module rebuilds it and checks all 2^26 words.
This file verifies the map for EVERY permitted distance and fibre pair,
so the geometric theorem applies to nonlinear codes too.

The added 4,194,304 weight-15 points are constructed here using a quaternary
[15,6,8] support code and an even binary Hamming [15,10,4] sign code.
All 4*binom(61,2) normalized weight-two points are also compatible.
No floating point calculation is used to accept an inequality.

GENERAL FAMILY. For any integer d>=4, put m=4d-4 and t=d/(4d-2).
Use tail magnitudes 1/sqrt(d), t/sqrt(d), and head coefficients
1/(2*sqrt(d)), sqrt((d-t*t)/(4*d*(d-1))). This embeds every binary
distance-d code of length 4d-2 into a 60-degree code in dimension 4d-3.
The two mixed-fibre inequalities reduce to
    (d-1)*R + 2*t <= d,    d*R - 2*t <= d,
where R=sqrt((d-t*t)/(d-1)). Both follow from R<=2d/(2d-1), since
    (2d/(2d-1))^2 - R^2 = d*(d-4)/(4*(d-1)*(2d-1)^2) >= 0.
The d=16 map at the top uses slightly simpler coefficients instead.

Run: python research/dim49plus/latitude_compression.py
"""

import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction as Q
from math import comb, sqrt
from pathlib import Path

import numpy as np

from binary_62_code import build_code, popcounts64, weight_enumerator


def quad_leq(p, q, r, bound=Q(1, 2)):
    """Exact decision p + q*sqrt(r) <= bound, for rational p,q and r>=0."""
    a, b = Q(p) - Q(bound), Q(q)
    if not b or not r:
        return a <= 0
    if b > 0:
        return a <= 0 and b * b * r <= a * a
    return a <= 0 or a * a <= b * b * r


def geometry(distance=16):
    """Check all 16 ordered pairs of labels and every allowed head distance."""
    m = 60
    tail = [Q(1, 4), Q(1, 16)]
    head2 = [Q(1, 64), Q(17, 1024)]
    checks, rows, bad = 0, [], []
    for j in range(2):
        assert m * head2[j] + tail[j] ** 2 == 1
    assert head2[0] != head2[1] and all(tail)
    for i, j in itertools.product(range(2), repeat=2):
        for k, l in itertools.product(range(2), repeat=2):
            b = (i != k) + (j != l)
            low = max(0, distance - b)
            for dh in range(low, m + 1):
                sign_dot = m - 2 * dh
                p = (-1 if i != k else 1) * tail[j] * tail[l]
                if j == l:
                    p += sign_dot * head2[j]
                    q = Q(0)
                else:
                    q = Q(sign_dot, 256)
                valid = quad_leq(p, q, 17)
                checks += 1
                if not valid:
                    bad.append([i, j, k, l, dh, str(p), str(q)])
                if dh == low:
                    rows.append({"first": [i, j], "second": [k, l],
                                 "removed_distance": b, "head_distance_min": low,
                                 "worst_ip": {"rational": str(p), "sqrt17": str(q)},
                                 "worst_ip_decimal": float(p) + float(q) * sqrt(17),
                                 "valid": valid})
    return {"valid": not bad, "ordered_distance_checks": checks,
            "worst_cases": rows, "violations": bad}


def general_family(d):
    """Exact parameter check plus the algebraic identity proving the family."""
    assert d >= 4
    t = Q(d, 4 * d - 2)
    r2 = (d - t * t) / (d - 1)
    rmax = Q(2 * d, 2 * d - 1)
    gap = Q(d * (d - 4), 4 * (d - 1) * (2 * d - 1) ** 2)
    assert rmax * rmax - r2 == gap >= 0
    assert (d - 1) * rmax + 2 * t == d
    assert d * rmax - 2 * t == d
    assert quad_leq(2 * t, d - 1, r2, d)
    assert quad_leq(-2 * t, d, r2, d)
    return {'distance': d, 'binary_length': 4 * d - 2,
            'spherical_dimension': 4 * d - 3, 't': str(t),
            'squared_ratio_gap': str(gap), 'valid': True}


def four_fibre_obstruction():
    """Sharp obstruction for the universal positive-head, one-tail template.

Try length 4d-1 -> dimension 4d-2, removing two bits as before. A fibre's
head coefficient h is at least 1/sqrt(4d), by its within-fibre constraint.
Its tail is consequently at most sqrt(3/(4d)) in magnitude. The two
Hamming-distance-two diagonals of the four labels must have opposite
tail signs, with magnitude products >=1/(4d). Hence precisely two
tails are positive and two negative. The two same-sign pairs have
label distance one, forcing their magnitude products <=1/(4d).
Multiplying gives equality throughout and forces all head coefficients
to equal 1/sqrt(4d). The same-sign products are then 3/(4d), a contradiction.
This excludes this universal template, not arbitrary codes or embeddings.
"""
    return {'excluded': 'universal four-fibre compression 4d-1 -> 4d-2',
            'scope': 'positive scalar head coefficient and one tail coordinate per label',
            'reason': 'distance-two products >=1/(4d) and same-sign distance-one '
                      'products <=1/(4d) force equality, contradicting the norms',
            'does_not_exclude': 'code-specific maps, additional fibres or tail coordinates'}


def gf4_support_code():
    """Build the cyclic [15,6,8]_4 code via its roots in GF(16)."""
    exp, log, a = [0] * 30, [0] * 16, 1
    for i in range(15):
        exp[i], log[a] = a, i
        a <<= 1
        if a & 16:
            a ^= 0b10011
    exp[15:] = exp[:15]

    def mul16(x, y):
        return 0 if not x or not y else exp[log[x] + log[y]]

    embedded = {0: 0, 1: 1, exp[5]: 2, exp[10]: 3}
    roots = {0, 1, 2, 3, 4, 5, 8, 10, 12}
    assert all((4 * r) % 15 in roots for r in roots)
    g = [1]
    for r in sorted(roots):
        ng = [0] * (len(g) + 1)
        for i, x in enumerate(g):
            ng[i] ^= mul16(x, exp[r])
            ng[i + 1] ^= x
        g = ng
    assert all(x in embedded for x in g)
    g = [embedded[x] for x in g]
    mul = np.array([[0, 0, 0, 0], [0, 1, 2, 3],
                    [0, 2, 3, 1], [0, 3, 1, 2]], dtype=np.uint8)
    messages = np.array(list(itertools.product(range(4), repeat=6)), dtype=np.uint8)
    words = np.zeros((4096, 15), dtype=np.uint8)
    for j in range(6):
        for i, coefficient in enumerate(g):
            words[:, i + j] ^= mul[messages[:, j], coefficient]
    weights = np.count_nonzero(words, axis=1)
    assert len(np.unique(words, axis=0)) == 4096
    assert weights[0] == 0 and weights[1:].min() == 8
    cells = 4 * np.arange(15)[None, :] + words
    supports = [sum(1 << int(x) for x in row) for row in cells]
    assert len(set(supports)) == 4096
    assert all(bin(x).count('1') == 15 for x in supports)
    # Exact pair overlap enumeration also checks the grid interpretation.
    worst, pairs = 0, 0
    for i, x in enumerate(supports):
        for y in supports[:i]:
            worst = max(worst, bin(x & y).count('1'))
            pairs += 1
    assert worst == 7
    return {"q": 4, "n": 15, "k": 6, "distance": 8,
            "weight_histogram": dict(Counter(map(int, weights))),
            "generator_polynomial_low_to_high": g,
            "supports": len(supports), "ambient_coordinates": 60,
            "maximum_overlap": worst, "pairs_checked": pairs}


def binary_sign_code():
    """Enumerate the even-weight subcode of the length-15 Hamming code."""
    words = []
    for x in range(1 << 15):
        if bin(x).count('1') % 2:
            continue
        syndrome = 0
        for i in range(15):
            if (x >> i) & 1:
                syndrome ^= i + 1
        if syndrome == 0:
            words.append(x)
    weights = [bin(x).count('1') for x in words]
    assert len(words) == 1024 and min(w for w in weights if w) == 4
    # The parity and syndrome equations are linear, so nonzero weights cover
    # the distances of ALL distinct pairs, not merely distances to zero.
    return {"q": 2, "n": 15, "k": 10, "distance": 4,
            "words": len(words), "weight_histogram": dict(Counter(weights))}


def sparse_geometry():
    checks = {}

    def record(name, value):
        checks[name] = bool(value)
        assert value, name

    record("weight15_same_support", Q(15 - 2 * 4, 15) <= Q(1, 2))
    record("weight15_distinct_supports", Q(7, 15) <= Q(1, 2))
    record("weight2_same_support", Q(2 - 2, 2) <= Q(1, 2))
    record("weight2_distinct_supports", Q(1, 2) <= Q(1, 2))
    record("weight15_vs_weight2_squared", Q(2, 15) <= Q(1, 4))
    record("dense_j0_vs_weight15_squared", Q(15, 64) <= Q(1, 4))
    record("dense_j1_vs_weight15_squared", Q(255, 1024) <= Q(1, 4))
    record("dense_j0_vs_weight2_two_head_squared", Q(2, 64) <= Q(1, 4))
    record("dense_j1_vs_weight2_two_head_squared", Q(34, 1024) <= Q(1, 4))
    record("dense_j0_vs_weight2_head_tail_squared", Q(9, 128) <= Q(1, 4))
    record("dense_j1_vs_weight2_head_tail_squared",
           quad_leq(Q(21, 2048), Q(4, 2048), 17, Q(1, 4)))
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    # Check the packed popcount on every 16-bit word and on patterns spanning
    # all eight bytes, independently against Python's binary representation.
    trials = np.concatenate([np.arange(1 << 16, dtype=np.uint64),
                             np.array([1 << i for i in range(64)] +
                                      [0xFFFFFFFFFFFFFFFF, 0xAAAAAAAAAAAAAAAA,
                                       0x5555555555555555, 0x0123456789ABCDEF], dtype=np.uint64)])
    assert popcounts64(trials).tolist() == [bin(int(x)).count('1') for x in trials]
    G = build_code()
    histogram = weight_enumerator(G)
    assert histogram[0] == 1 and histogram[1:16].sum() == 0 and histogram[16] > 0
    assert histogram.sum() == 1 << 26
    print('Binary [62,26,16]: all 67,108,864 words enumerated.', flush=True)
    dense = geometry()
    assert dense['valid']
    negative = geometry(15)
    assert not negative['valid']
    print('Dense map: all %d ordered fibre/distance cases pass exactly.' %
          dense['ordered_distance_checks'], flush=True)
    support = gf4_support_code()
    signs = binary_sign_code()
    sparse = sparse_geometry()
    total = (1 << 26) + support['supports'] * signs['words'] + 4 * comb(61, 2)
    previous = 60220542
    report = {
        'dimension': 61,
        'status': 'Explicit construction; exact geometry and exhaustive binary-code verification',
        'map': {'j0': '(s/8, (-1)^i/4)',
                'j1': '(sqrt(17)*s/32, (-1)^i/16)', 'head_coordinates': 60},
        'binary_code_input': {
            'parameters': [62, 26, 16], 'status': 'built and exhaustively verified',
            'words_enumerated': int(histogram.sum()),
            'weight_histogram': {str(i): int(x) for i, x in enumerate(histogram) if x},
            'source_chain': 'https://www.codetables.de/BKLC/BKLC.php?q=2&n=62&k=26'},
        'dense_geometry': dense,
        'negative_control_distance15': {
            'rejected': not negative['valid'],
            'violations': len(negative['violations'])},
        'support_code': support, 'sign_code': signs, 'sparse_geometry': sparse,
        'counts': {'dense': 1 << 26, 'weight15': support['supports'] * signs['words'],
                   'weight2': 4 * comb(61, 2), 'total': total,
                   'previous_research_bound': previous, 'gain': total - previous,
                   'relative_gain': total / previous - 1},
        'arithmetic': 'integers and fractions; exact quadratic-radical comparisons',
        'general_family': {
            'theorem': 'K(4d-3) >= A_binary(4d-2,d) for every integer d>=4',
            'proof': 'See GENERAL FAMILY in the module docstring; all inequalities '
                     'follow from the displayed nonnegative rational identity.',
            'parameter_checks': [general_family(d) for d in [4, 5, 8, 16, 32, 64, 128]]},
        'nearby_obstruction': four_fibre_obstruction(),
        'injectivity': 'different j have different nonzero head magnitudes; '
                       'different i have opposite nonzero last coordinate; '
                       'different s have different head sign vectors',
        'novelty': 'not established by this verifier'
    }
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report['counts'], indent=2))
    print('Written:', args.output)


if __name__ == '__main__':
    main()
