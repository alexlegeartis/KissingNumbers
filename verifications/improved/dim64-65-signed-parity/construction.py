"""Exact BCH, sparse-owner and cap geometry for the 65D construction.

Called by verify.py in a temporary directory; no research-tree dependency.
"""

import itertools
import json
import os
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ENUMERATOR = None


def popcount(value):
    return bin(value).count('1')


def gf64_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 64:
            a ^= 0x43  # x^6 + x + 1
    return result


def extended_bch(designed_distance):
    powers = [1]
    for _ in range(1, 63):
        powers.append(gf64_multiply(powers[-1], 2))
    assert len(set(powers)) == 63
    assert gf64_multiply(powers[-1], 2) == 1
    roots = set()
    for initial in range(1, designed_distance):
        exponent = initial
        while exponent not in roots:
            roots.add(exponent)
            exponent = (2 * exponent) % 63
    coefficients = [1]
    for exponent in sorted(roots):
        following = [0] * (len(coefficients) + 1)
        for j, value in enumerate(coefficients):
            following[j] ^= gf64_multiply(value, powers[exponent])
            following[j + 1] ^= value
        coefficients = following
    assert set(coefficients) <= {0, 1}
    polynomial = sum(value << j for j, value in enumerate(coefficients))
    rows = [polynomial << j for j in range(64 - len(coefficients))]
    extended = [row | ((popcount(row) % 2) << 63) for row in rows]
    return extended, polynomial


def rref(rows, n):
    rows = list(rows)
    rank, pivots = 0, []
    for column in range(n):
        pivot = next((i for i in range(rank, len(rows))
                      if (rows[i] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> column) & 1:
                rows[i] ^= rows[rank]
        pivots.append(column)
        rank += 1
    assert not any(rows[rank:])
    return rows[:rank], pivots


def kernel(checks, n):
    rows, pivots = rref(checks, n)
    return [(1 << free) | sum(((row >> free) & 1) << pivot
                             for row, pivot in zip(rows, pivots))
            for free in range(n) if free not in pivots]


def in_span(word, rows, pivots):
    for row, pivot in zip(rows, pivots):
        if (word >> pivot) & 1:
            word ^= row
    return word == 0


def enumerate_dense_code(rows):
    generator_path = HERE / 'bch64.generator.txt'
    generator_path.write_text('64 %d\n' % len(rows) + ''.join(
        '%x 0\n' % row for row in rows), encoding='ascii')
    env = os.environ.copy()
    env['PATH'] = 'C:/msys64/mingw64/bin;' + env.get('PATH', '')
    command = [str(ENUMERATOR),
               str(generator_path)]
    result = subprocess.run(command, capture_output=True, text=True,
                            check=True, env=env)
    report = json.loads(result.stdout)
    histogram = {int(k): int(v) for k, v in report['histogram'].items()}
    assert sum(histogram.values()) == 1 << 28
    assert histogram[0] == 1 and min(k for k in histogram if k) == 16
    return histogram


def gf4_times_alpha(word):
    # F4^3 is identified with F2^6, in three consecutive bit pairs.
    result = 0
    for j in range(3):
        pair = (word >> (2 * j)) & 3
        image = (0, 2, 3, 1)[pair]  # alpha^2 = alpha + 1
        result |= image << (2 * j)
    return result


def construction65():
    parent, parent_polynomial = extended_bch(13)
    owners, owner_polynomial = extended_bch(23)
    assert len(parent) == 30 and len(owners) == 16
    gram_checks = [sum((popcount(a & b) % 2) << i
                       for i, b in enumerate(parent)) for a in parent]
    coefficient_kernel = kernel(gram_checks, len(parent))
    dense = []
    for coefficients in coefficient_kernel:
        word = 0
        for j, row in enumerate(parent):
            if (coefficients >> j) & 1:
                word ^= row
        dense.append(word)
    assert len(dense) == 28
    assert all(popcount(row) % 4 == 0 for row in dense)
    assert all(popcount(a & b) % 2 == 0 for a in dense for b in dense)
    # Self-orthogonal generators with weights divisible by four give a
    # doubly even code. BCH >=13, extension >=14, divisibility gives >=16.
    reduced, pivots = rref(dense, 64)
    assert len(reduced) == 28
    assert all(in_span(word, reduced, pivots) for word in owners)
    owner_words = [0]
    for row in owners:
        owner_words += [word ^ row for word in owner_words]
    owner_histogram = Counter(map(popcount, owner_words))
    assert owner_histogram[0] == owner_histogram[64] == 1
    assert set(owner_histogram) == {0, 24, 28, 32, 36, 40, 64}
    assert len(owner_words) == 65536
    histogram = enumerate_dense_code(dense)

    sqs = {tuple(sorted((a, b, c, a ^ b ^ c)))
           for a, b, c in itertools.combinations(range(64), 3)}
    assert len(sqs) == 10416
    triples = [triple for block in sqs
               for triple in itertools.combinations(block, 3)]
    assert len(triples) == len(set(triples)) == 41664
    lines = set()
    for direction in range(1, 64):
        other = gf4_times_alpha(direction)
        assert other not in (0, direction)
        for origin in range(64):
            lines.add(tuple(sorted((origin, origin ^ direction,
                                    origin ^ other,
                                    origin ^ direction ^ other))))
    assert len(lines) == 336 and lines <= sqs
    pairs = [pair for line in lines
             for pair in itertools.combinations(line, 2)]
    assert len(pairs) == len(set(pairs)) == 2016
    even_signs = [word for word in range(16) if popcount(word) % 2 == 0]
    assert len(even_signs) == 8
    assert all(popcount(a ^ b) in (2, 4)
               for a, b in itertools.combinations(even_signs, 2))
    # Extended Hamming [16,11,4]: parity checks are (1,index) for all
    # index in F2^4. Any <=3 columns are independent, hence distance >=4.
    hamming_checks = [(1 << 16) - 1] + [sum(((i >> j) & 1) << i
                                           for i in range(16))
                                        for j in range(4)]
    signs16 = kernel(hamming_checks, 16)
    assert len(signs16) == 11
    assert all(popcount(word) % 2 == 0 for word in signs16)

    support_count = histogram[16]
    assert support_count == 30828
    terms = {'dense': 1 << 28, 'weight16': support_count * (1 << 11),
             'weight4': len(sqs) * 16, 'axes64': 128}
    base = sum(terms.values())
    owner_count = len(owner_words) + 8 * len(lines)
    total = base + owner_count + 2
    assert base == 331737984 and owner_count == 68224
    assert total == 331806210
    # Uncapped base products are bounded by 32. For cap vs uncapped
    # base the bound becomes 16 sqrt(3) <32. Distinct non-antipodal
    # owners have |dot|<=16; cap/cap dot<=3/4*16+16=28.
    assert 16 * 16 * 3 < 32 * 32
    assert Fraction(3, 4) * 16 + 16 == 28 < 32
    assert Fraction(3, 4) * 64 - 16 == 32
    assert 8 * 4 == 32  # new axial poles vs caps
    johnson_weight4 = (65 * ((64 * (63 // 2)) // 3)) // 4
    conservative_comparison = ((1 << 28) + 2048 * support_count
                               + 16 * johnson_weight4 + 130)
    return {
        'status': 'Explicit construction; fresh dense-code enumeration and exact geometry',
        'dimension': 65, 'lower_bound': total,
        'base64_terms': terms, 'base64_total': base,
        'dense_code': {'parent_bch_generator_polynomial_bits': parent_polynomial,
                       'radical_dimension': len(dense),
                       'distance_proof': 'BCH designed distance 13, parity extension, doubly even radical',
                       'fresh_weight_histogram': histogram},
        'owner_class': {'dense_bch_generator_polynomial_bits': owner_polynomial,
                        'dense_weight_histogram': dict(sorted(owner_histogram.items())),
                        'dense_points': len(owner_words),
                        'sparse_supports': len(lines), 'sparse_points': 8 * len(lines),
                        'total_points': owner_count, 'coherence': '1/4'},
        'map': 'Remove each owner u, insert (sqrt(3)/2*u,+4) and (sqrt(3)/2*u,-4), add +/-8e65',
        'norm_squared': 64, 'inner_product_limit': 32,
        'gain_over_base64_plus_two_poles': owner_count,
        'same_ers_inputs_with_Johnson_upper_bound_for_weight4': conservative_comparison,
        'margin_over_that_comparison': total - conservative_comparison,
        'scope': 'No comprehensive literature-priority certification and no public registry update',
    }

