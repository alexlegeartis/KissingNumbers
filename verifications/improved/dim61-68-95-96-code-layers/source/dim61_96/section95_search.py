#!/usr/bin/env python3
"""Select a hyperplane section of the enlarged explicit 96D configuration.

Compare coordinate equalities with sum-zero sections on two complete
four-cell rows.  The latter use a weight-eight dual relation of the
dense binary code.  Sign permutations may be selected separately for
each sparse support, since different supports already meet in <=12 cells.
"""

import itertools
import json
from math import comb
from time import perf_counter

import numpy as np

from weighted_layers import (HERE, all_binary_words, dense96, field_mul,
                             kernel, matrix, packed, popcounts64,
                             quaternary24, rref, sign24)


def partitions(total, slots, maximum=None):
    if not slots:
        if total == 0:
            yield ()
        return
    maximum = min(total, maximum if maximum is not None else total)
    for first in range(maximum, -1, -1):
        for tail in partitions(total - first, slots - 1, first):
            yield (first,) + tail


def best_rs_section():
    M16 = np.array([[field_mul(a, b, 4, 0b10011) for b in range(16)] for a in range(16)], dtype=np.uint8)
    G = np.ones((4, 6), dtype=np.uint8)
    for i in range(1, 4):
        G[i] = M16[G[i - 1], np.arange(1, 7)]
    words = np.zeros((16 ** 4, 6), dtype=np.uint8)
    for i, row in enumerate(G):
        size = 16 ** i
        for a in range(1, 16):
            words[a * size:(a + 1) * size] = words[:size] ^ M16[a, row]
    assert np.min(np.count_nonzero(words[1:], axis=1)) == 3
    best = (-1, None, None)
    # j=6 uses the odd sign-code coset, for which all twenty balanced patterns survive.
    sign_count = np.array([32, 0, 16, 0, 12, 0, 20])
    for distribution in partitions(8, 6):
        k = np.count_nonzero(words < np.array(distribution, dtype=np.uint8)[None, :], axis=1)
        histogram = np.bincount(k, minlength=7)
        count = int(histogram @ sign_count)
        if count > best[0]:
            best = count, distribution, histogram.tolist()
    return {'points': best[0], 'selected_cell_counts_per_GF16_row': list(best[1]),
            'intersection_histogram': best[2], 'retained_signs_per_intersection': sign_count.tolist(),
            'selected_symbols': '0,...,count-1 in each row; permute these eight cells to the section coordinates'}


def main():
    start = perf_counter()
    q4, q4_report = quaternary24()
    sign_G, _ = sign24()
    G96 = dense96(sign_G, q4)
    # Every pair of quaternary columns is independent: all sixteen symbol pairs occur.
    normalized_columns = []
    M4 = np.array([[field_mul(a, b, 2, 0b111) for b in range(4)] for a in range(4)], dtype=np.uint8)
    for column in q4.T:
        first = int(column[np.flatnonzero(column)[0]])
        inverse = next(a for a in range(1, 4) if M4[first, a] == 1)
        normalized_columns.append(tuple(M4[inverse, column].tolist()))
    assert len(set(normalized_columns)) == 24
    data = np.load(HERE / 'support_doubling.npz')
    supports = [int(a) | (int(b) << 64) for a, b in zip(data['added_supports_low'], data['added_supports_high'])]
    incidence = np.array([[(s >> c) & 1 for c in range(96)] for s in supports], dtype=np.int32)
    degrees = incidence.sum(axis=0)
    pairs = incidence.T @ incidence
    best_pair = (-1, None)
    for a, b in itertools.combinations(range(96), 2):
        if a // 4 == b // 4:
            continue
        count = int((2 * len(supports) - 2 * degrees[a] - 2 * degrees[b] + 3 * pairs[a, b]) * 8192)
        if count > best_pair[0]:
            best_pair = count, (a, b)
    # Check the dense equality constraint is nonzero on the constructed binary code.
    a, b = best_pair[1]
    assert np.any(G96[:, a] != G96[:, b])
    terms_pair = {'dense_section': 1 << 32,
                  'transversal24_section': 147456 * 16384 + 16384 * 8192,
                  'new_support24_section': best_pair[0],
                  'support6_section': 1847296, 'ordinary_poles': 188}
    # The dual of the 24-bit sign code has minimum eight. An eight-bit dual support
    # can be moved onto any selected eight positions by a permutation of that code.
    dual_G = matrix(kernel(packed(sign_G), 24), 24)
    dual_words = all_binary_words(dual_G)
    dual_weights = popcounts64(dual_words)
    assert np.min(dual_weights[1:]) == 8
    sign_dual8 = int(dual_words[np.flatnonzero(dual_weights == 8)[0]])
    sign_words = all_binary_words(sign_G)
    sign_counts = []
    for j in range(9):
        mask = sign_dual8 if j == 8 else (1 << j) - 1
        sign_counts.append(0 if j % 2 else int(np.count_nonzero(popcounts64(sign_words & np.uint64(mask)) == j // 2)))
    best8 = (-1, None, None)
    row_counts = incidence.reshape(len(supports), 24, 4).sum(axis=2)
    for ra, rb in itertools.combinations(range(24), 2):
        histogram = np.bincount(row_counts[:, ra] + row_counts[:, rb], minlength=9)
        count = int(histogram @ np.array(sign_counts))
        if count > best8[0]:
            best8 = count, (ra, rb), histogram.tolist()
    coordinates = list(range(4 * best8[1][0], 4 * best8[1][0] + 4))
    coordinates += list(range(4 * best8[1][1], 4 * best8[1][1] + 4))
    projection = rref([sum(int(row[c]) << j for j, c in enumerate(coordinates)) for row in G96], 8)[0]
    projected_words = all_binary_words(matrix(projection, 8))
    assert len(projection) == 7 and np.all(popcounts64(projected_words) % 2 == 0)
    dense8 = int(np.count_nonzero(popcounts64(projected_words) == 4)) * (1 << (33 - 7))
    assert dense8 == 4697620480
    rs = best_rs_section()
    terms8 = {'dense_section': dense8, 'transversal24_section': (4 ** 9) * (1 << 13),
              'new_support24_section': best8[0], 'support6_section': rs['points'],
              'ordinary_poles': 176}
    pair_report = {'coordinates': list(best_pair[1]), 'equation': 'x_a=x_b',
                   'terms': terms_pair, 'total': sum(terms_pair.values())}
    eight_report = {'coordinates': coordinates, 'equation': 'sum of these eight coordinates is zero',
                    'terms': terms8, 'total': sum(terms8.values()),
                    'new_support_intersection_histogram': best8[2],
                    'retained_signs_per_intersection': sign_counts,
                    'sign24_dual8_support': [i for i in range(24) if (sign_dual8 >> i) & 1],
                    'rs6_section': rs}
    winner = 'dual8_section' if eight_report['total'] > pair_report['total'] else 'coordinate_equality'
    total = max(eight_report['total'], pair_report['total'])
    report = {'status': 'Exact projection and sign counts in a 96D configuration',
              'coordinate_equality': pair_report, 'dual8_section': eight_report,
              'best_method': winner, 'dimension95_total': total,
              'repository_baseline': 6473065046, 'gain': total - 6473065046,
              'seconds': perf_counter() - start,
              'geometry_proof': 'Central sections and orthogonal identifications preserve the verified 96D distances'}
    np.savez_compressed(HERE / 'section95_search.npz', dense_generator96=G96,
                        support_generator_GF4=q4, sign_generator24=sign_G,
                        selected_coordinates=np.array(coordinates, dtype=np.uint8))
    (HERE / 'section95_search.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
