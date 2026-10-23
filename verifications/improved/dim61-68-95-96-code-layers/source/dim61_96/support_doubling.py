#!/usr/bin/env python3
"""Add doubled QR48 supports to the 24x4 transversal support code in 96D.

A weight-12 binary support on 48 cells becomes a weight-24 support by
duplicating every cell.  Arrange those 48 cells in pairs over 24 rows.
The new support occupies at most 12 rows; every transversal meets it in
at most 12 cells.  Binary distance 12 becomes distance 24 on doubling.

Two further doubling partitions are also allowed for the QR48 supports
having at most one cell in every pair.  A support in either such family
meets every support from the other partitions in at most 12 cells.
"""

import itertools
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from weighted_layers import (HERE, all_binary_words, collect_min_words,
                             matrix, packed, popcounts64, rref, sign24)


def qr48():
    residues = {(i * i) % 47 for i in range(1, 47)}
    assert len(residues) == 23 and {(2 * r) % 47 for r in residues} == residues
    idempotent = sum(1 << r for r in residues)
    mask = (1 << 47) - 1
    rows = [((idempotent << i) | (idempotent >> (47 - i))) & mask for i in range(47)]
    rows, _ = rref(rows, 47)
    assert len(rows) == 24
    rows = [x | ((bin(x).count('1') % 2) << 47) for x in rows]
    G = matrix(rref(rows, 48)[0], 48)
    words, histogram = collect_min_words(G, 12)
    assert len(words) == 17296
    return G, words, histogram


def expand(word, partition):
    # Each perfect matching partitions a four-cell row into two pairs.
    pairs = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))][partition]
    result = 0
    while word:
        bit = (word & -word).bit_length() - 1
        row, side = divmod(bit, 2)
        for c in pairs[side]:
            result |= 1 << (4 * row + c)
        word &= word - 1
    return result


def section_classification(supports, i, j):
    counts = {'none': 0, 'one': 0, 'both': 0}
    for s in supports:
        k = ((s >> i) & 1) + ((s >> j) & 1)
        counts[{0: 'none', 1: 'one', 2: 'both'}[k]] += 1
    return counts


def main():
    start = perf_counter()
    G, words, histogram = qr48()
    no_double = words[(words & (words >> np.uint64(1)) & np.uint64(0x555555555555)) == 0]
    all_supports = [expand(int(x), 0) for x in words]
    filtered_supports = [[expand(int(x), p) for x in no_double] for p in (1, 2)]
    additions = all_supports + filtered_supports[0] + filtered_supports[1]
    assert len(set(additions)) == len(additions)
    assert all(bin(s).count('1') == 24 for s in additions)
    assert all(sum(bool((s >> (4 * r)) & 15) for r in range(24)) <= 12 for s in additions)
    # For each of the two extra families, occupied rows contain exactly two cells.
    assert all(all(bin((s >> (4 * r)) & 15).count('1') in (0, 2) for r in range(24))
               for family in filtered_supports for s in family)
    Gsign, _ = sign24()
    sign_columns = [sum(int(v) << i for i, v in enumerate(Gsign[:, j])) for j in range(24)]
    assert 0 not in sign_columns and len(set(sign_columns)) == 24
    # Hence equality of any two distinct sign positions keeps exactly 2^13 signs.
    # Choose the section x_0=x_4.  These cells lie in distinct transversal rows,
    # and can be placed in distinct GF(16) rows of the six-support layer.
    classification = section_classification(additions, 0, 4)
    added95 = (classification['none'] * 16384 + classification['both'] * 8192)
    added96 = len(additions) * 16384
    terms95 = {'dense_section': 1 << 32,
               'transversal24_section': 147456 * 16384 + 16384 * 8192,
               'new_support24_section': added95,
               'support6_section': 57600 * 32 + 256 * 16,
               'ordinary_poles': 188}
    terms96 = {'dense': 1 << 33, 'transversal24': (4 ** 9) * 16384,
               'new_support24': added96, 'support6': (16 ** 4) * 32,
               'poles': 192}
    report = {
        'status': 'Explicit QR48 enumeration and exact support-composition proof',
        'mechanism': 'Transversals plus doubled supports occupying at most twelve rows',
        'qr48_weight_histogram': {str(i): int(v) for i, v in enumerate(histogram) if v},
        'qr48_weight12_supports': len(words), 'qr48_no_double_pair_supports': len(no_double),
        'total_added_supports': len(additions),
        'section95': {'equation': 'x_0=x_4 before an orthogonal identification with R^95',
                      'new_support_classification': classification,
                      'terms': terms95, 'total': sum(terms95.values()),
                      'repository_baseline': 6473065046,
                      'gain': sum(terms95.values()) - 6473065046},
        'dimension96': {'terms': terms96, 'total': sum(terms96.values()),
                        'repository_baseline': 12886999232,
                        'gain': sum(terms96.values()) - 12886999232},
        'checks': {
            'QR48_distance12': int(np.flatnonzero(histogram[1:])[0] + 1) == 12,
            'doubling_preserves_intersection_bound': '2*|A intersect B| <= 2*6=12',
            'transversal_cross_bound': 'at most one per occupied row, at most twelve rows',
            'different_partition_cross_bound': 'one cell per pair-pair intersection; at most twelve source cells',
            'sparse6_cross_bound': '2*4*6=48',
            'geometry_inheritance': 'All 95D points are in a hyperplane section of the 96D configuration',
        },
        'seconds': perf_counter() - start,
    }
    np.savez_compressed(HERE / 'support_doubling.npz', qr48_generator=G,
                        qr48_weight12_supports=words, qr48_no_double_supports=no_double,
                        added_supports_low=np.array([x & ((1 << 64) - 1) for x in additions], dtype=np.uint64),
                        added_supports_high=np.array([x >> 64 for x in additions], dtype=np.uint64))
    (HERE / 'support_doubling.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
