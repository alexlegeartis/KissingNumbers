#!/usr/bin/env python3
"""Central sections of the explicit [64,28,16] kissing construction.

For a weight-eight dual word T, the projected binary code on T is the
even-weight length-eight code.  The hyperplane sum_{i in T} x_i=0 keeps
70 of its 128 patterns, rather than half the dense population.

The signs on each weight-16 support may be permuted separately.  Put its
intersection with T first in the extended Hamming sign code, making an
eight-point intersection an affine half of that code's sixteen positions.
Every section inherits the exact sphere and distance constraints.
"""

import json
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np

from weighted_layers import (HERE, all_binary_words, build_code64_rows,
                             collect_min_words, kernel, matrix, packed,
                             popcounts64, rref, sqs64)


def positions(x, n=64):
    return [i for i in range(n) if (x >> i) & 1]


def project_rows(G, coordinates):
    return [sum(int(row[c]) << j for j, c in enumerate(coordinates)) for row in G]


def dense_section_count(G, sections):
    union = sorted(set(c for s in sections for c in s))
    rowspace = rref(project_rows(G, union), len(union))[0]
    projected = all_binary_words(matrix(rowspace, len(union)))
    surviving = np.ones(len(projected), dtype=bool)
    for s in sections:
        mask = sum(1 << union.index(c) for c in s)
        surviving &= popcounts64(projected & np.uint64(mask)) == len(s) // 2
    return int(surviving.sum()) * (1 << (len(G) - len(rowspace))), len(rowspace)


def main():
    start = perf_counter()
    cache = HERE / 'code64_supports.npz'
    if cache.exists():
        data = np.load(cache)
        G, words, dual = data['generator'], data['words'], data['dual_generator']
    else:
        G = matrix(build_code64_rows(), 64)
        words, _ = collect_min_words(G, 16)
        dual = matrix(kernel(packed(G), 64), 64)
    assert np.array_equal(G, matrix(build_code64_rows(), 64))
    assert len(words) == 30828 and len(np.unique(words)) == len(words)
    assert np.all(popcounts64(words) == 16)
    remainder = words.copy()
    code_rows, code_pivots = rref(packed(G), 64)
    for row, pivot in zip(code_rows, code_pivots):
        selected = ((remainder >> np.uint64(pivot)) & np.uint64(1)).astype(bool)
        remainder[selected] ^= np.uint64(row)
    assert not np.any(remainder)
    # The dual basis and its pair sums suffice to exhibit disjoint weight-eight words.
    dual_rows = packed(dual)
    small = set(x for x in dual_rows if bin(x).count('1') == 8)
    small.update(x ^ y for i, x in enumerate(dual_rows) for y in dual_rows[:i]
                 if bin(x ^ y).count('1') == 8)
    T1, T2 = next((a, b) for a in sorted(small) for b in sorted(small) if a < b and not (a & b))
    assert all(bin(row & T).count('1') % 2 == 0 for row in packed(G) for T in (T1, T2))
    hamming_checks = [(1 << 16) - 1]
    hamming_checks += [sum(((j >> i) & 1) << j for j in range(16)) for i in range(4)]
    sign_G = matrix(kernel(hamming_checks, 16), 16)
    sign_words = all_binary_words(sign_G)
    assert np.min(popcounts64(sign_words)[1:]) == 4
    base_sqs = sqs64()
    reports = {}
    for sections_masks in [(T1,), (T1, T2)]:
        sections = [positions(T) for T in sections_masks]
        dim = 64 - len(sections)
        dense, projection_rank = dense_section_count(G, sections)
        intersection = [popcounts64(words & np.uint64(T)) for T in sections_masks]
        support_sign_counts, shape_hist = {}, {}
        sparse_count = 0
        for shape in zip(*intersection):
            shape = tuple(sorted((int(j) for j in shape), reverse=True))
            shape_hist[shape] = shape_hist.get(shape, 0) + 1
        for shape, frequency in shape_hist.items():
            valid, offset = np.ones(len(sign_words), dtype=bool), 0
            for j in shape:
                if j % 2:
                    valid[:] = False
                else:
                    mask = ((1 << j) - 1) << offset
                    valid &= popcounts64(sign_words & np.uint64(mask)) == j // 2
                offset += j
            retained = int(valid.sum())
            support_sign_counts[shape] = retained
            sparse_count += frequency * retained
        # Relabel SQS64 so disjoint T's are consecutive affine 3-spaces
        # {0,...,7}, {8,...,15}. All other coordinates occupy the remaining positions.
        section4 = 0
        block_shape = {}
        for block in base_sqs:
            shape = tuple(sum(8 * i <= c < 8 * (i + 1) for c in block)
                          for i in range(len(sections)))
            block_shape[shape] = block_shape.get(shape, 0) + 1
            if all(j % 2 == 0 for j in shape):
                signs = 2 ** (4 - sum(shape))
                for j in shape:
                    signs *= comb(j, j // 2)
                section4 += signs
        poles = 2 * (64 - 8 * len(sections))
        terms = {'dense_section': dense, 'weight16_section': sparse_count,
                 'weight4_section': section4, 'ordinary_poles': poles}
        baseline = {62: 71310732, 63: 138419844}[dim]
        reports[str(dim)] = {
            'hyperplane_supports': sections, 'projection_rank': projection_rank,
            'terms': terms, 'total': sum(terms.values()), 'repository_baseline': baseline,
            'gain': sum(terms.values()) - baseline,
            'weight16_intersection_shapes': {str(s): v for s, v in sorted(shape_hist.items())},
            'retained_signs_per_shape': {str(s): v for s, v in sorted(support_sign_counts.items())},
            'weight4_intersection_shapes': {str(s): v for s, v in sorted(block_shape.items())},
            'sign_order': 'Put larger intersection first, then the smaller intersection, then the rest',
        }
    assert reports['63']['terms']['dense_section'] == (1 << 21) * 70
    assert reports['63']['total'] == 182407428
    report = {'status': 'Exact binary projection counts and exact hyperplane sections',
              'dimensions': reports, 'seconds': perf_counter() - start,
              'geometry_proof': 'Take central sections of a norm-64 ERS configuration; orthogonal identification preserves all norms and inner products.'}
    (HERE / 'dual_sections.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    np.savez_compressed(HERE / 'dual_sections.npz', generator64=G,
                        weight16_supports=words, sign_generator16=sign_G,
                        section_supports=np.array([positions(T1), positions(T2)], dtype=np.uint8))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
