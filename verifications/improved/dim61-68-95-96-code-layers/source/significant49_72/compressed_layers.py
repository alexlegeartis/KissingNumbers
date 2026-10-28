"""Exact 61D parity layers and 62D hexagonal compression.

Rebuild and enumerate the [64,28,16] code.  No search output is a proof
input.  Coordinates are zero based; the retained head has 60 coordinates.
All geometric comparisons use rational numbers and quadratic radicals.
"""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
from math import comb
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dim61_96'))
from weighted_layers import (build_code64_rows, collect_min_words, impose_zeros,
                             kernel, matrix, packed, popcounts64, rref, sqs64)
sys.path.insert(0, str(HERE.parent / 'dim49plus'))
from latitude_compression import geometry, quad_leq

FOUR = (0, 1, 12, 52)
HEX_LABELS = (0, 1, 3, 7, 6, 4)
# Coordinates (x, y*sqrt(3)); multiply both by 1/8.
HEX = ((2, 0), (1, 1), (-1, 1), (-2, 0), (-1, -1), (1, -1))
# Python 3.9 is included in the repository's CI matrix.
bits = getattr(int, 'bit_count', lambda word: bin(word).count('1'))


def words(rows):
    result = [0]
    for row in rows:
        result += [word ^ row for word in result]
    return result


def reduce_word(word, rows, pivots):
    for row, pivot in zip(rows, pivots):
        if word >> pivot & 1:
            word ^= row
    return word


def best_four(shell):
    """Reproduce the entire 635,376-choice shortening search with bitsets."""
    masks = []
    for coordinate in range(64):
        bitmap = bytearray((len(shell)+7)//8)
        for i, support in enumerate(shell):
            if not (support >> coordinate & 1):
                bitmap[i//8] |= 1 << (i % 8)
        masks.append(int.from_bytes(bitmap, 'little'))
    histogram, best, choice = Counter(), -1, None
    for a, b, c, d in itertools.combinations(range(64), 4):
        count = bits(masks[a] & masks[b] & masks[c] & masks[d])
        histogram[count] += 1
        if count > best:
            best, choice = count, (a, b, c, d)
    assert sum(histogram.values()) == comb(64, 4)
    assert choice == FOUR and best == 9673
    return dict(choices=comb(64, 4), maximum=best,
                maximizing_choices=histogram[best],
                count_histogram={str(k): histogram[k] for k in sorted(histogram)})


def hex_geometry(labels=HEX_LABELS):
    gram, violations = [], []
    for a, (x, y) in enumerate(HEX):
        assert x*x+3*y*y == 4
        line = []
        for b, (u, v) in enumerate(HEX):
            tail_dot = x*u+3*y*v
            line.append(tail_dot)
            removed_distance = bits(labels[a] ^ labels[b])
            # For distinct codewords head distance is >= 16-removed_distance.
            worst_numerator = 60-2*(16-removed_distance)+tail_dot
            if worst_numerator > 32:
                violations.append([a, b, removed_distance, worst_numerator])
        gram.append(line)
    return gram, violations


def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are required.')
    rows, pivots = rref(build_code64_rows(), 64)
    assert len(rows) == 28
    # This orthogonality is the extra fact missed by the absolute-value bound.
    assert all(bits(a & b) % 2 == 0 for a in rows for b in rows)
    G = matrix(rows, 64)
    raw_shell, histogram = collect_min_words(G, 16)
    shell = sorted(map(int, raw_shell))
    assert len(shell) == len(set(shell)) == 30828
    assert histogram.sum() == 1 << 28 and not any(histogram[1:16])
    assert all(reduce_word(s, rows, pivots) == 0 for s in shell)
    with np.load(HERE.parent/'dim61_96/code64_supports.npz', allow_pickle=False) as archive:
        assert shell == sorted(map(int, archive['words']))
        assert rref(packed(archive['generator']), 64)[0] == rows
    search = best_four(shell)
    removed = sum(1 << c for c in FOUR)
    supports = [s for s in shell if not (s & removed)]
    assert len(supports) == 9673
    head = [c for c in range(64) if c not in FOUR]
    assert len(head) == 60
    assert len(rref([sum((row >> c & 1) << j for j, c in enumerate(FOUR))
                     for row in rows], 4)[0]) == 4
    shortened = impose_zeros(rows, FOUR[2:])
    assert len(shortened) == 26
    assert len(impose_zeros(rows, FOUR)) == 24
    # Also check support intersections directly, independently of linear distance.
    pair_count, maximum_intersection = 0, 0
    array = np.array(supports, dtype=np.uint64)
    for i, support in enumerate(supports[:-1]):
        intersections = popcounts64(array[i+1:] & np.uint64(support))
        maximum_intersection = max(maximum_intersection, int(intersections.max()))
        pair_count += len(intersections)
    assert pair_count == len(supports)*(len(supports)-1)//2
    assert maximum_intersection == 8
    # Extended Hamming [16,11,4].  Its single-bit translate has ODD parity.
    checks = [65535] + [sum(((i >> j) & 1) << i for i in range(16)) for j in range(4)]
    even_signs = words(kernel(checks, 16))
    odd_signs = [word ^ 1 for word in even_signs]
    assert len(even_signs) == len(set(odd_signs)) == 2048
    assert min(bits(word) for word in even_signs if word) == 4
    assert all(bits(word) % 2 == 0 for word in even_signs)
    assert all(bits(word) % 2 == 1 for word in odd_signs)
    # Expanded coordinates independently check the translated sign distances.
    signs = np.array([[1-2*(word >> c & 1) for c in range(16)]
                      for word in odd_signs], dtype=np.int16)
    dots = signs @ signs.T
    assert np.all(np.diag(dots) == 16)
    np.fill_diagonal(dots, -16)
    assert int(dots.max()) == 8
    # Dense 61D words have even restriction to every selected support.
    # A sparse odd-parity word must disagree at least once on that support.
    assert all(bits(row & s) % 2 == 0 for row in shortened for s in supports)
    assert geometry()['valid']
    assert Q(14, 32) <= Q(1, 2)
    assert quad_leq(0, Q(14, 128), 17)
    # With the parity reversed, matching all signs is a genuine counterexample.
    witness = next(row for row in shortened if row >> FOUR[1] & 1)
    support = supports[0]
    restriction = sum((witness >> c & 1) << j
                      for j, c in enumerate(c for c in head if support >> c & 1))
    assert bits(restriction) % 2 == 0
    bad_even_coset = [word ^ restriction for word in even_signs]
    assert restriction in bad_even_coset
    assert not quad_leq(0, Q(16, 128), 17)
    gram, violations = hex_geometry()
    assert not violations
    bad_labels = list(HEX_LABELS)
    bad_labels[1], bad_labels[2] = bad_labels[2], bad_labels[1]
    assert hex_geometry(bad_labels)[1]
    assert 60+4 == 64  # Every hexagonal dense point has norm one.
    assert Q(16, 8*4) == Q(1, 2)  # 62D dense against weight sixteen.
    # All four-support signs, on affine planes restricted to the new ambient space.
    sqs = sqs64()
    quadruples = {d: [block for block in sqs if max(block) < d] for d in (61, 62)}
    for dimension, blocks in quadruples.items():
        triples = [t for block in blocks for t in itertools.combinations(block, 3)]
        assert len(triples) == len(set(triples)) == 4*len(blocks)
    assert [len(quadruples[d]) for d in (61, 62)] == [8555, 9145]
    # 61D dense against four-supports: either zero or one tail coordinate.
    for j in (0, 1):
        for tail_count in (0, 1):
            heads = 4-tail_count
            if j == 0:
                assert Q(heads, 16)+Q(tail_count, 8) <= Q(1, 2)
            else:
                assert quad_leq(Q(tail_count, 32), Q(heads, 64), 17)
    # 62D dense against four-supports, including each subset of the two tails.
    for x, y in HEX:
        for subset in ((), (0,), (1,), (0, 1)):
            rational = Q(4-len(subset), 16)
            radical = Q(0)
            if 0 in subset:
                rational += Q(abs(x), 16)
            if 1 in subset:
                radical += Q(abs(y), 16)
            assert quad_leq(rational, radical, 3)
    # Within and between all sparse layers, and against the coordinate poles.
    assert Q(16-2*4, 16) == Q(8, 16) == Q(1, 2)
    assert Q(4-2, 4) == Q(2, 4) == Q(1, 2)
    assert Q(4, 4*2) == Q(1, 2)
    assert Q(1, 4) <= Q(1, 2) and Q(1, 2) <= Q(1, 2)
    assert quad_leq(0, Q(1, 32), 17)  # largest 61D head coordinate
    assert Q(1, 4) <= Q(1, 2)  # all tail coordinates, both constructions
    dimensions = {}
    previous = {61: 72376856, 62: 100437548}
    for dimension in (61, 62):
        terms = dict(dense=(1 << 26) if dimension == 61 else 6*(1 << 24),
                     weight16=len(supports)*len(even_signs),
                     weight4=len(quadruples[dimension])*16, poles=2*dimension)
        total = sum(terms.values())
        dimensions[str(dimension)] = dict(total=total, previous_bound=previous[dimension],
            gain=total-previous[dimension], terms=terms)
    report = dict(status='Verified: exhaustive parent enumeration and exact compressed-layer geometry',
        dimensions=dimensions, removed_coordinates=list(FOUR), head_coordinates=head,
        parent_words_enumerated=1 << 28,
        parent_weight_histogram={str(i): int(v) for i, v in enumerate(histogram) if v},
        parent_generator_sha256=hashlib.sha256(G.tobytes()).hexdigest(),
        supports=len(supports),
        supports_sha256=hashlib.sha256(np.array(supports, dtype='<u8').tobytes()).hexdigest(),
        support_pairs_checked=pair_count, maximum_support_intersection=maximum_intersection,
        sign_words_per_support=2048, sign_pairs_checked=2048*2047//2,
        shortened_dense_dimension=26, four_coordinate_kernel_dimension=24,
        hex_labels=list(HEX_LABELS), hex_scaled_gram=gram,
        shortening_search=search,
        negative_controls=['even-parity sparse coset with matching dense signs violates the angle',
                           'misordered hexagon labels violate the angle'],
        geometry='Unit vectors; every distinct pair has inner product at most 1/2. '
                 'Dense, weight16, weight4 and axial layers are disjoint by their numbers of nonzero coordinates.')
    (HERE/'compressed_layers.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
