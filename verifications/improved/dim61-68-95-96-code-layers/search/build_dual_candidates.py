"""Rebuild the candidate pools for the pinned dual-obstruction support walk.

Usage: python build_dual_candidates.py /tmp/dual-walk
Only the final certificate verifiers establish geometric correctness.
"""
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
bits = getattr(int, 'bit_count', lambda word: bin(word).count('1'))


def multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


def row_mask(word):
    return sum(1 << r for r in range(24) if (word >> (2*r)) & 3)


def error_syndromes(columns, weight):
    buckets = {}
    for support in combinations(range(24), weight):
        mask = sum(1 << r for r in support)
        for values in product(range(3), repeat=weight):
            syndrome = word = 0
            for r, a in zip(support, values):
                syndrome ^= columns[r][a]
                word |= (a+1) << (2*r)
            buckets.setdefault(syndrome, []).append((word, mask))
    return buckets


def join_errors(left, right, same=False):
    result = set()
    for syndrome, group in left.items():
        pairs = combinations(group, 2) if same else product(group, right.get(syndrome, []))
        for (a, ma), (b, mb) in pairs:
            if not ma & mb:
                result.add(a ^ b)
    return sorted(result)


def restricted_duals(columns):
    def half(coordinates):
        syndrome = np.zeros(1, dtype=np.uint32)
        words = np.zeros(1, dtype=np.uint64)
        weights = np.zeros(1, dtype=np.uint8)
        for r in coordinates:
            masks = [columns[r][a-1] for a in (2, 3)]
            syndrome = np.r_[syndrome, syndrome ^ masks[0], syndrome ^ masks[1]]
            words = np.r_[words, words | np.uint64(2 << (2*r)), words | np.uint64(3 << (2*r))]
            weights = np.r_[weights, weights+1, weights+1]
        return syndrome, words, weights
    a, wa, ha = half(range(12))
    b, wb, hb = half(range(12, 24))
    order = np.argsort(b)
    b, wb, hb = b[order], wb[order], hb[order]
    low, high = np.searchsorted(b, a, 'left'), np.searchsorted(b, a, 'right')
    counts = high-low
    left = np.repeat(np.arange(len(a)), counts)
    right = np.repeat(low, counts) + np.arange(int(counts.sum())) - np.repeat(np.cumsum(counts)-counts, counts)
    weight = ha[left]+hb[right]
    keep = (weight >= 6) & (weight <= 10)
    result = set(map(int, wa[left[keep]] | wb[right[keep]]))
    assert len(result) == 11152
    return result


def main():
    target = Path(sys.argv[1])
    target.mkdir(parents=True, exist_ok=True)
    with np.load(HERE/'dual_obstruction_start.npz', allow_pickle=False) as data:
        generator, start = data['quaternary_generator'], data['nonordinary_supports']
    columns = [[sum(multiply(a, int(v)) << (2*j) for j, v in enumerate(generator[:, r]))
                for a in (1, 2, 3)] for r in range(24)]
    triples, quads = error_syndromes(columns, 3), error_syndromes(columns, 4)
    d6 = join_errors(triples, triples, same=True)
    d7 = join_errors(triples, quads)
    d8 = join_errors(quads, quads, same=True)
    assert [len(d6), len(d7), len(d8)] == [180, 2880, 20925]
    full = d6+d7+d8
    pool13 = restricted_duals(columns)
    for d in full:
        if sum(((d >> (2*r)) & 3) == 1 for r in range(24)) <= 2:
            pool13.add(d)
    pool13 = sorted(pool13)
    assert len(pool13) == 21296
    annotated = [(d, row_mask(d), sum(1 << r for r in range(24)
                                     if ((d >> (2*r)) & 3) == 1)) for d in full]
    pairs = []
    for index, (d, md, ad) in enumerate(annotated):
        if bits(md) > 7:
            continue
        for e, me, ae in annotated[index+1:]:
            if bits(md)+bits(me) <= 14 and not md & me and bits(ad | ae) <= 4:
                pairs.append((d, e))
    assert len(pairs) == 247800
    (target/'duals13.txt').write_text(str(len(pool13))+'\n'+''.join(f'{d}\n' for d in pool13))
    (target/'pairs14.txt').write_text(str(len(pairs))+'\n'+''.join(f'{d} {e}\n' for d, e in pairs))
    (target/'start.txt').write_text(str(len(start))+'\n'+''.join(' '.join(map(str, row))+'\n' for row in start))
    print({'start_supports':len(start), 'dual_pool13':len(pool13), 'dual_pairs14':len(pairs)})


if __name__ == '__main__':
    main()
