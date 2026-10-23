#!/usr/bin/env python3
"""Build and exhaustively verify a binary linear [62,26,16] code.

Follow Grassl's cyclic [73,36,16] -> puncture -> Construction B2 ->
extend -> shorten chain, with an explicit weight-nine dual word fixing B2.
The final minimum distance is checked over ALL 2^26 words. Thus the
verification does not assume the cited starting code's minimum distance.

Source: https://www.codetables.de/BKLC/BKLC.php?q=2&n=62&k=26
"""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np


EXPONENTS = [37, 36, 34, 33, 32, 27, 25, 24, 22, 21, 19, 18,
             15, 11, 10, 8, 7, 5, 3, 0]
DUAL_SUPPORT = [21, 37, 38, 49, 57, 63, 67, 70, 72]


def rref(rows, n):
    rows = list(rows)
    rank, pivots = 0, []
    for c in range(n):
        k = next((i for i in range(rank, len(rows)) if (rows[i] >> c) & 1), None)
        if k is None:
            continue
        rows[rank], rows[k] = rows[k], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> c) & 1:
                rows[i] ^= rows[rank]
        pivots.append(c)
        rank += 1
        if rank == len(rows):
            break
    assert all(x == 0 for x in rows[rank:])
    return rows[:rank], pivots


def impose_zeros(rows, coordinates):
    rows = list(rows)
    for c in coordinates:
        k = next((i for i, x in enumerate(rows) if (x >> c) & 1), None)
        if k is None:
            continue
        pivot = rows.pop(k)
        rows = [x ^ pivot if (x >> c) & 1 else x for x in rows]
    assert all(not ((x >> c) & 1) for x in rows for c in coordinates)
    return rows


def keep_coordinates(rows, keep):
    return [sum(((x >> c) & 1) << j for j, c in enumerate(keep)) for x in rows]


def build_code64_rows():
    g = sum(1 << e for e in EXPONENTS)
    dividend, quotient = (1 << 73) | 1, 0
    while dividend.bit_length() >= g.bit_length():
        shift = dividend.bit_length() - g.bit_length()
        quotient ^= 1 << shift
        dividend ^= g << shift
    assert dividend == 0
    original = [g << i for i in range(36)]
    assert len(rref(original, 73)[0]) == 36
    dual = sum(1 << c for c in DUAL_SUPPORT)
    assert len(DUAL_SUPPORT) == 9 and 0 not in DUAL_SUPPORT
    assert all(bin(x & dual).count('1') % 2 == 0 for x in original)
    # Puncture coordinate 0, and shorten the nine coordinates in DUAL_SUPPORT.
    # Their projection has rank eight, hence only eight dimensions are lost.
    shortened = impose_zeros(original, DUAL_SUPPORT)
    assert len(shortened) == 28
    remaining = [i for i in range(73) if i != 0 and i not in DUAL_SUPPORT]
    assert len(remaining) == 63
    code63 = keep_coordinates(shortened, remaining)
    assert len(rref(code63, 63)[0]) == 28
    extended = [x | ((bin(x).count('1') % 2) << 63) for x in code63]
    extended, _ = rref(extended, 64)
    return extended


def build_code():
    extended = build_code64_rows()
    code62 = impose_zeros(extended, [62, 63])
    assert len(code62) == 26
    code62 = keep_coordinates(code62, list(range(62)))
    code62, pivots = rref(code62, 62)
    assert len(code62) == 26
    G = np.array([[(x >> c) & 1 for c in range(62)] for x in code62], dtype=np.uint8)
    assert G[:, pivots].tolist() == np.eye(26, dtype=np.uint8).tolist()
    return G


def popcounts64(x):
    """Exact vectorized SWAR popcount, with unsigned wraparound multiplication."""
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return ((x * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.uint8)


def weight_enumerator(G, batch_bits=18):
    rows, pivots = rref([sum(int(v) << i for i, v in enumerate(row)) for row in G], G.shape[1])
    assert len(rows) == len(G) == len(pivots) and G.shape[1] <= 64
    low_bits = min(batch_bits, len(rows))
    low = np.zeros(1 << low_bits, dtype=np.uint64)
    for i in range(low_bits):
        low[1 << i:1 << (i + 1)] = low[:1 << i] ^ np.uint64(rows[i])
    high = [0]
    for row in rows[low_bits:]:
        high += [x ^ row for x in high]
    histogram = np.zeros(G.shape[1] + 1, dtype=np.int64)
    for x in high:
        histogram += np.bincount(popcounts64(low ^ np.uint64(x)), minlength=len(histogram))
    assert histogram.sum() == 1 << len(G) and histogram[0] == 1
    return histogram


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.npz'))
    args = parser.parse_args()
    start = perf_counter()
    G = build_code()
    histogram = weight_enumerator(G)
    distance = int(np.flatnonzero(histogram[1:])[0] + 1)
    assert distance == 16
    np.savez_compressed(args.output, generator=G, weight_histogram=histogram)
    report = {'parameters': [62, 26, distance], 'status': 'built and exhaustively verified',
              'words_enumerated': int(histogram.sum()),
              'weight_histogram': {str(i): int(x) for i, x in enumerate(histogram) if x},
              'construction': {'cyclic_length': 73, 'generator_exponents': EXPONENTS,
                               'punctured_coordinate': 0, 'dual_weight9_support': DUAL_SUPPORT,
                               'shortened_after_extension': [62, 63], 'coordinate_indexing': 'zero based'},
              'generator_sha256': hashlib.sha256(G.tobytes()).hexdigest(),
              'source_chain': 'https://www.codetables.de/BKLC/BKLC.php?q=2&n=62&k=26',
              'seconds': perf_counter() - start}
    args.output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
