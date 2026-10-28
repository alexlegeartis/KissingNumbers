"""Check the complete support family using packed integer intersections.

This replay is independent of the coordinate matrix products in
mixed_pair_certificate.py. It needs no search state or code membership claim.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def bit_counts(words):
    """Parallel population counts, entirely in unsigned 64-bit arithmetic."""
    words = words - ((words >> np.uint64(1)) & np.uint64(0x5555555555555555))
    words = ((words & np.uint64(0x3333333333333333)) +
             ((words >> np.uint64(2)) & np.uint64(0x3333333333333333)))
    words = (words + (words >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    # Add the eight byte counts explicitly; this also avoids modular products.
    total = np.zeros(words.shape, dtype=np.uint8)
    for shift in range(0, 64, 8):
        total += ((words >> np.uint64(shift)) & np.uint64(255)).astype(np.uint8)
    return total


def check_obstructions(data, masks, occupied):
    # Independent polynomial-basis GF(4) multiplication (x^2=x+1).
    def mul(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & 4:
                a ^= 7
        return value
    def trace(a):
        value = a ^ mul(a, a)
        assert value in (0, 1)
        return value
    generator = data['quaternary_generator']
    assert generator.shape == (9, 24) and generator.dtype == np.uint8
    assert np.all(generator <= 3)
    # Enumerate the outer code independently to check size and distance.
    words = np.zeros((4**9, 24), dtype=np.uint8)
    for j, row in enumerate(generator):
        for a in (1, 2, 3):
            multiple = np.array([mul(a, int(v)) for v in row], dtype=np.uint8)
            words[a*4**j:(a+1)*4**j] = words[:4**j] ^ multiple
    weights = np.count_nonzero(words, axis=1)
    assert np.count_nonzero(weights == 0) == 1 and min(weights[1:]) == 12
    witness = data['dual_witnesses']
    assert witness.dtype == np.uint64 and witness.shape == (len(masks), 2)
    assert np.all(witness < 2**48)
    total = bound = 0
    for index, mask in enumerate(masks):
        used = set()
        misses = 0
        for encoded in witness[index]:
            encoded = int(encoded)
            if not encoded:
                continue
            d = [(encoded >> (2*r)) & 3 for r in range(24)]
            active = {r for r in range(24) if d[r]}
            assert not used & active
            used |= active
            for row in generator:
                product = 0
                for a, b in zip(d, row):
                    product ^= mul(a, int(b))
                assert product == 0
            parity = 0
            for r in active:
                bits = {trace(mul(d[r], c)) for c in range(4)
                        if (mask >> (4*r+c)) & 1}
                assert len(bits) == 1
                parity ^= bits.pop()
            assert parity == 1
            misses += 1
        assert occupied[index]-misses <= 12
        bound = max(bound, occupied[index]-misses)
        total += misses
    return dict(occupied_row_histogram=np.bincount(occupied, minlength=15).tolist(),
                dual_obstruction_supports=sum(v > 12 for v in occupied),
                disjoint_dual_checks=total, maximum_transversal_intersection_bound=bound)


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    path = HERE / 'mixed_pair_certificate.npz'
    with np.load(path, allow_pickle=False) as data:
        data = dict(data)
        low, high = data['supports_low'], data['supports_high']
    assert low.dtype == high.dtype == np.uint64
    assert low.ndim == high.ndim == 1 and low.shape == high.shape and len(low) > 1
    assert np.all(high < 2**32)
    assert np.all(bit_counts(low) + bit_counts(high) == 24)
    masks = [int(lo) | (int(hi) << 64) for lo, hi in zip(low, high)]
    assert len(set(masks)) == len(masks)
    occupied = [sum(bool(mask & (15 << (4*r))) for r in range(24)) for mask in masks]
    assert max(occupied) <= 14
    obstruction_report = check_obstructions(data, masks, occupied)
    histogram = [0]*9
    for mask in masks:
        histogram[bin(mask & 255).count('1')] += 1
    maximum = pairs = 0
    block = 256
    for i in range(0, len(low), block):
        end_i = min(i + block, len(low))
        for j in range(i, len(low), block):
            end_j = min(j + block, len(low))
            intersections = (bit_counts(low[i:end_i, None] & low[None, j:end_j]) +
                             bit_counts(high[i:end_i, None] & high[None, j:end_j]))
            if i == j:
                intersections = intersections[np.triu_indices(end_i-i, 1)]
            if intersections.size:
                maximum = max(maximum, int(intersections.max()))
                assert np.all(intersections <= 12), (i, j, maximum)
                pairs += intersections.size
    assert pairs == len(low)*(len(low)-1)//2
    report = dict(status='PASS: complete packed integer support replay',
                  supports=len(low), support_pairs_checked=pairs,
                  maximum_distinct_support_intersection=maximum,
                  maximum_occupied_rows=max(occupied),
                  section_intersection_histogram=histogram,
                  certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    report.update(obstruction_report)
    (HERE/'support_bits_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
