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


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    path = HERE / 'mixed_pair_certificate.npz'
    with np.load(path, allow_pickle=False) as data:
        low, high = data['supports_low'], data['supports_high']
    assert low.dtype == high.dtype == np.uint64
    assert low.ndim == high.ndim == 1 and low.shape == high.shape and len(low) > 1
    assert np.all(high < 2**32)
    assert np.all(bit_counts(low) + bit_counts(high) == 24)
    masks = [int(lo) | (int(hi) << 64) for lo, hi in zip(low, high)]
    assert len(set(masks)) == len(masks)
    occupied = [sum(bool(mask & (15 << (4*r))) for r in range(24)) for mask in masks]
    assert max(occupied) <= 12
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
    (HERE/'support_bits_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
