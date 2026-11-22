"""Bound a two-copy, four-head-coordinate lift of the repository's 39D QR construction."""
import json
from pathlib import Path

import numpy as np


SOURCE = Path(__file__).resolve().parents[2] / "verifications/improved/dim39-qr-three-coordinate-heads"


def quotient(certificate, parent):
    pivots = {}
    for generator in parent["generators_hex"]:
        v = int(generator, 16)
        while v:
            i = v.bit_length() - 1
            if i not in pivots:
                pivots[i] = v
                break
            v ^= pivots[i]
    assert len(pivots) == 15
    free = sorted(set(range(39)) - pivots.keys())

    def syndrome(v):
        for i in sorted(pivots, reverse=True):
            if v & (1 << i):
                v ^= pivots[i]
        return sum(((v >> i) & 1) << j for j, i in enumerate(free))

    marks = certificate["marked_axes"]
    columns = [syndrome(1 << j) for j in range(39) if j not in marks]
    basis = {}

    def label(v):
        result = 0
        for i in sorted(basis, reverse=True):
            if v & (1 << i):
                row, code = basis[i]
                v ^= row
                result ^= code
        return v, result

    extra = []
    for column in columns:
        row, code = label(column)
        if row:
            basis[row.bit_length() - 1] = row, code ^ (1 << len(basis))
        else:
            extra.append(code)
    assert len(basis) == 24 and len(extra) == 12

    body_classes = []
    for pattern in range(8):
        marks_word = sum(((pattern >> i) & 1) << j for i, j in enumerate(marks))
        for word in parent["representatives_hex"]:
            remainder, code = label(syndrome(marks_word ^ int(word, 16)))
            assert remainder == 0
            body_classes.append(code)
    assert len(set(body_classes)) == 80 and len(certificate["extra_syndromes"]) == 4
    return np.array(body_classes, dtype=np.int32), extra


def leader_weights(extra):
    count = 1 << 24
    distances = np.zeros(count, dtype=np.uint8)
    for i in range(24):
        size = 1 << i
        distances[size:2 * size] = 1 + distances[:size]
    indices = np.arange(count, dtype=np.int32)
    for column in extra:
        distances = np.minimum(distances, 1 + distances[indices ^ column])
    return distances


def main():
    assert __debug__
    certificate = json.loads((SOURCE / "certificate.json").read_text())
    parent = json.loads((SOURCE / "parent39.json").read_text())
    classes, extra = quotient(certificate, parent)
    distance = leader_weights(extra)
    differences, counts = np.unique(classes[:, None] ^ classes, return_counts=True)
    assert len(differences) == 368 and counts.min() == 16

    # At Hamming distance at most five, any two central heads of length at most two conflict.
    # Edges with the same difference form a matching because the cosets in each copy are distinct.
    offsets = np.arange(len(distance), dtype=np.int32)
    for delta in differences:
        offsets = offsets[distance[offsets ^ delta] >= 6]
    assert len(offsets) == 64

    # Distance six requires antipodal heads of length exactly two. Each of the eight distinct
    # heads in one copy can match at most one of the eight heads in the other copy.
    cover_bounds = []
    for offset in offsets:
        critical = distance[offset ^ classes[:, None] ^ classes] == 6
        tier_counts = critical.reshape(8, 10, 8, 10).sum(axis=(1, 3))
        remaining = int(critical.sum() - tier_counts.max(axis=1).sum())
        max_degree = int(max(critical.sum(axis=0).max(), critical.sum(axis=1).max()))
        cover_bounds.append((remaining + max_degree - 1) // max_degree)
    assert min(cover_bounds) == 34
    max_fibres = 2 * (80 + len(certificate["extra_syndromes"])) - int(counts.min())
    assert max_fibres == 152
    print(json.dumps(dict(offsets=len(distance), exceptional_offsets=len(offsets),
                          min_exceptional_deletions=min(cover_bounds), maximum_fibres=max_fibres,
                          dense_point_upper_bound=max_fibres * (1 << 12)), indent=2))


if __name__ == "__main__":
    main()
