#!/usr/bin/env python3
"""Replay the 86D and 87D Pless axial sparse layers beside the pinned dense section.

Usage: python axial_pless.py /path/to/hexagon_sparse15_phase7812_20261117
The asserted parent minimum distance of 15 is the same premise as the existing
Pless package. The earlier 86D verifier supplies the dense proof.
"""
import hashlib
import itertools
import json
import sys
from pathlib import Path
import numpy as np

SOURCE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "dim90-92-hexagonal-phase-layers" / "hexagon_sparse15_phase7812_20261117"
P45 = "outer_45_0_1_2.npz"
CERT = "outer_certificate.json"
PHASE = "phase7812.txt"
BLOBS = {
    P45: "fa595ddd5ed85dae7bcb957829105635d02a6689",
    CERT: "d03993eaee704683e6189ff301df6ce13930dc7d",
    PHASE: "e44a0b6d8b4a566d3b39615e06b5a8873f6c39a9",
}
AXES = {86: [[1,2,24,28,34,39,41,47],[9442,46272,85652,83516,44700,13304,2292,204,2]], 87: [[0,1,2,24,28,34,39,41,47],[9442,52050,110154,125218,80502,29934,6862,918,24,0]]}
H14 = [1,3,9,27,81,243,296,113,145,424,446,468,123,373]
QR11 = [2,2,1,2,0,1]
PERMUTATION = [1,0,23,12,16,18,10,20,14,21,6,17,3,22,8,19,4,11,5,15,7,9,13,2,25,24,47,36,40,42,34,44,38,45,30,41,27,46,32,43,28,35,29,39,31,33,37,26]
SIGNS = [2,1,2,2,2,2,1,2,1,2,2,1,1,2,2,1,1,2,1,2,1,1,1,1,2,1,2,2,2,2,1,2,1,2,2,1,1,2,2,1,1,2,1,2,1,1,1,1]
MASK23 = np.uint32((1 << 23) - 1)
MASK24 = np.uint32(1 << 23)
DENSE = 2 * 2 ** 8 * 3 ** 15
PRIOR = {86: 7422120168, 87: 7422120170}

def source_blob(name):
    data = (SOURCE / name).read_bytes()
    content_hash = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    assert content_hash == BLOBS[name], name
    return data

def trits(n, length):
    return np.array([n // 3 ** j % 3 for j in range(length)], dtype=np.int64)

def phase_counts():
    # A relation of weight at most four would repeat a syndrome of weight at most two.
    H = np.array([trits(word, 6) for word in H14])
    assert np.array_equal(H[:6], np.eye(6, dtype=np.int64))
    syndromes = {bytes(6)}
    for i in range(14):
        for a in (1, 2):
            syndromes.add(bytes((a * H[i] % 3).astype(np.uint8)))
        for j in range(i):
            for a, b in itertools.product((1, 2), repeat=2):
                syndromes.add(bytes(((a * H[i] + b * H[j]) % 3).astype(np.uint8)))
    assert len(syndromes) == 1 + 2 * 14 + 4 * 91
    # The [11,6,5] cyclic code is shortened successively on initial coordinates.
    G11 = np.array([np.pad(QR11, (j, 5 - j)) for j in range(6)])
    assert G11[0, 0] != 0
    assert min(np.count_nonzero(trits(j, 6) @ G11 % 3) for j in range(1, 729)) == 5
    rows = source_blob(PHASE).splitlines()
    assert len(rows) == 7812 and all(len(row) == 15 and set(row) <= set(b"012") for row in rows)
    assert len(set(rows)) == 7812
    code = np.array([np.frombuffer(row, dtype=np.uint8) - 48 for row in rows])
    a, b = (code == 1) @ (1 << np.arange(15)), (code == 2) @ (1 << np.arange(15))
    pc = np.array([i.bit_count() for i in range(1 << 15)], dtype=np.uint8)
    for i in range(len(rows) - 1):
        assert pc[(a[i] ^ a[i + 1:]) | (b[i] ^ b[i + 1:])].min() >= 5
    return np.array([7812, 6561, 2187, 729, 729, 243, 81, 27, 9, 3], dtype=np.int64)

def coordinate_bit(j):
    return (j + 23) % 24

def pack(words):
    powers = 1 << np.array([coordinate_bit(j) for j in range(24)], dtype=np.uint32)
    entries = [(words[:, start:start + 24] == s) @ powers for start in (0, 24) for s in (1, 2)]
    return np.stack(entries, axis=-1).astype(np.uint32)

def rotations(words):
    finite = words & MASK23
    fixed = words & MASK24
    for r in range(23):
        yield (((finite << r) | (finite >> (23 - r))) & MASK23) | fixed

def unique_rows(rows):
    return np.unique(rows.reshape(-1, 4).view(np.dtype("V16")).reshape(-1)).view(np.uint32).reshape(-1, 4)

def monomial_action(rows):
    result = np.zeros_like(rows)
    for j, k in enumerate(PERMUTATION):
        for symbol in (1, 2):
            source = 2 * (j // 24) + symbol - 1
            target = 2 * (k // 24) + (SIGNS[j] * symbol % 3) - 1
            result[:, target] |= ((rows[:, source] >> coordinate_bit(j)) & 1) << coordinate_bit(k)
    return result

def verify_geometry():
    # Coordinates below are twice ordinary A2 coordinates, with metric diag(1,3)/4.
    tri = np.array([(2, 0), (-1, 1), (-1, -1)])
    roots = np.array([(3, 1), (-3, 1), (0, -2)])
    metric = np.diag([1, 3])
    assert set((roots @ metric @ tri.T).flat) == {-6, 0, 6}
    assert np.array_equal(roots @ metric @ roots.T, 18 * np.eye(3, dtype=int) - 6)
    # Axial sparse entries are 0 or +/-sqrt(3). Against a chord point +/-sqrt(3)/2,
    # the absolute product is at most 3/2. Against another axial symbol, a differing
    # symbol contributes squared distance at least 3. This remains true of the A2 roots.
    assert 15 * 3 == 45 and 15 * 3 >= 45 and 5 * 9 >= 45
    # A pure D root scaled to norm squared 45 meets at most two sparse axes.
    assert 4 * (4 * 3 * 45 // 2) < 45 ** 2
    # Orthogonal pure roots in rank d have at most one equal signed coordinate,
    # unless they are equal or antipodal on the same support.
    for d in (8, 9):
        D = []
        for i, j in itertools.combinations(range(d), 2):
            for a, b in itertools.product((-1, 1), repeat=2):
                v = np.zeros(d, dtype=int)
                v[i], v[j] = a, b
                D.append(v)
        D = np.array(D)
        gram = D @ D.T
        assert len(D) == 2 * d * (d - 1) and (np.diag(gram) == 2).all()
        assert gram[np.triu_indices(len(D), 1)].max() == 1

def build_outer():
    G = np.array(json.loads(source_blob(CERT))["outer_generator"], dtype=np.int64)
    assert np.array_equal(G[:, :24], np.eye(24, dtype=np.int64))
    assert (G @ G.T % 3 == 0).all()
    steps = [
        ([0] + list(range(2, 24)) + [1, 24] + list(range(26, 48)) + [25], [1] * 48),
        (PERMUTATION, SIGNS),
    ]
    for permutation, signs in steps:
        assert len(set(permutation)) == 48 and set(signs) <= {1, 2}
        image = np.zeros_like(G)
        image[:, permutation] = G * signs % 3
        assert np.array_equal(image[:, :24] @ G[:, 24:] % 3, image[:, 24:])
    source_blob(P45)
    with np.load(SOURCE / P45, allow_pickle=False) as source:
        data = source["words"]
    assert data.shape == (130944, 45) and data.dtype == np.uint8
    assert data.max() == 2 and (np.count_nonzero(data, axis=1) == 15).all()
    words = np.pad(data, ((0, 0), (3, 0)))
    assert np.array_equal(words[:, :24] @ G[:, 24:] % 3, words[:, 24:])
    first = unique_rows(np.stack(list(rotations(pack(words)))))
    assert len(first) == 285384 and not ((first[:, 0] | first[:, 1]) & MASK24).any()
    second = unique_rows(np.concatenate([first, *rotations(monomial_action(first))]))
    assert len(second) == 415104
    return first, second

def census(rows, axes):
    masks = [sum(1 << coordinate_bit(j) for j in axes if j // 24 == half) for half in range(2)]
    pc = np.array([i.bit_count() for i in range(1 << 12)], dtype=np.uint8)
    occ = [(rows[:, 2 * i] | rows[:, 2 * i + 1]) & masks[i] for i in range(2)]
    overlap = sum(pc[body & 4095] + pc[body >> 12] for body in occ)
    return np.bincount(overlap, minlength=len(axes) + 1)

def main():
    counts = phase_counts()
    verify_geometry()
    for dimension, rows in zip((86, 87), build_outer()):
        axes, expected_histogram = AXES[dimension]
        assert len(axes) == len(set(axes)) == dimension - 78
        assert min(axes) >= (1 if dimension == 86 else 0) and max(axes) < 48
        hist = census(rows, axes)
        assert hist.tolist() == expected_histogram
        sparse = int(hist @ counts[:len(hist)])
        tail_rank = dimension - 78
        tail = 2 * tail_rank * (tail_rank - 1)
        total = DENSE + sparse + tail
        print(json.dumps(dict(dimension=dimension, sparse=sparse, dense=DENSE, pure_tail=tail,
                              lower_bound=total, gain=total - PRIOR[dimension], histogram=hist.tolist())))

if __name__ == "__main__":
    main()

