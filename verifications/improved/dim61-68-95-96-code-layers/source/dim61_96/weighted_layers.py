#!/usr/bin/env python3
"""Explicit weighted dense codes with sparse layers in dimensions 62, 63, 95.

All vectors have squared norm N=4D.  A binary distance-D code of length
n=N-u is mapped to signs of (1,...,1,sqrt(u+1)).  Its inner products are
at most N-2D=N/2.  Weight-D sparse vectors with entries +-2 avoid the
last coordinate, so their cross inner products are at most 2D=N/2.
This removes the usual requirement n >= 4D for adjoining that layer.

The program constructs every input code and verifies its relevant minimum
distance exhaustively, or inherits it through a proved concatenation and
shortening.  Pairwise geometric bounds are checked with integer arithmetic.
No world-record or priority claim is made.
"""

import hashlib
import itertools
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dim49plus'))
from binary_62_code import (build_code, build_code64_rows, impose_zeros,
                            keep_coordinates, popcounts64, rref,
                            weight_enumerator)

SOURCE = 'https://www.codetables.de/BKLC/BKLC.php?q=2&n=95&k=32'


def matrix(rows, n):
    return np.array([[(x >> c) & 1 for c in range(n)] for x in rows], dtype=np.uint8)


def packed(G):
    return [sum(int(v) << c for c, v in enumerate(row)) for row in G]


def all_binary_words(G):
    words = np.zeros(1 << len(G), dtype=np.uint64)
    for i, row in enumerate(packed(G)):
        words[1 << i:1 << (i + 1)] = words[:1 << i] ^ np.uint64(row)
    return words


def collect_min_words(G, weight):
    rows = packed(G)
    low_bits = min(18, len(rows))
    low = all_binary_words(G[:low_bits])
    high = [0]
    for row in rows[low_bits:]:
        high += [x ^ row for x in high]
    selected = []
    histogram = np.zeros(G.shape[1] + 1, dtype=np.int64)
    for x in high:
        words = low ^ np.uint64(x)
        weights = popcounts64(words)
        histogram += np.bincount(weights, minlength=len(histogram))
        selected.append(words[weights == weight])
    assert histogram.sum() == 1 << len(G) and histogram[0] == 1
    assert np.flatnonzero(histogram[1:])[0] + 1 == weight
    return np.concatenate(selected), histogram


def kernel(checks, n):
    checks, pivots = rref(checks, n)
    free = [i for i in range(n) if i not in pivots]
    return [sum((((row >> f) & 1) << p) for row, p in zip(checks, pivots)) | (1 << f)
            for f in free]


def sqs64():
    # Affine 2-planes in F_2^6: every triple lies in its unique fourth-point block.
    blocks = {tuple(sorted((a, b, c, a ^ b ^ c)))
              for a, b, c in itertools.combinations(range(64), 3)}
    assert len(blocks) == 10416
    triples = [t for block in blocks for t in itertools.combinations(block, 3)]
    assert len(triples) == len(set(triples)) == 41664
    return sorted(blocks)


def field_mul(a, b, degree, modulus):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= modulus
    return result


M4 = np.array([[field_mul(a, b, 2, 0b111) for b in range(4)]
               for a in range(4)], dtype=np.uint8)


def gf4_rref(rows):
    rows = np.array(rows, dtype=np.uint8, copy=True)
    rank, pivots = 0, []
    for c in range(rows.shape[1]):
        nz = np.flatnonzero(rows[rank:, c])
        if not len(nz):
            continue
        k = rank + int(nz[0])
        rows[[rank, k]] = rows[[k, rank]]
        inverse = next(a for a in range(1, 4) if M4[rows[rank, c], a] == 1)
        rows[rank] = M4[inverse, rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i, c]:
                rows[i] ^= M4[rows[i, c], rows[rank]]
        pivots.append(c)
        rank += 1
        if rank == len(rows):
            break
    assert not np.any(rows[rank:])
    return rows[:rank], pivots


def gf4_shorten(G, coordinates):
    rows = np.array(G, dtype=np.uint8, copy=True)
    for c in coordinates:
        nz = np.flatnonzero(rows[:, c])
        if not len(nz):
            continue
        k = int(nz[0])
        pivot = rows[k].copy()
        inverse = next(a for a in range(1, 4) if M4[pivot[c], a] == 1)
        pivot = M4[inverse, pivot]
        rows = np.delete(rows, k, axis=0)
        for i in range(len(rows)):
            rows[i] ^= M4[rows[i, c], pivot]
    assert not np.any(rows[:, coordinates])
    keep = [c for c in range(G.shape[1]) if c not in coordinates]
    return gf4_rref(rows[:, keep])[0]


def quaternary24():
    # ord_29(2)=28, so Phi_29 is irreducible over F_2.  Its root alpha has
    # order 29 in F_(2^28).  The quadratic-residue Gaussian period eta obeys
    # eta^2=eta+1 and identifies the F_4 subfield without an external library.
    modulus = (1 << 29) - 1
    mul = lambda a, b: field_mul(a, b, 28, modulus)
    alpha_powers, x = [], 1
    for _ in range(29):
        alpha_powers.append(x)
        x = mul(x, 2)
    assert x == 1 and len(set(alpha_powers)) == 29
    assert len({pow(2, i, 29) for i in range(28)}) == 28
    residues = sorted({(i * i) % 29 for i in range(1, 29)})
    eta = 0
    for r in residues:
        eta ^= alpha_powers[r]
    assert mul(eta, eta) == (eta ^ 1)
    coefficients = [1]
    for r in residues:
        new = [0] * (len(coefficients) + 1)
        for j, v in enumerate(coefficients):
            new[j] ^= mul(v, alpha_powers[r])
            new[j + 1] ^= v
        coefficients = new
    subfield = {0: 0, 1: 1, eta: 2, eta ^ 1: 3}
    g = [subfield[x] for x in coefficients]
    G29 = np.zeros((15, 29), dtype=np.uint8)
    for i in range(15):
        G29[i, i:i + 15] = g
    G30 = np.column_stack((G29, np.bitwise_xor.reduce(G29, axis=1)))
    G24 = gf4_shorten(G30, list(range(24, 30)))
    assert G24.shape == (9, 24)
    words = np.zeros((4 ** 9, 24), dtype=np.uint8)
    for i, row in enumerate(G24):
        size = 4 ** i
        for a in (1, 2, 3):
            words[a * size:(a + 1) * size] = words[:size] ^ M4[a, row]
    weights = np.count_nonzero(words, axis=1)
    hist = np.bincount(weights, minlength=25)
    assert hist[0] == 1 and np.flatnonzero(hist[1:])[0] + 1 == 12
    pair = next((a, b) for a, b in itertools.combinations(range(24), 2)
                if len(gf4_rref(G24[:, [a, b]].T)[0]) == 2)
    survivors = (words[:, pair[0]] != 0) & (words[:, pair[1]] != 0)
    assert int(survivors.sum()) == 147456
    return G24, {'generator_polynomial_GF4_low_first': g,
                 'weight_histogram': {str(i): int(v) for i, v in enumerate(hist) if v},
                 'excluded_zero_cells_rows': list(pair),
                 'retained_supports_on_94_cells': int(survivors.sum())}


def sign24():
    # [23,14,5] stored matrix on Grassl's page, extended by binary parity.
    # Preserve the full matrix: its pivots are not the first fourteen columns.
    G23 = np.array(json.loads((HERE / 'source_sign23.json').read_text()), dtype=np.uint8)
    assert G23.shape == (14, 23) and len(rref(packed(G23), 23)[0]) == 14
    G24 = np.column_stack((G23, np.bitwise_xor.reduce(G23, axis=1)))
    weights = popcounts64(all_binary_words(G24))
    hist = np.bincount(weights, minlength=25)
    assert hist[0] == 1 and np.flatnonzero(hist[1:])[0] + 1 == 6
    return G24, {str(i): int(v) for i, v in enumerate(hist) if v}


def dense96(Gbinary24, Gquaternary24):
    # Nested inner codes <1111> subset even(4) subset F_2^4.
    # Binary outer code: d=6; quaternary outer code: d=12;
    # final binary repetition outer code: d=24.  Distance is >= min(4*6,2*12,24).
    rows = [sum(int(v) * 15 << (4 * i) for i, v in enumerate(row)) for row in Gbinary24]
    for row in Gquaternary24:
        for scalar in (1, 2):
            symbols = M4[scalar, row]
            rows.append(sum((((int(s) & 1) * 3) ^ (((int(s) >> 1) & 1) * 5)) << (4 * i)
                            for i, s in enumerate(symbols)))
    rows.append(sum(1 << (4 * i) for i in range(24)))
    assert len(rref(rows, 96)[0]) == 33
    return matrix(rref(rows, 96)[0], 96)


def dense95(Gbinary24, Gquaternary24):
    rows = impose_zeros(packed(dense96(Gbinary24, Gquaternary24)), [95])
    assert len(rows) == 32 and all(x < (1 << 95) for x in rows)
    return matrix(rref(rows, 95)[0], 95)


def rs6():
    M16 = np.array([[field_mul(a, b, 4, 0b10011) for b in range(16)]
                   for a in range(16)], dtype=np.uint8)
    G = np.zeros((4, 6), dtype=np.uint8)
    G[0] = 1
    for i in range(1, 4):
        G[i] = M16[G[i - 1], np.arange(1, 7)]
    words = np.zeros((16 ** 4, 6), dtype=np.uint8)
    for i, row in enumerate(G):
        size = 16 ** i
        for a in range(1, 16):
            words[a * size:(a + 1) * size] = words[:size] ^ M16[a, row]
    hist = np.bincount(np.count_nonzero(words, axis=1), minlength=7)
    assert hist[0] == 1 and np.flatnonzero(hist[1:])[0] + 1 == 3
    retained = int(np.count_nonzero(words[:, 0]))
    assert retained == 61440
    return G, retained


def geometric_checks(N, u, w, min_sign_distance, lower_weight):
    D, n = N // 4, N - u
    assert N % 4 == 0 and n - 1 >= D
    checks = {
        'dense_norm_squared': (n - 1) + (u + 1) == N,
        'dense_pair_ip': N - 2 * D <= N // 2,
        'dense_sparse_ip': 2 * D <= N // 2,
        'weightD_different_support_ip': 4 * (D // 2) <= N // 2,
        'weightD_same_support_ip': 4 * (D - 2 * min_sign_distance) <= N // 2,
        'sparse_two_levels_ip': 4 * N * lower_weight <= (N // 2) ** 2,
        'pole_dense_ip_squared': N * (u + 1) <= (N // 2) ** 2,
        'pole_weightD_ip_squared': 4 * N <= (N // 2) ** 2,
        'pole_lower_ip_squared': N * N <= lower_weight * (N // 2) ** 2,
        # sqrt(N/lower_weight)*(lower_weight-1+sqrt(u+1)) <= N/2.
        # Here N/lower_weight=16 in both applications, so no irrational test remains.
        'dense_lower_ip': N // lower_weight == 16 and
            (N // 8 - lower_weight + 1) >= 0 and
            u + 1 <= (N // 8 - lower_weight + 1) ** 2,
        'lower_different_support_ip': 16 * (lower_weight // 2) <= N // 2,
        'lower_same_support_ip': 16 * (lower_weight - 2 * (1 if lower_weight == 4 else 2)) <= N // 2,
    }
    assert all(checks.values()), checks
    return checks


def main():
    start = perf_counter()
    G64 = matrix(build_code64_rows(), 64)
    hist64 = weight_enumerator(G64)
    assert np.flatnonzero(hist64[1:])[0] + 1 == 16
    rows63 = impose_zeros(packed(G64), [63])
    assert len(rows63) == 27
    G63 = matrix(keep_coordinates(rows63, list(range(63))), 63)
    G62 = build_code()
    min_words, hist62 = collect_min_words(G62, 16)
    assert len(min_words) == 17249
    occurrences = np.array([np.count_nonzero((min_words >> np.uint64(c)) & np.uint64(1))
                            for c in range(62)])
    best_coordinate = int(np.argmin(occurrences))
    retained61 = int(len(min_words) - occurrences[best_coordinate])
    assert retained61 >= 12798
    checks16 = [sum(1 << j for j in range(16))]
    checks16 += [sum(((j >> i) & 1) << j for j in range(16)) for i in range(4)]
    sign16 = matrix(kernel(checks16, 16), 16)
    assert sign16.shape == (11, 16)
    assert np.min(popcounts64(all_binary_words(sign16))[1:]) == 4
    sqs = sqs64()
    results = {}
    for n, dense_G, supports in [(62, G62, retained61), (63, G63, len(min_words))]:
        blocks = [b for b in sqs if max(b) < n]
        assert len(blocks) == {62: 9145, 63: 9765}[n]
        terms = {'dense': 1 << len(dense_G), 'weight16': supports * 2048,
                 'weight4': len(blocks) * 16, 'poles': 2 * n}
        baseline = {62: 71310732, 63: 138419844}[n]
        results[str(n)] = {'total': sum(terms.values()), 'terms': terms,
                           'repository_baseline': baseline,
                           'gain': sum(terms.values()) - baseline,
                           'norm_squared': 64, 'heavy_coordinate_squared': 65 - n,
                           'support_count_weight16': supports,
                           'checks': geometric_checks(64, 64 - n, 16, 4, 4)}
    Gq4, q4_report = quaternary24()
    Gsign24, sign24_hist = sign24()
    G95 = dense95(Gsign24, Gq4)
    Grs6, support6 = rs6()
    # The even-weight length-six sign code has dimension five and distance two.
    sign6 = matrix([(1 << i) | (1 << 5) for i in range(5)], 6)
    assert np.min(popcounts64(all_binary_words(sign6))[1:]) == 2
    terms95 = {'dense': 1 << 32, 'weight24': 147456 * 16384,
               'weight6': support6 * 32, 'poles': 190}
    results['95'] = {'total': sum(terms95.values()), 'terms': terms95,
                     'repository_baseline': 6473065046,
                     'gain': sum(terms95.values()) - 6473065046,
                     'norm_squared': 96, 'heavy_coordinate_squared': 2,
                     'checks': geometric_checks(96, 1, 24, 6, 6)}
    artifact = HERE / 'weighted_layers.npz'
    np.savez_compressed(artifact, dense_generator62=G62, dense_generator63=G63,
                        dense_generator95=G95, support_weight16_on62=min_words,
                        deleted_coordinate_for_support61=best_coordinate,
                        sign_generator16=sign16, sign_generator24=Gsign24,
                        support_generator_GF4_24=Gq4, support_generator_GF16_6=Grs6,
                        sign_generator6=sign6, sqs64=np.array(sqs, dtype=np.uint8))
    report = {'status': 'Explicit input codes and exact layer geometry verified',
              'construction': 'One weighted dense coordinate; weight-D sparse layer avoids it',
              'dimensions': results,
              'binary64_weight_histogram': {str(i): int(v) for i, v in enumerate(hist64) if v},
              'binary62_weight_histogram': {str(i): int(v) for i, v in enumerate(hist62) if v},
              'extra_shortening_coordinate_for_support61': best_coordinate,
              'support61_count': retained61,
              'quaternary24': q4_report, 'binary_sign24_weight_histogram': sign24_hist,
              'distance95_proof': 'Nested inner distances 4,2,1 times outer distances 6,12,24; shorten one coordinate',
              'code_source': SOURCE,
              'generator95_sha256': hashlib.sha256(G95.tobytes()).hexdigest(),
              'artifact_sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
              'seconds': perf_counter() - start}
    (HERE / 'weighted_layers.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
