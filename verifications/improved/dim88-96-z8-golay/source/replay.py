"""Exact compressed-construction proof for the Z8 Golay sections in 88--96D.

The 2^37-word binary parent is due to Duursma--Greferath--Litsyn--Schmidt
(2001). No enumeration of that parent or its selected sections is claimed.
The minimum-distance proof reduces to a complete Z4 enumeration and 759
small coset exclusions. All extra support pairs are checked explicitly.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
bits = getattr(int, 'bit_count', lambda x: bin(x).count('1'))


def pc(a):
    a = np.asarray(a, dtype=np.uint64)
    a = a-((a >> np.uint64(1)) & np.uint64(0x5555555555555555))
    a = (a & np.uint64(0x3333333333333333))+((a >> np.uint64(2)) & np.uint64(0x3333333333333333))
    a = (a+(a >> np.uint64(4))) & np.uint64(0x0f0f0f0f0f0f0f0f)
    out = np.zeros(a.shape, dtype=np.uint8)
    for shift in range(0, 64, 8):
        out += ((a >> np.uint64(shift)) & np.uint64(255)).astype(np.uint8)
    return out


def binary_words(G):
    rows = [sum(int(x) << j for j, x in enumerate(row)) for row in G]
    words = np.zeros(2**len(rows), dtype=np.uint64)
    for j, row in enumerate(rows):
        words[2**j:2**(j+1)] = words[:2**j] ^ np.uint64(row)
    assert np.count_nonzero(words == 0) == 1
    return words


def binary_kernel(G):
    rows = [sum(int(x) << j for j, x in enumerate(row)) for row in G]
    rank = 0
    pivots = []
    for j in range(G.shape[1]):
        found = next((i for i in range(rank, len(rows)) if rows[i] >> j & 1), None)
        if found is None:
            continue
        rows[rank], rows[found] = rows[found], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> j & 1:
                rows[i] ^= rows[rank]
        pivots.append(j)
        rank += 1
    result = []
    for free in set(range(G.shape[1]))-set(pivots):
        word = 1 << free
        for row, pivot in zip(rows, pivots):
            word |= ((row >> free) & 1) << pivot
        result.append([(word >> j) & 1 for j in range(G.shape[1])])
    return np.array(result, dtype=np.uint8)


def ring_words(G, modulus):
    result = np.zeros((modulus**len(G), G.shape[1]), dtype=np.uint8)
    for i, row in enumerate(G):
        size = modulus**i
        for a in range(1, modulus):
            result[a*size:(a+1)*size] = (result[:size]+a*row) % modulus
    return result


def unit_rows(G):
    """Unit-pivot elimination over Z4; nonunit rows remain below the pivots."""
    M = np.asarray(G, dtype=np.int64).copy() % 4
    rank = 0
    pivots = []
    for j in range(M.shape[1]):
        nz = np.flatnonzero(M[rank:, j] % 2)
        if not len(nz):
            continue
        k = rank+int(nz[0])
        M[[rank, k]] = M[[k, rank]]
        M[rank] = M[rank]*int(M[rank, j]) % 4
        for i in range(len(M)):
            if i != rank:
                M[i] = (M[i]-M[i, j]*M[rank]) % 4
        pivots.append(j)
        rank += 1
    return M, pivots


def canonical(V, M, pivots):
    V, M = V.copy() % 4, M.copy()
    for i, j in enumerate(pivots):
        V = (V-V[:, [j]]*M[[i]]) % 4
    rank = len(pivots)
    even_pivots = []
    for j in range(M.shape[1]):
        nz = np.flatnonzero(M[rank:, j])
        if not len(nz):
            continue
        k = rank+int(nz[0])
        M[[rank, k]] = M[[k, rank]]
        assert M[rank, j] == 2
        for i in range(rank+1, len(M)):
            if M[i, j]:
                M[i] = (M[i]-M[rank]) % 4
        V = (V-(V[:, [j]]//2)*M[[rank]]) % 4
        even_pivots.append(j)
        rank += 1
    assert not np.any(M[rank:])
    return V, even_pivots


def small_lee_errors():
    rows = [np.zeros(16, dtype=np.int64)]
    for weight in (1, 2, 3):
        for support in combinations(range(16), weight):
            for values in product((1, 2, 3), repeat=weight):
                if sum(min(v, 4-v) for v in values) > 3:
                    continue
                row = np.zeros(16, dtype=np.int64)
                row[list(support)] = values
                rows.append(row)
    result = np.array(rows)
    assert len(result) == 5489
    return result


def verify_z8(data):
    G = np.array(data['generator'], dtype=np.int64)
    g = np.array(data['generator_polynomial_low_first'], dtype=np.int64)
    h = np.array(data['complementary_factor_low_first'], dtype=np.int64)
    assert G.shape == (12, 24) and np.all((G >= 0) & (G < 8))
    assert np.array_equal(np.convolve(g, h) % 8, np.ones(23, dtype=np.int64))
    expected = np.zeros((12, 24), dtype=np.int64)
    for i in range(12):
        expected[i, i:i+12] = g
    expected[:, 23] = -expected[:, :23].sum(axis=1) % 8
    assert np.array_equal(G, expected)
    assert not np.any(G @ G.T % 8)
    # This triangular unit block proves freeness and uniform projection on
    # the first m ring coordinates, for every m <=12.
    assert not np.any(np.tril(G[:, :12], -1)) and np.all(np.diag(G[:, :12]) % 2 == 1)
    grey = np.array([[(a >> 2 & 1) ^ ((a & 1)*x) ^ ((a >> 1 & 1)*y)
                      for x, y in product((0, 1), repeat=2)] for a in range(8)], dtype=np.uint8)
    homogeneous = np.array([0, 2, 2, 2, 4, 2, 2, 2])
    for a, b in product(range(8), repeat=2):
        assert np.count_nonzero(grey[a] != grey[b]) == homogeneous[(a-b) % 8]
    assert np.all(grey.sum(axis=1) % 2 == 0) and len(np.unique(grey, axis=0)) == 8
    G4 = (G % 4).astype(np.uint8)
    low = ring_words(G4[:8], 4)
    high = ring_words(G4[8:], 4)
    histogram = np.zeros(49, dtype=np.int64)
    lee = np.array([0, 1, 2, 1], dtype=np.uint8)
    for value in high:
        words = (low+value) % 4
        histogram += np.bincount(lee[words].sum(axis=1).astype(np.int64), minlength=49)
    assert histogram.sum() == 2**24 and histogram[0] == 1
    assert np.flatnonzero(histogram[1:])[0]+1 == 12
    residue = ring_words((G % 2).astype(np.uint8), 2)
    coeff = ring_words(np.eye(12, dtype=np.uint8), 2)
    weights = residue.sum(axis=1)
    assert set(map(int, weights)) == {0, 8, 12, 16, 24}
    assert np.count_nonzero(weights == 0) == 1
    octads = np.flatnonzero(weights == 8)
    assert len(octads) == 759
    errors = small_lee_errors()
    for index in octads:
        outside = np.flatnonzero(residue[index] == 0)
        u = coeff[index].astype(np.int64) @ G % 8
        assert not np.any(u[outside] % 2)
        target = u[outside]//2
        M, pivots = unit_rows(G4[:, outside])
        assert len(pivots) == 11
        values, even_pivots = canonical(np.vstack([target, errors]), M, pivots)
        assert len(even_pivots) == 1
        assert not np.any(np.all(values[1:] == values[0], axis=1))
    # Even residues give 2*C4 and homogeneous weight 2*Lee(C4)>=24.
    # Residues of weight >=12 already contribute >=24. For each octad,
    # the exclusions above force >=8 more units outside its eight odd cells.
    return G, dict(z8_words=2**36, binary_union_words=2**37, binary_distance=24,
                   z4_words_enumerated=2**24, z4_lee_histogram=histogram.tolist(),
                   binary_residue_histogram=np.bincount(weights.astype(int), minlength=25).tolist(),
                   octads_checked=759, low_lee_errors_per_octad=5489,
                   low_lee_coset_exclusions=759*5489,
                   generator_sha256=hashlib.sha256(G.astype(np.uint8).tobytes()).hexdigest(),
                   grey_map=grey.tolist(), independent_uniform_ring_coordinates=list(range(12)))


def binary_field_multiply(a, b, modulus):
    value = 0
    degree = modulus.bit_length()-1
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= modulus
    return value


def field_tables(q):
    add = np.empty((q, q), dtype=np.uint8)
    mul = add.copy()
    for a, b in product(range(q), repeat=2):
        if q in (4, 16):
            add[a, b] = a ^ b
            mul[a, b] = binary_field_multiply(a, b, 7 if q == 4 else 19)
        elif q == 9:
            add[a, b] = (a % 3+b % 3) % 3+3*((a//3+b//3) % 3)
            mul[a, b] = ((a % 3)*(b % 3)-(a//3)*(b//3)) % 3+3*(((a % 3)*(b//3)+(a//3)*(b % 3)) % 3)
        else:
            assert q in (11, 13)
            add[a, b], mul[a, b] = (a+b) % q, (a*b) % q
    return add, mul


def qwords(G, q, add, mul):
    words = np.zeros((q**len(G), G.shape[1]), dtype=np.uint8)
    for i, row in enumerate(G):
        for a in range(1, q):
            words[a*q**i:(a+1)*q**i] = add[words[:q**i], mul[a, row]]
    return words


def verify_supports(data):
    lo, hi = data['supports_low'], data['supports_high']
    assert lo.dtype == hi.dtype == np.uint64 and lo.ndim == hi.ndim == 1 and lo.shape == hi.shape
    assert np.all(hi < 2**32) and np.all(pc(lo)+pc(hi) == 24)
    assert len(np.unique(np.column_stack([lo, hi]), axis=0)) == len(lo)
    I = np.column_stack([((lo >> np.uint64(c)) & 1) for c in range(64)]+
                        [((hi >> np.uint64(c)) & 1) for c in range(32)]).astype(np.float32)
    maximum = 0
    with threadpool_limits(limits=4):
        for start in range(0, len(I), 512):
            stop = min(start+512, len(I))
            products = I[start:stop] @ I[:stop].T
            products[np.arange(stop-start), start+np.arange(stop-start)] = 0
            maximum = max(maximum, int(products.max()))
            assert np.all(products <= 12)
    G = data['quaternary_generator']
    assert G.shape == (9, 24) and np.all(G <= 3)
    add, mul = field_tables(4)
    words = qwords(G, 4, add, mul)
    weights = np.count_nonzero(words, axis=1)
    assert np.count_nonzero(weights == 0) == 1 and min(weights[1:]) == 12
    witness = data['dual_witnesses']
    assert witness.dtype == np.uint64 and witness.shape == (len(lo), 2) and np.all(witness < 2**48)
    occupied = np.count_nonzero(I.reshape(-1, 24, 4).sum(axis=2), axis=1)
    trace = np.array([0, 0, 1, 1], dtype=np.uint8)
    checks = 0
    for i, pair in enumerate(witness):
        used = np.zeros(24, dtype=bool)
        misses = 0
        for encoded in pair:
            if not encoded:
                continue
            dual = ((encoded >> (2*np.arange(24, dtype=np.uint64))) & 3).astype(np.uint8)
            active = dual != 0
            assert not np.any(used & active)
            used |= active
            assert not np.any(np.bitwise_xor.reduce(mul[G, dual], axis=1))
            parity = 0
            for r in np.flatnonzero(active):
                cells = np.flatnonzero(I[i, 4*r:4*r+4])
                assert len(cells)
                values = trace[mul[dual[r], cells]]
                assert np.all(values == values[0])
                parity ^= int(values[0])
            assert parity == 1
            misses += 1
        assert occupied[i]-misses <= 12
        checks += misses
    sign = data['sign_generator']
    assert sign.shape == (14, 24) and np.all(sign <= 1)
    signs = binary_words(sign)
    assert min(pc(signs[1:])) == 6
    dual = binary_words(binary_kernel(sign))
    assert min(pc(dual[1:])) == 8
    octad = int(dual[np.flatnonzero(pc(dual) == 8)[0]])
    retained = [0 if j % 2 else int(np.count_nonzero(pc(signs & np.uint64(octad if j == 8 else (1 << j)-1)) == j//2)) for j in range(9)]
    section95 = int(np.bincount(pc(lo & np.uint64(255)), minlength=9) @ np.array(retained))
    return signs, dict(supports=len(lo), support_pairs_checked=len(lo)*(len(lo)-1)//2,
                       maximum_support_intersection=maximum, disjoint_dual_checks=checks,
                       occupied_row_histogram=np.bincount(occupied, minlength=15).tolist(),
                       sign_words=len(signs), sign_distance=6, quaternary_words=len(words),
                       quaternary_distance=12, section95_support_points=section95)


def sparse_sections(data, signs):
    lookup = np.array([bits(i) for i in range(16)], dtype=np.int8)
    profiles = Counter(tuple(int(lookup[(int(s) >> (4*r)) & 15]) for r in range(8)) for s in data['supports_low'])
    table = {(start, w): w-2*lookup[(signs >> np.uint64(start)) & np.uint64((1 << w)-1)]
             for start in range(25) for w in range(min(4, 24-start)+1)}
    zero_counts, equal_counts = [0]*9, [0]*5
    for profile, multiplicity in profiles.items():
        at = 0
        zero = np.ones(len(signs), dtype=bool)
        equal = zero.copy()
        origin = None
        for r, weight in enumerate(profile):
            values = table[at, weight]
            at += weight
            if origin is None:
                origin = values
            zero &= values == 0
            equal &= values == origin
            zero_counts[r+1] += multiplicity*int(zero.sum())
            if r < 4:
                equal_counts[r+1] += multiplicity*int(equal.sum())
    return dict(occupancy_profiles=len(profiles), zero_counts=zero_counts, equal_counts=equal_counts)


def rs_layers():
    result = {}
    for q in (9, 11, 13, 16):
        add, mul = field_tables(q)
        G = np.ones((4, 6), dtype=np.uint8)
        for i in range(1, 4):
            G[i] = mul[G[i-1], np.arange(1, 7)]
        W = qwords(G, q, add, mul)
        weights = np.count_nonzero(W, axis=1)
        assert np.count_nonzero(weights == 0) == 1 and min(weights[1:]) == 3
        result[str(q)] = dict(supports=q**4, signs=32, points=q**4*32, distance=3)
        if q == 16:
            j = np.count_nonzero(W < np.array([2, 2, 1, 1, 1, 1]), axis=1)
            hist = np.bincount(j, minlength=7)
            result['section95'] = dict(selected_counts=[2, 2, 1, 1, 1, 1],
                histogram=hist.tolist(), points=int(hist @ np.array([32, 0, 16, 0, 12, 0, 20])))
    assert result['section95']['points'] == 1316444
    return result


def section_matrix(d):
    m = 96-d if d <= 92 else 97-d if d <= 94 else 2 if d == 95 else 0
    u = np.zeros((m, 96), dtype=np.int8)
    for i in range(m):
        u[i, 4*i:4*i+4] = 1
    if d <= 92:
        A = u
    elif d <= 94:
        A = u[1:]-u[0]
    elif d == 95:
        A = u[:1]+u[1:]
    else:
        A = u
    packed = [sum(int(v % 2) << j for j, v in enumerate(row)) for row in A]
    basis = {}
    for value in packed:
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    assert len(basis) == len(A) == 96-d
    return m, A


def bounds(support, sparse, rs):
    result = {}
    for d in range(88, 97):
        m, A = section_matrix(d)
        if d <= 92:
            dense = F(2**36*6**m, 8**m)
            transversal = 0
            extra = sparse['zero_counts'][m]
            q = 9 if m == 8 else 11 if m >= 5 else 13
            small, pure = rs[str(q)]['points'], 192+4*m
            section = 'each of the first m four-cell row sums is zero'
        elif d <= 94:
            dense = F(2**36*(6**m+2), 8**m)+F(2**37, 2**m)
            transversal = 2**(33-m)
            extra = sparse['equal_counts'][m]
            q = 13
            small, pure = rs[str(q)]['points'], 192+4*m
            section = 'the first m four-cell row sums are equal'
        elif d == 95:
            dense, transversal, extra = 70*2**30, 2**31, support['section95_support_points']
            small, pure = rs['section95']['points'], 176
            section = 'the sum of the first eight coordinates is zero'
        else:
            dense, transversal, extra = 2**37, 2**32, support['supports']*16384
            small, pure = 2**21, 192
            section = 'no section'
        assert F(dense).denominator == 1
        if d <= 94:
            assert 6*q <= 96-4*m
        terms = dict(dense=int(dense), transversal24=transversal, support24=extra,
                     support6=small, pure=pure)
        # All comparisons use norm squared 96 and threshold 48.
        assert 96-2*24 == 4*(24-2*6) == 4*12 == 48
        assert 16*3 == 4*2*6 == 48 and 16*(6-2*2) <= 48
        assert (4*96 <= 48**2) and (16*96 <= 48**2)
        assert 4*48 <= 48**2 and 16*48 <= 48**2  # A3/dense and A3/weight24
        result[str(d)] = dict(total=sum(terms.values()), terms=terms, section=section,
                             selected_four_cell_rows=m, constraint_rows=A.tolist())
    return result


def main():
    if not __debug__:
        raise SystemExit('Run without -O: exact assertions are required.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads((ROOT/'data/z8_generator.json').read_text())
    _, parent = verify_z8(raw)
    with np.load(ROOT/'data/supports.npz', allow_pickle=False) as z:
        data = dict(z)
    signs, support = verify_supports(data)
    sparse = sparse_sections(data, signs)
    rs = rs_layers()
    report = dict(status='PASS', representation='canonical Z4 cosets and exact coordinate products',
        z8=parent, support=support, sparse_sections=sparse, reed_solomon=rs,
        dimensions=bounds(support, sparse, rs),
        support_sha256=hashlib.sha256((ROOT/'data/supports.npz').read_bytes()).hexdigest(),
        scope='Explicit compressed constructions and finite arithmetic proofs; no enumeration of 2^37 binary words. The parent code is published prior work.')
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', dimensions={d:r['total'] for d,r in report['dimensions'].items()})))


if __name__ == '__main__':
    main()
