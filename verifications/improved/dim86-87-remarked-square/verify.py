#!/usr/bin/env python3
"""Exact replay for the re-marked 86–87D axial Pless, square and D_n union."""
import hashlib
import importlib.util
import itertools as it
import json
import math
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
IMP = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def quadruples():
    # The 91 data blocks are the 14-point La Jolla design (as relabeled by Sage).
    B = [tuple(map(int, x.split())) for x in (HERE / 'sqs14.txt').read_text().splitlines()]
    assert len(B) == 91 and all(0 <= a < b < c < d < 14 for a, b, c, d in B)
    small_triples = [t for q in B for t in it.combinations(q, 3)]
    assert len(small_triples) == len(set(small_triples)) == math.comb(14, 3)
    # Hanani's 3n-2 construction, then delete infinity.
    def r(k, x):
        return (k % 3) * 13 + x
    C = []
    for a, b, c, d in B:
        if d < 13:
            for i, j, k in it.product(range(3), repeat=3):
                C.append((r(i, a), r(j, b), r(k, c), r(-i-j-k, d)))
        else:
            for i, j in it.product(range(3), repeat=2):
                C.append((r(i, a), r(j, b), r(-i-j, c), 39))
            for x, y, z in ((a, b, c), (a, c, b), (b, c, a)):
                for i in range(3):
                    C.append((r(i, x), r(i, y), r(i+1, z), r(i+2, z)))
    for a, b in it.combinations(range(13), 2):
        for i in range(3):
            C.append((r(i, a), r(i, b), r(i+1, a), r(i+1, b)))
    for a in range(13):
        C.append((r(0, a), r(1, a), r(2, a), 39))
    C = [tuple(sorted(q)) for q in C]
    assert len(C) == len(set(C)) == math.comb(40, 3) // 4
    assert all(0 <= a < b < c < d < 40 for a, b, c, d in C)
    triples = [t for q in C for t in it.combinations(q, 3)]
    assert len(triples) == len(set(triples)) == math.comb(40, 3)
    D = [q for q in C if q[-1] < 39]
    assert len(D) == 2223
    # At least one uncovered triple on every pair in any order-39 packing;
    # therefore this support count is optimal among triple-disjoint quadruples.
    assert len(D) == (math.comb(39, 3) - math.comb(39, 2) // 3) // 4
    return D


def square_count():
    sups = quadruples()
    signs = [s for s in it.product((-1, 1), repeat=8) if math.prod(s) == 1]
    assert len(set(signs)) == len(signs) == 128
    assert max(sum(x*y for x, y in zip(a, b)) for a, b in it.combinations(signs, 2)) == 4
    # Dot with sparse <=4 sqrt(3) cos(15 degrees); squaring gives 24+12 sqrt(3)<45.
    assert 12**2 * 3 < (45 - 24)**2
    assert 4**2 < 45  # square / dense
    assert Q(2)**2 / 2 == Q(2*4, 4)  # root / square
    assert 2*(2*3) < Q(2*45, 4)  # root / prior dense or sparse (Cauchy)
    assert 15-9 > 4  # distinguish squares, sparse and dense by support
    return len(sups) * len(signs)


def main():
    assert __debug__, 'Run without -O'
    axial = load('axial_parent', IMP / 'dim86-87-axial-pless' / 'axial_pless.py')
    chord = load('chord_parent', IMP / 'dim86-87-ternary-chord-sections' / 'verify.py')
    assert hashlib.sha256(chord.QR.read_bytes()).hexdigest() == chord.SOURCES[chord.QR]
    G = chord.np.asarray(json.loads(chord.QR.read_text())['generator_Q48'], dtype='i8')
    assert G.shape == (24, 48) and set(G.flat) == {0, 1, 2}
    assert chord.rank3(G) == 24 and (G @ G.T % 3 == 0).all()
    norms = 39 + 8 * Q(3, 4)
    assert norms == 45 and axial.DENSE == 2 * 2**8 * 3**(24-9)
    assert all(-norms/2 <= norms-Q(3,2)*d <= norms/2 for d in range(15,48,3))
    axial.verify_geometry()
    phases = axial.phase_counts()
    squares = square_count()
    assert squares == 284544
    outer = axial.build_outer()
    rows = json.loads((HERE / 'claims.json').read_text())
    for d, data, result in zip((86, 87), outer, rows):
        axes, hist = result['axes'], result['histogram']
        marked = set(axes) | {0}
        assert len(marked) == 9 and (0 in axes) == (d == 87)
        assert chord.rank3(G[:,sorted(marked)]) == 9 and 2*39+len(axes) == d
        counts = axial.census(data, axes)
        assert counts.tolist() == hist
        sparse = int(counts @ phases[:len(counts)])
        roots = 4 * math.comb(d, 2)
        # Normal form for D_d: any distinct signed two-subsets have dot <= 1.
        # On equal supports a distinct sign negates at least one contribution;
        # on different supports there is at most one coordinate in common.
        assert result['dense'] == axial.DENSE and result['sparse'] == sparse
        assert result['squares'] == squares and result['root_tail'] == roots
        assert result['lower_bound'] == axial.DENSE + sparse + squares + roots
        assert result['gain'] == result['lower_bound'] - result['previous']
        print(json.dumps(result))


if __name__ == '__main__':
    main()
