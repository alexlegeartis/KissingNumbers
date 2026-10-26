"""Exact, standard-library verification of layered_kissing_v2 certificates.

All comparisons are integral. Antipodal lines denote both signs, so absolute
dot products are used. No claim of optimality is made.
"""
import itertools
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
pc = lambda x: bin(x).count('1')


def verify(path):
    started = time.monotonic()
    c = json.loads(path.read_text())
    assert c['format'] == 'layered_kissing_v2'
    n = c['dimension']
    supports = c['full_parity_supports']
    assert len({s for s, p in supports}) == len(supports)
    for s, p in supports:
        assert 0 <= s < 1 << n and pc(s) == 8 and p in (0, 1)
    for (s, _), (t, _) in itertools.combinations(supports, 2):
        assert pc(s & t) <= 4
    # Same-support parity words differ in at least two signs, so their
    # unscaled dot is <=4. Different supports have <=4 shared positions.
    coordinates = c['dense_coordinates']
    assert len(coordinates) == len(set(coordinates)) == 32
    assert all(0 <= j < n for j in coordinates)
    generators = c['dense_generator']
    words = [0]
    for generator in generators:
        assert 0 <= generator < 1 << 32
        words += [w ^ generator for w in words]
    assert len(words) == len(set(words)) == 131072
    min_distance = min(pc(w) for w in words[1:])
    assert min_distance >= 8
    shift = c['dense_coset_shift']
    assert 0 <= shift < 1 << 32
    # Dense/middle <=16, root/middle <=16, dense/root <=8, root/root <=16
    # from their coordinate shapes; all these base vectors have norm^2 32.
    removed = [tuple(r) for r in c['removed_root_lines']]
    assert len(removed) == len(set(removed))
    for a, b, sign in removed:
        assert 0 <= a < b < n and sign in (-1, 1)
    roots = [(a, b, sign) for a, b in itertools.combinations(range(n), 2)
             for sign in (-1, 1) if (a, b, sign) not in removed]
    heads = c['head_antipodal_lines']
    assert all(len(h) == n and all(type(x) is int for x in h) for h in heads)
    norms = [sum(x*x for x in h) for h in heads]
    assert all(norm > 0 for norm in norms)
    maxima = dict(dense=0, middle=0, roots=0, heads=0)
    for h, norm in zip(heads, norms):
        for a, b, sign in roots:
            dot = abs(4*(h[a]+sign*h[b]))
            assert 4*dot*dot <= 32*norm
            maxima['roots'] = max(maxima['roots'], dot)
        for mask, parity in supports:
            values = [h[j] for j in range(n) if mask >> j & 1]
            dot = sum(abs(x) for x in values)
            if sum(x < 0 for x in values) % 2 != parity:
                dot -= 2*min(abs(x) for x in values)
            # Weight eight and fixed parity are invariant under negation.
            dot *= 2
            assert 4*dot*dot <= 32*norm
            maxima['middle'] = max(maxima['middle'], dot)
        projected = [h[j] for j in coordinates]
        assert all(abs(x) == 1 for x in projected)
        minus = sum((x < 0) << j for j, x in enumerate(projected))
        maximum = max(abs(32-2*pc(w ^ shift ^ minus)) for w in words)
        assert 4*maximum*maximum <= 32*norm
        maxima['dense'] = max(maxima['dense'], maximum)
    for i, j in itertools.combinations(range(len(heads)), 2):
        dot = abs(sum(a*b for a, b in zip(heads[i], heads[j])))
        assert 4*dot*dot <= norms[i]*norms[j]
        maxima['heads'] = max(maxima['heads'], dot)
    points = len(words)+128*len(supports)+2*len(roots)+2*len(heads)
    assert points == c['points']
    return dict(dimension=n, points=points, dense_points=len(words),
                middle_points=128*len(supports), root_points=2*len(roots),
                head_points=2*len(heads), head_squared_norms=sorted(set(norms)),
                dense_min_distance=min_distance, maximum_absolute_dots=maxima,
                seconds=round(time.monotonic()-started, 3), verified=True)


if __name__ == '__main__':
    records = [verify(path) for path in sorted(ROOT.glob('certificate*.json'))]
    (ROOT/'verification.json').write_text(json.dumps(records, indent=2)+'\n')
    print(json.dumps(records, indent=2))
