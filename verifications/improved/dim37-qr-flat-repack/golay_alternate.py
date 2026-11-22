#!/usr/bin/env python3
"""Independent 576-point construction on the same 12-coordinate 37D flat."""
from collections import Counter
from itertools import combinations

stencil = sum(1 << i for i in (0, 1, 5, 6, 7, 9, 11))
code = {0}
for i in range(12):
    row = stencil << i
    row |= (row.bit_count() & 1) << 23
    code |= {c ^ row for c in code}
assert len(code) == 4096
assert min(c.bit_count() for c in code if c) == 8

# 23 encodes infinity. Fractional linear map s -> -1/s exchanges 0/infinity.
perm = [23] + [-pow(i, -1, 23) % 23 for i in range(1, 23)] + [0]
assert all(perm[perm[i]] == i != perm[i] for i in range(24))
def transform(c):
    return sum(((c >> i) & 1) << perm[i] for i in range(24))

pairs = [(i, j) for i, j in enumerate(perm) if i < j]
assert len(pairs) == 12
points, fixed = set(), set()
hist = Counter()
for c in code:
    d = transform(c)
    assert d in code
    dist = (c ^ d).bit_count()
    hist[dist] += 1
    values = sum(((c >> i) & 1) << k for k, (i,j) in enumerate(pairs))
    if dist == 0:
        fixed.add(values)
    elif dist == 8:
        mask = sum(int(((c >> i) ^ (c >> j)) & 1 == 0) << k for k, (i,j) in enumerate(pairs))
        points.add((mask, values & mask))

assert hist == {0:64, 8:960, 12:2048, 16:960, 24:64}
assert len(points) == 480
assert len(fixed) == 64
assert Counter(Counter(m for m, b in points).values()) == {32:15}
for mask in (2367, 1991, 3832):
    signs = {mask & c for c in fixed}
    orbit = {(mask, b) for b in signs}
    assert len(orbit) == 32 and not (points & orbit)
    points |= orbit
assert len(points) == 576
assert all(m.bit_count() == 8 and m & a == a for m,a in points)
for (m,a),(n,b) in combinations(points,2):
    q = m & n
    assert q.bit_count() - 2*((a ^ b) & q).bit_count() <= 4
print("576 local signed octads; all 165600 pairs verified")
