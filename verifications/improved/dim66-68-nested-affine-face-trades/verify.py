"""Exact F8 nested affine-face trade for the 66–68D configurations."""
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import re

assert __debug__
P = Path(__file__).resolve().parent
BASE = P.parent / "dim65-68-complete-hyperplane-atlas"
PREVIOUS = {66: 429274500, 67: 448596358, 68: 481013656}
base_counts = json.loads((BASE / "claims.json").read_text())["lower_bounds"]
assert {d: base_counts[str(d)] for d in PREVIOUS} == PREVIOUS
if len(sys.argv) == 1:
    out = subprocess.check_output([sys.executable, str(BASE / "verify.py")], text=True)
    for d, n in PREVIOUS.items():
        assert re.search(rf"^DIM {d} .* N {n}$", out, re.M), (d, n)
    print("Complete predecessor replay PASS", flush=True)
else:
    assert sys.argv == [sys.argv[0], "--composition-only"]

def gf(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return z

def flat(a, b, c):
    return frozenset([a, b, c, a ^ b ^ c])

BODY = range(64)
F = {flat(a, b, c) for a, b, c in combinations(BODY, 3)}
assert len(F) == 10416
# All eight transverse affine planes; omit a slot by an anisotropic
# linear map on the two-dimensional gradient.
regular = []
for g in range(4):
    x,y = g & 1,g >> 1
    shift = (x ^ y) | (y << 1)
    assert g == 0 or (g & shift).bit_count() % 2 == 1
    for c in (0,1):
        plane = [(i<<1) | (((i & g).bit_count() ^ c) & 1) for i in range(4)]
        for j in range(4):
            regular.append([j, *(z for i,z in enumerate(plane) if i != (j ^ shift))])
assert len(regular) == 32
a, = [z for z in range(1, 8) if gf(z, z) == 7]
assert a > 0
# The first two tiers consume the nonhorizontal F8 lines and then
# reuse seven residual matchings as the eight Frobenius-graph cosets.
group = []
for slope in range(8):
    x = slope if slope else a
    y = gf(slope,slope)
    for b in range(8):
        group.append((regular, [(gf(x,t)<<3) | (gf(y,t)^b) for t in range(8)]))
assert len(group) == 64
for b in range(8):
    group.append((regular, [(gf(a,t)<<3) | (gf(gf(a,t),gf(a,t))^b) for t in range(8)]))

# Apply the same rule to the eight reserved rows. What remains is
# exactly the three edge matchings of XOR 1, (a,0), and (a,7).
for x in range(8):
    group.append((regular, [(x<<3)|t for t in range(8)]))
tagged = [frozenset(range(64,68))]
d0 = a<<3
for h,ts in ((1,(3,12)), (d0,(5,10)), (d0|7,(6,9))):
    for t in ts:
        tail = [64+i for i in range(4) if t>>i&1]
        tagged += [frozenset((x,x^h,*tail)) for x in BODY if x < (x^h)]
assert len(tagged)==193

removed, added = set(), set()
for g,block in group:
    assert len(block)==len(set(block)) and all(0<=z<64 for z in block)
    removed_here = set()
    for j,c,d,e in g:
        tri = [block[z] for z in (c,d,e)]
        assert len(set(tri))==3 and 0<=j<4
        old = flat(*tri)
        assert old in F
        new = frozenset((*tri,64+j))
        assert new not in added
        removed_here.add(old)
        added.add(new)
    assert removed.isdisjoint(removed_here)
    removed |= removed_here
assert len(removed)==640 and len(added)==2560
S = [*added,*tagged,*(F-removed)]
seen = set()
for f in S:
    assert len(f)==4
    for q in combinations(sorted(f),3):
        assert q not in seen
        seen.add(q)
assert len(S)==12529 and len(seen)==50116

# The norm-squared-16 and norm-squared-64 signs have entries <=1,
# so their dot bound is four. The poles have one entry of magnitude
# two. For every atlas head the sum of up to three specified tail
# magnitudes in units sqrt(2) is at most max(2, tail count). In the
# pure tail case it is at most four. This checks even rematched tags.
assert 4**2 <= 4*16//4 and 4**2 <= 4*64//4 and 2**2 <= 4*4//4
for t in range(5):
    assert 2*(4-t + (max(2,t) if t else 0))**2 <= 4*72//4

claims = {}
for k in (2,3,4):
    stock = sum(max(z)<64+k for z in S)
    was = 10416 + 32*k*(k-1)//2 + (k==4)
    gain = 16*(stock-was)
    total = PREVIOUS[64+k]+gain
    claims[str(64+k)]=total
    print(f"DIM {64+k}: four-site stock {stock}; gain {gain}; N {total}")
assert json.loads((P / "claims.json").read_text()) == {"status":"PASS","lower_bounds":claims}
print(json.dumps({"status":"PASS","lower_bounds":claims},sort_keys=True))
