#!/usr/bin/env python3
# Usage: python verify.py own-witness.npz onepencil.json multipencil.json
# Exact 33D dense-only replays. No full-code improvement is claimed.
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json, sys
import numpy as np

base = Path(sys.argv[1])
assert sha256(base.read_bytes()).hexdigest() == "bda549debc08b487691747376cb8a4309c9f382bced78cb0b47b29979e0c0960"
with np.load(base, allow_pickle=False) as npz:
    rows = list(map(int,npz["dense_generator"]))
    assert int(npz["dense_shift"]) == 0
def words(basis):
    H = [0]
    for w in basis:
        H += [e^w for e in H]
    return H
C = set(words(rows))
assert len(C) == 131072 and min(x.bit_count() for x in C if x) == 8

# a=19/25 and eta=sqrt(8809/1250); comparisons use exact rational squares.
assert 31*722 + 2*8809 == 32*1250
assert 62*722 - 32*1250 < 16*1250
# eta < 16 - 17a; eta + 7a < 8; eta + a < 4
assert 400-17*19 > 0 and 8809 < 2*(400-17*19)**2
assert 200-7*19 > 0 and 8809 < 2*(200-7*19)**2
assert 100-19 > 0 and 8809 < 2*(100-19)**2
def above(c, v, r): # v*sqrt(c) > r, for integers
    if v == 0: return r < 0
    if v > 0: return r <= 0 or c*v*v > r*r
    return r < 0 and c*v*v < r*r
def dot(v,w):
    return sum(x*y for x,y in zip(v,w))

p = json.loads(Path(sys.argv[2]).read_text())
assert p["pair"] == [0,18] and len(set(p["physicalDenseAxes"])) == 32
assert set(p["physicalDenseAxes"]) == set(range(1,33))
u0=p["bodyBit"]
assert p["physicalDenseAxes"][u0] == 18
K = words(p["basis"])
assert len(K)==len(set(K))==64 and all(k in C for k in K)
assert p["translate"] in C and p["translate"] not in K
assert min((x & (0xffffffff^(1<<u0))).bit_count() for x in K if x) == 14
assert 32*1250 - 2*14*722 < 16*1250
print("one pencil exact +128 dense, not a full-code bound")

q=json.loads(Path(sys.argv[3]).read_text())
H=words(q["H0"])
assert H==q["H"] and len(H)==len(set(H))==256 and all(x in C for x in H)
assert len(q["V"])==256 and set(a>>1 for a in q["V"])<=set(range(32))
R=[]
for c,m in zip(H,q["V"]):
    u, l = (m>>1)+1, m&1
    s = [0 if j in (0,u) else (-1 if c>>(j-1)&1 else 1) for j in range(33)]
    for sign in (-1,1):
        tail=[0]*33
        tail[0]=sign
        tail[u]=sign*(-1 if l else 1)
        R.append((s,tail))

for (s,h),(t,k) in combinations(R,2):
    U,V,W=dot(s,t),dot(s,k)+dot(h,t),dot(h,k)
    assert not above(17618, 19*V, 16*1250-722*U-8809*W)
print("256 replaced parents / 512 new, 130816 exact new/new checks; old/new bound 17a+eta")
print("total dense =", 131072-256+512, "; sparse retention not claimed")
