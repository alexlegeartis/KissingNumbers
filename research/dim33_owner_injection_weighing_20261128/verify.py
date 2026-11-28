#!/usr/bin/env python3
# Usage: python verify.py witness.npz lower132_33d_clean_le3_exact_owner_assignment_20261128.json
# This proves only the clean <=3-original-owner full-bank nonaugmentation theorem.
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json, sys
import numpy as np

payload = Path(sys.argv[1]).read_bytes()
assert sha256(payload).hexdigest() == "bda549debc08b487691747376cb8a4309c9f382bced78cb0b47b29979e0c0960"
with np.load(sys.argv[1], allow_pickle=False) as z:
    original = tuple(map(int, z["supports"]))
    rows = z["sparse"].astype(np.int64)
assert len(original) == len(set(original)) == 1837 and rows.shape == (4348,33)
assert all(m < 2**33 and m.bit_count() == 8 for m in original)
assert (np.sum(rows**2, axis=1)==32).all()
C = [[comb(n,k) for k in range(6)] for n in range(33)]
mask_bits = lambda m: tuple(i for i in range(33) if (m >> i) & 1)
owners = [-1] * comb(33,5)
for j,m in enumerate(original):
    for T in combinations(mask_bits(m), 5):
        q = sum(C[a][i+1] for i,a in enumerate(T))
        assert owners[q] < 0
        owners[q] = j
guards = set()
for row in rows:
    s = sum(1<<i for i in range(33) if row[i])
    h = sum(1<<i for i in range(33) if abs(row[i])==4)
    if s.bit_count() > 2:
        assert all(abs(int(x)) in (0,2,4) for x in row)
        guards.add((s,h))
assert len(guards)==44
inc = {}
old = set(original)
hist = [0]*4
def gen(a, start, left, m, own, next4):
    if not left:
        if m in old or any((m&s).bit_count() + (m&h).bit_count()>4 for s,h in guards):
            return
        inc[m] = own
        hist[len(own)] += 1
        return
    for j in range(start, 34-left):
        U = set(own)
        for x in next4:
            q = owners[x+C[j][5]]
            if q >= 0:
                U.add(q)
                if len(U)>3: break
        if len(U)>3: continue
        extra = [C[x][1]+C[y][2]+C[z][3]+C[j][4] for x,y,z in combinations(a,3)]
        gen(a+(j,),j+1,left-1,m|(1<<j),frozenset(U),next4+extra)
gen((),0,8,0,frozenset(),[])
assert hist==[0,109,659,1601]
G = {i:set() for i in range(len(original))}
for m,V in inc.items():
    for i in V: G[i].add(m)
def compatible(x,y):
    return (x&y).bit_count()<=4
Q = {m:[{n for n in G[i] if compatible(m,n)} for i in V] for m,V in inc.items()}
live=set(inc)
def good(m):
    S=sorted((q & live for q in Q[m]),key=len)
    if any(not q for q in S): return False
    def walk(i,A):
        if i==len(S): return True
        for w in S[i]:
            if w in A:
                if walk(i+1,A): return True
            elif all(compatible(w,z) for z in A) and walk(i+1,A+(w,)):
                return True
        return False
    return walk(0,())
peel=[]
while True:
    kill={m for m in live if not good(m)}
    if not kill: break
    peel.append(len(kill))
    live-=kill
assert peel==[948,156,66,7,4,1] and len(live)==1187
data=json.loads(Path(sys.argv[2]).read_text())
pairs=data["assignment"]
assert len(pairs)==len(live)
assignment={int(m):int(j) for m,j in pairs}
assert set(assignment)==live
bin={}
for m,j in assignment.items():
    assert j in inc[m]
    assert all(not compatible(m,s) for s in bin.get(j,[]))
    bin.setdefault(j,[]).append(m)
print("verified: 2369 clean supports, 1187 certified by owner injection")
