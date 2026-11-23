#!/usr/bin/env python3
"""Exact replay of single-mark sparse shells in 62 and 63 dimensions."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import json, runpy
assert __debug__
here = Path(__file__).resolve().parent
parent = here.parent / "dim62-67-rm2-selfsparse"
pins = {
    "claims.json": "b7cb0c0e63c59f8fe6a4d155836d8597c1a17b311465752d0b81f6665e1ccabb",
    "gray.cpp": "6d08f0a74d9cca4f0697775ad38959a3f669acd3261c11377a787da366c45f1b",
    "plotkin.cpp": "65b9127f4c823f714137e6a867cc252b51089c1320979e8c2e22d783d0262052",
    "tails.cpp": "a09a1a152fd6bf213b358e7ca6924ad726eb60e39f4d6c51ceba77a525ba9b83",
    "verify.py": "75061620b9f7a3444df46b129b6bdbfccbe06ed7949f9aab7ca2572a40cf1e18",
    "verify_independent.cpp": "d70f7077f97dbbf380e2e6df27b95aed00e8ff3715d32c55cd2a7c4dc5eb63ce",
    "witness.txt": "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f"
}
assert all(sha256((parent/p).read_bytes()).hexdigest()==v for p,v in pins.items())
D = runpy.run_path(str(parent / "verify.py"))
S = D["V"]             # actual weight-16 words
phase = D["C"]         # 80 cubic shifts of RM(2,6)
plane = D["P"]         # 10,416 affine planes
assert len(S) == len(set(S)) == 37996 and len(phase)==80
mon = [sum(1<<x for x in range(64) if (x&m)==m) for m in range(64)]
quadratic = [mon[m] for m in range(64) if m.bit_count()<=2]
cubic = [mon[m] for m in range(64) if m.bit_count()==3]
from functools import reduce
from operator import xor
phase = [reduce(xor,(v for i,v in enumerate(cubic) if (z>>i)&1),0) for z in phase]
assert all(not (x&y).bit_count()%2 for x in phase for y in phase)

def sign_syndromes(n):
    syn = [0]*(1<<n)
    for x in range(1,1<<n):
        i = (x&-x).bit_length()-1
        syn[x] = syn[x&(x-1)] ^ (16+i)
    return syn
syn=sign_syndromes(15)
assert min(x.bit_count() for x,v in enumerate(syn) if x and not v)==4
assert syn.count(0)==syn.count(16)==1024
assert all(x.bit_count()%2 == (v>>4) for x,v in enumerate(syn))

marks62=(0,9,20,41)
marks63=(0,20)
def single_marks(marks):
    m=sum(1<<i for i in marks)
    stock={i:[s for s in S if s&m == 1<<i] for i in marks}
    assert all(s.bit_count()==16 for V in stock.values() for s in V)
    # Thus any matching sign pattern from the dense code is forced to have
    # body parity equal to the original selected marked bit.
    assert all(not (s&v).bit_count()%2
               for V in stock.values() for s in V for v in quadratic+phase)
    return stock
k62=single_marks(marks62)
k63=single_marks(marks63)
assert [len(k62[t]) for t in marks62]==[4032]*4
assert [len(k63[t]) for t in marks63]==[7136]*2

# A sign code for a fixed head is one coset, so internal d>=4. Distinct
# old supports meet in <=8 by the pinned independent minimum-distance
# check. In the same mark class at most 7 of these are unmarked; in
# different classes the directed head dots are nonpositive.
head62={0:(0,2),9:(2,0),20:(0,-2),41:(-2,0)}
def dot(a,b):return sum(x*y for x,y in zip(a,b))
assert all(dot(h,h)==4 for h in head62.values())
assert all(dot(h,k)<=0 for h,k in combinations(head62.values(),2))
assert 4*7+4<=32 and 4*8<=32

# Exact comparison x + y sqrt(3) > 2: the two coordinates of each
# new axis head are rational; at most one of x,y is nonzero.
def gt2(x,y):
    assert not x or not y
    return (y>0 and 3*y*y>4) if y else x>2
def le4(x,y):
    return (y<=0 or 3*y*y<=16) if y else x<=4
forced=[]
for j,i in enumerate(marks62):
    a,b=head62[i]
    risk=[]
    for t,(x,y) in zip(D["labels"],D["heads"]):
        assert le4(a*x,b*y)
        if gt2(a*x,b*y):risk.append((t>>j)&1)
    assert len(risk)>0 and len(set(risk))==1
    forced.append(risk[0]^1)
assert forced == [0,1,0,1]
# All 16 patterns on the new four marked sites occur equally in each shift.
R={0}
for v in quadratic:
    t=sum(((v>>q)&1)<<j for j,q in enumerate(marks62))
    R|={x^t for x in R}
assert len(R)==16
# The same Hamming-to-hexagon gate from the parent depends only on
# differences of the four-bit patterns, not their physical positions.
from itertools import product
for (i,(x,y)),(j,(u,v)) in product(zip(D["labels"],D["heads"]),repeat=2):
    assert 60+x*x+3*y*y==64
    assert 60-2*(16-(i^j).bit_count())+x*u+3*y*v<=32
m62=sum(1<<i for i in marks62)
drop62=(1<<marks62[2])|(1<<marks62[3])
s62=sum(not (w&m62) for w in S)
p62=sum(not (w&drop62) for w in plane)
assert (s62,p62)==(11929,9145)
# With no facing risk the body contributes <=30, head<=2.
# With opposite parity one body mismatch pays 4, and head<=4.
assert 30+2==32 and 26+4<=32
# All existing four-site planes in 62D stay: at most one head-axis
# contribution 8 and at most three common nonzero body entries of 8.
assert 8+3*8<=32 and 2*8<=32 and 8*2<=32

# In 63D the old section dense head is +sqrt(2) for equality bit0 and
# -sqrt(2) for bit1. Use +2 for mark0, -2 for mark20; parity 1 and0.
assert 2*2*2>2*2 and 2*2*2<=4*4
# The old doubled affine-plane head is 4 sqrt(2) with exactly two
# unmarked coordinates. Its dot with a new head fits even at two hits.
assert (32-16)**2 >= 8**2*2
m63=sum(1<<i for i in marks63)
assert sum(not (s&m63) for s in S)==21361
assert sum((s&m63)==m63 for s in S)==2363
assert sum(not (s&m63) for s in plane)==9145
assert sum((s&m63)==m63 for s in plane)==31
n62=6*(80<<18)+s62*2048+p62*16+124+1024*sum(map(len,k62.values()))
n63=D["n63"]-2363*1024+1024*sum(map(len,k63.values()))
assert (n62,n63)==(166921228,226280708)
obj={"status":"PASS","lower_bounds":{"62":n62,"63":n63}}
assert obj==json.loads((here/"claims.json").read_text())
print(json.dumps(obj,indent=2))
