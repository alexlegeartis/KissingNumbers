#!/usr/bin/env python3
"""Punctured RM80 with antipodal two-scale heads; exact standalone."""
from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from hashlib import sha256
import json,runpy,sys
assert __debug__
root=Path(__file__).resolve().parent
pp=Path(sys.argv[1]) if len(sys.argv)>1 else root.parent/'dim62-67-rm2-selfsparse'
assert all(sha256((pp/f).read_bytes()).hexdigest()==sha for f,sha in {
 'claims.json':'b7cb0c0e63c59f8fe6a4d155836d8597c1a17b311465752d0b81f6665e1ccabb',
 'gray.cpp':'6d08f0a74d9cca4f0697775ad38959a3f669acd3261c11377a787da366c45f1b',
 'plotkin.cpp':'65b9127f4c823f714137e6a867cc252b51089c1320979e8c2e22d783d0262052',
 'tails.cpp':'a09a1a152fd6bf213b358e7ca6924ad726eb60e39f4d6c51ceba77a525ba9b83',
 'verify.py':'75061620b9f7a3444df46b129b6bdbfccbe06ed7949f9aab7ca2572a40cf1e18',
 'verify_independent.cpp':'d70f7077f97dbbf380e2e6df27b95aed00e8ff3715d32c55cd2a7c4dc5eb63ce',
 'witness.txt':'5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f',
 }.items())
with redirect_stdout(StringIO()):
    data=runpy.run_path(str(pp/'verify.py'))
def cosgate(dot,n,m):
    return dot <= 0 or 4*dot*dot<=n*m
def finite_syndrome(n):
    gate=[0]*(1<<n)
    for v in range(1,1<<n):
        j=(v&-v).bit_length()-1
        gate[v]=gate[v&(v-1)]^(16+j)
    count=Counter(gate)
    assert count[0]==count[16]==2**(n-5)
    assert all(t.bit_count()%2==(x>>4) for t,x in enumerate(gate))
    assert min(v.bit_count() for v in range(1,1<<n) if gate[v]==0)==4
    return count
sig16,sig15=finite_syndrome(16),finite_syndrome(15)
M0=[m for m in range(64) if m.bit_count()<=2]
M1=[m for m in range(64) if m.bit_count()==3]
def row(m):return sum(1<<x for x in range(64) if (x&m)==m)
R=[row(m) for m in M0]
v1=[row(m) for m in M1]
C=[]
for s in data['C']:
    w=0
    for i,v in enumerate(v1):
        if (s>>i)&1:w^=v
    C.append(w)
V=data['V']
P=data['P']
assert all((p&q).bit_count()%2==0 for p in R+C for q in R+C)
def run(n,marks,headmarks=()):
    assert len(marks)==len(set(marks))==4
    bits=[1<<x for x in marks]
    mask=sum(bits)
    proj=lambda w:sum(1<<j for j,i in enumerate(bits) if w&i)
    span={0}
    for w in R:
        span|={s^proj(w) for s in span}
    assert len(span)==16
    family= [(0,0),(1,0)] if n==60 else [(0,Q(1,2)),(1,2),(2,-Q(1,2)),(3,-2)]
    for (k,h),(w,t) in combinations_with_replacement(family,2):
        # Parent d>=16, loss when deleting marked sites exactly wt(k^w).
        assert cosgate(28+2*(k^w).bit_count()+h*t,60+h*h,60+t*t)
        assert cosgate(abs(h),60+h*h,1)
        assert cosgate(1,60+h*h,1) and cosgate(4,60+h*h,4)
    zero=[s for s in V if not s&mask]
    keep_axis=bits[0] if n==61 else 0
    planes=[r for r in P if not r&(mask^keep_axis)]
    assert cosgate(3+1,16,4) # marked sparse versus a plane using its tail
    for z in range(1+(n==61)):
        assert all(cosgate(4-z+z*abs(h),60+h*h,4) for _,h in family)
    # Every retained ordinary support is orthogonal to dense parent,
    # so odd (sparse) versus even (dense) supplies one disagreement.
    assert all(all(not ((s&c).bit_count()%2) for c in R+C) for s in zero)
    assert all(cosgate(14,16,60+h*h) for _,h in family)
    assert cosgate(8,16,16) and cosgate(4,16,4) and cosgate(1,16,1)
    assert cosgate(2,4,4) and cosgate(1,4,1)
    plus=0
    for k,v in zip(headmarks,(1,-1)):
        j=marks.index(k); part=[s for s in V if (s&mask)==bits[j]]
        assert all(all(not ((s&c).bit_count()%2) for c in R+C) for s in part)
        threatening=[label for label,h in family if not cosgate(15+v*h,16,60+h*h)]
        assert len(threatening)==1
        required=(threatening[0]>>j&1)^1
        for label,h in family:
            d=15+v*h-2*(label in threatening)
            assert cosgate(d,16,60+h*h)
            if label in threatening: assert (label>>j&1)!=required
        # Same marked class shares the removed coordinate, so ordinary
        # overlap <=7. Distinct classes have opposite tail directions.
        assert cosgate(7+1,16,16) and cosgate(8-1,16,16)
        assert cosgate(15+1-2*4,16,16)
        plus+=len(part)*sig15[16*required]
    assert (n==60 and not headmarks) or len(set(headmarks))==len(headmarks)==2
    dense=len(C)*2**18*len(family)
    out={'n':n,'marks':marks,'headmarks':headmarks,'dense':dense,'ordinary_supports':len(zero),
         'marked_points':plus,'planes':len(planes),'poles':2*n}
    out['total']=dense+len(zero)*sig16[16]+plus+len(planes)*16+2*n
    return out
x=run(60,(0,20,41,61))
y=run(61,(0,20,41,61),(0,20))  # Can substitute optimized marks.
assert x['total']==66657320 and y['total']==116711210
print(json.dumps([x,y],indent=2))
