#!/usr/bin/env python3
"""Exact 62D joint support extension with orthogonal two-mark heads."""
from itertools import combinations
from collections import Counter,defaultdict
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from hashlib import sha256
import base64,gzip,struct,json,runpy,sys
assert __debug__
here=Path(__file__).resolve().parent
parent=Path(sys.argv[1]) if len(sys.argv)>1 else here.parent
if len(sys.argv)==1 and not (parent/'dim62-63-guarded-rephase').is_dir():
    parent=here.parent/'prior'
p=parent/'dim62-63-guarded-rephase'
pins={'verify.py':'29957c51f761edff380d832dece779fd9b26a44581a85103cb873f61e92c81b9',
'witness62.bin.gz.b64':'321e121bbc6c3834f35d63e79399b3302a197898988cfa75edb0724001beb7a7',
'witness63.bin.gz.b64':'78abe2f27b448e7c973e4a1b4bd5e798151bb704c48f0e9a407835524a91081e'}
for name,v in pins.items():assert sha256((p/name).read_bytes()).hexdigest()==v
with redirect_stdout(StringIO()):
    D=runpy.run_path(str(p/'verify.py'))
old=D['sparse']; rows=D['dense']
assert D['res']['62']==170279948
M=(0,9,20,41); mask=sum(1<<j for j in M)
SQ=((0,1),(1,0),(0,-1),(-1,0))
hexagon=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
lab=(0,1,3,7,6,4)
def fmt(s):
    t=[i for i,k in enumerate(M) if s>>k&1]
    assert s.bit_count()==16 and len(t)<=2
    a,b=(sum(SQ[i][j] for i in t) for j in range(2))
    assert a*a+b*b==len(t)
    # Exact truth value A+B sqrt(3)>len(t).
    def unsafe(A,B,h):
        if B==0:return A>h
        if B>0:return A>=h or 3*B*B>(h-A)**2
        return A>h and (A-h)**2>3*B*B
    risk=[l for l,(x,y) in zip(lab,hexagon) if unsafe(a*x,b*y,len(t))]
    assert not risk if not t else len(risk)>0
    v={sum((l>>i)&1 for i in t)%2 for l in risk}
    assert len(v)<=1
    # If full ordinary signs align on an orthogonal support, their parity
    # is the sum of the removed dense bits. Choose the other sign class.
    parity=(next(iter(v))^1) if v else 0
    assert all(not unsafe(a*x,b*y,len(t)+2) for x,y in hexagon)
    return s&~mask,a,b,parity,len(t)
P=sorted(s for s in old if (s&mask).bit_count()<=1)
syn=defaultdict(list)
G=set()
B=gzip.decompress(base64.b64decode((p/'witness62.bin.gz.b64').read_text()))
# Replay exact incumbent local checks and their assignment. Its checker
# above already independently enumerates and verifies incumbent angles.
j=0
while j<len(B):
    s,n=struct.unpack_from('<QB',B,j);j+=9
    H=fmt(s)[0]
    pos=[1<<k for k in range(64) if H>>k&1]
    for _ in range(n):
        i,t=struct.unpack_from('<HH',B,j);j+=4
        f=sum(pos[k] for k in range(len(pos)) if t>>k&1)
        syn[s].append((f,0));syn[P[i]].append((f,1))
    G.add(s)
P=set(P)|G
variant=here/'merged62.txt.gz.b64'
blocks=gzip.decompress(base64.b64decode(variant.read_text())).splitlines()
n,e=map(int,blocks[0].split())
assert len(blocks)==1+n+e
new=[int(s) for s in blocks[1:1+n]]; add=set(new)
assert len(add)==len(new) and not P&add
stock=P|add
data={s:fmt(s) for s in stock}
for s in stock:
    if s&mask:
        assert all(not ((s&r).bit_count()&1) for r in rows)
    syn[s].append((s,data[s][3]))
E={}
for line in blocks[1+n:]:
    a,b,f,t=map(int,line.split())
    assert a in add and b in stock and a!=b and t in (0,1)
    assert f.bit_count()==8 and f&(a&b)==f
    key=(min(a,b),max(a,b));assert key not in E;E[key]=f
    syn[a].append((f,t^1));syn[b].append((f,t))
# Verify extendability and all affine assignments on full 16 columns.
for s,L in syn.items():
    R={0:0}
    for g,t in L:
        if g in R:assert R[g]==t
        else:R.update({a^g:x^t for a,x in list(R.items())})
    assert len(R) in (2,4,8,16,32) and all(g in (0,s) or g.bit_count()==8 for g in R)
    # Minimum dual weight 8 guarantees full rank on 0–2 marked dummy
    # bits; fixing them all to zero keeps the prescribed affine syndrome.
def cnt(a,b):
    x,h,k,_,_=data[a];y,p,q,_,_=data[b]
    z=(x&y).bit_count()+h*p+k*q
    if z<=8:return 0
    assert z<=10
    f=E[(min(a,b),max(a,b))]
    assert (f&~mask)&(x&y)==(f&~mask)
    return 1
total=sum(cnt(x,y) for x in new for y in P)+sum(cnt(x,y) for x,y in combinations(new,2))
assert total==len(E)
sizes=Counter(data[s][4] for s in new)
# Planes: |s_i|<=1 in every physical coordinate therefore cross dot<=4;
# poles safe by the same fact. A sign coset has internal distance >=4.
assert all(max(map(abs,data[s][1:3]))<=1 for s in stock)
gain=sum((2048>>k)*v for k,v in sizes.items())
assert gain==1854464
claims={'status':'PASS','lower_bounds':{'62':170279948+gain}}
assert claims==json.loads((here/'claims.json').read_text())
print(json.dumps(dict(status='PASS',old=170279948,additional_banks=sizes,additional_edges=e,added=gain,lower_bound=170279948+gain)))
