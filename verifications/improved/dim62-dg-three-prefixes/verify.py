#!/usr/bin/env python3
"""Exact checks for three omitted prefixes in the 62D Delsarte–Goethals mixer."""
import contextlib,json,runpy,sys
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from itertools import combinations_with_replacement
from pathlib import Path
assert __debug__
P=Path(__file__).resolve().parent
H=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
Q=H/'dim62-63-dg-spike-banks'
checks={
'verify.py':'4991b47dfb003c2e6c4cc3ede21d3c40b5eb73f5f4602be98d6c4a3b6eebf4b3',
'verify62.py':'675a8ac01c3bed1360bd9a6d604358977ba218d7741afcf576005e240bdea3a9',
'dg.py':'eefd262c41ae6d93c3407cf41222b0231d40e1e2cf14573b0691c0fece8f0eba',
'sign62.json':'1eac36c60851b4e076d30bf17f7bc355233f2ac3ed7be6a6d6f923cdbed085b6',
'sign63.json':'10178e84e469a71c18f9d3bdf195b9d7d9e9b54b9b0edc788270ebeb7e068f14',
'short62.txt.gz.b64':'d7ff4f18cc13d2e68e4573edb16cd8a4fc77f18bc189702a9441839359cdbe29',
'claims.json':'e76ef18cd44a45db056166cf4ba460b54eb54d710b32c1559e462e740dad9589'}
assert all(sha256((Q/f).read_bytes()).hexdigest()==v for f,v in checks.items())
sys.path.insert(0,str(Q))
with contextlib.redirect_stdout(StringIO()):
    E=runpy.run_path(str(Q/'verify.py'))

claimed=json.loads((Q/'claims.json').read_text())
assert claimed['lower_bounds']['62']==172734988
M=(0,9,20,41); mask=sum(1<<i for i in M)
G,A,CS=E['C'],E['A'],E['cube'][:60]
assert len(G)==1024 and len(A)==128 and len(CS)==len(set(CS))==60
assert len({a&mask for a in A})==16
counts={t:sum((a&mask)==t for a in A) for t in {a&mask for a in A}}
assert set(counts.values())=={8}
labels=(0,1,3,7,6,4)
hexhead=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
south=((0,1),(1,0),(0,-1),(-1,0))
PFX={5:(F(5),F(1,5),F(-3,10)),
     9:(F(29,5),F(1),F(7,5)),
     12:(F(29,5),F(1),-F(7,5))}
S=lambda a,x,y:59+a*a+x*x+y*y
valid=lambda d,n,m: d<=0 or 4*d*d<=n*m
lo,hi=F(1732,1000),F(1733,1000)
assert lo**2<3<hi**2
def upper(a,b):return a+b*(hi if b>=0 else lo)
def lower(a,b):return a+b*(lo if b>=0 else hi)
banks=[]
for ms in range(16):
    ix=[j for j in range(4) if ms>>j&1]
    k=len(ix)
    if k>2 or (k==2 and sum(ix)%2==0):continue
    x,y=[sum(south[j][p] for j in ix) for p in range(2)]
    assert x*x+y*y==k
    q=set()
    for l,(u,v) in zip(labels,hexhead):
        q0,q1=lower(x*u,y*v),upper(x*u,y*v)
        assert q0>k or q1<=k
        if q0>k:q.add((l&ms).bit_count()%2)
    assert len(q)==(1 if k else 0)
    parity=1-next(iter(q)) if k else 1
    banks.append((ms,k,x,y,parity))
assert len(banks)==9
n_new=0
for pref,(a,x,y) in PFX.items():
    c=sum(1<<j for i,j in enumerate(M) if pref>>i&1)
    assert pref not in labels and a>=1
    assert all(counts[(v^q^c)&mask]==8 for v in CS for q in G)
    n_new+=len(CS)*len(G)*8
    n=S(a,x,y)
    for l,(u,v) in zip(labels,hexhead):
        d=(l^pref).bit_count()
        assert valid(60-2*(16-d)+a-1+upper(x*u,y*v),n,64)
    for ms,k,u,v,parity in banks:
        mismatch=(pref&ms).bit_count()%2!=parity
        assert valid(16-k+a-1-2*mismatch+x*u+y*v,n,16)
    # Affine plane and surviving tail pole comparisons, subdivided by
    # whether 0, 1, or 2 of the old unit tail axes are on the plane.
    for d in (3+a,2+a+max(abs(x),abs(y)),1+a+abs(x)+abs(y)):
        assert valid(d,n,4)
    assert valid(max(abs(x),abs(y)),n,1)
for (p,(a,x,y)),(q,(b,u,v)) in combinations_with_replacement(PFX.items(),2):
    d=(p^q).bit_count(); h=x*u+y*v
    n,m=S(a,x,y),S(b,u,v)
    assert valid(60-2*(24-d)+(a*b-1)+h,n,m)
    assert valid(60-2*(16-d)+(a-1)+(b-1)+h,n,m)
assert n_new==1474560
gain=n_new-491520-120
out={'status':'PASS','lower_bounds':{'62':claimed['lower_bounds']['62']+gain},'gain':gain}
assert out=={'status':'PASS','lower_bounds':{'62':173717908},'gain':982920}
assert out==json.loads((P/'claims.json').read_text())
print(json.dumps(out))
