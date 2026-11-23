#!/usr/bin/env python3
"""Exact incremental Kerdock prefix check, plus optional complete ancestor replay."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement,combinations
from hashlib import sha256
from contextlib import redirect_stdout
from io import StringIO
import argparse,json,runpy,sys
from dg import code
assert __debug__
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--siblings',type=Path)
args=ap.parse_args()
pins={'claims.json':'6ffb1e9d2fc0cfc73cf3324e057aa7b33c3e7ae5666adb4efe7604a30923c451',
      'verify.py':'aa1be33e9c667c135beaf4ef2c7432a8b9efc99ec2984678d819bd8c0695f214'}
H=args.siblings.resolve() if args.siblings else P.parent
target=H/'dim62-dg-three-prefixes'
full=target.is_dir()
if args.siblings:assert full
S0=target if full else P/'source'
assert all(sha256((S0/f).read_bytes()).hexdigest()==h for f,h in pins.items())
assert sha256((P/'dg.py').read_bytes()).hexdigest()=='eefd262c41ae6d93c3407cf41222b0231d40e1e2cf14573b0691c0fece8f0eba'
prior=json.loads((S0/'claims.json').read_text())['lower_bounds']['62'];assert prior==173717908
if full:
 with redirect_stdout(StringIO()) as capture:
  sys.argv=[str(target/'verify.py'),str(H)]
  inherited=runpy.run_path(str(target/'verify.py'))
 assert inherited['out']['lower_bounds']['62']==prior
else:
 assert sha256((P/'source/replay.log').read_bytes()).hexdigest()=='77fdd45052781a2eb787538d622e028160b3062bc0a005845c2e1f4e0472ab6e'
M=(0,9,20,41);mask=sum(1<<i for i in M)
G,A=code()
assert len(G)==1024 and len(A)==128
K=G[::32]
assert len(K)==32 and len(set(K))==32
assert min((r^s^t).bit_count()for r,s in combinations(K,2)for t in A)==28
if full:
 assert tuple(inherited['G'])==tuple(G) and tuple(inherited['A'])==tuple(A)
 assert len(inherited['CS'])==60 and len(set(inherited['CS']))==60
 assert inherited['PFX']=={5:(F(5),F(1,5),F(-3,10)),9:(F(29,5),F(1),F(7,5)),12:(F(29,5),F(1),-F(7,5))}
cts={t&mask:sum((z&mask)==(t&mask)for z in A)for t in A}
assert len(cts)==16 and set(cts.values())=={8}
labels=(0,1,3,7,6,4)
hexhead=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
south=((0,1),(1,0),(0,-1),(-1,0))
raw={5:(51692,3848,-7524),9:(52444,6667,11547),12:(52444,6667,-11547),13:(64829,3721,0)}
PFX={j:tuple(F(a,10000)for a in r)for j,r in raw.items()}
S=lambda a,x,y:59+a*a+x*x+y*y
def valid(d,n,m):
 assert d<=0 or 4*d*d<=n*m,(float(d),float((n*m)**.5/2))
lo,hi=F(1732,1000),F(1733,1000)
assert lo**2<3<hi**2
def upper(a,b):return a+b*(hi if b>=0 else lo)
def lower(a,b):return a+b*(lo if b>=0 else hi)
banks=[]
for ms in range(16):
 ix=[j for j in range(4)if ms>>j&1];k=len(ix)
 if k>2 or (k==2 and sum(ix)%2==0):continue
 x,y=[sum(south[j][p] for j in ix) for p in range(2)]
 assert x*x+y*y==k
 q=set()
 for l,(u,v) in zip(labels,hexhead):
  q0,q1=lower(x*u,y*v),upper(x*u,y*v)
  assert q0>k or q1<=k
  if q0>k:q.add((l&ms).bit_count()%2)
 assert len(q)==int(k>0)
 banks.append((ms,k,x,y,1-next(iter(q))if k else 1))
assert len(banks)==9
if full:assert banks==inherited['banks']
for p,(a,x,y) in PFX.items():
 assert p not in labels and a>=1
 n=S(a,x,y)
 for l,(u,v) in zip(labels,hexhead):
  d=(p^l).bit_count()
  valid(60-2*(16-d)+a-1+upper(x*u,y*v),n,64)
 for ms,k,u,v,parity in banks:
  mismatch=(p&ms).bit_count()%2!=parity
  valid(16-k+a-1-2*mismatch+x*u+y*v,n,16)
 for d in (3+a,2+a+max(abs(x),abs(y)),1+a+abs(x)+abs(y)):
  valid(d,n,4)
 valid(max(abs(x),abs(y)),n,1)
for (p,(a,x,y)),(q,(b,u,v)) in combinations_with_replacement(PFX.items(),2):
 d=(p^q).bit_count();h=x*u+y*v;n,m=S(a,x,y),S(b,u,v)
 # Different prefixes make any coincident quadratic choice differ affinely.
 # Their unpunctured affine distance is then 32, so the 24 bound is still valid.
 gap=28 if p==q==13 else 24
 valid(60-2*(gap-d)+a*b-1+h,n,m)
 valid(60-2*(16-d)+(a-1)+(b-1)+h,n,m)
gain=60*len(K)*8
assert gain==15360
out={'status':'PASS','lower_bounds':{'62':prior+gain},'gain':gain}
assert out==json.loads((P/'claims.json').read_text())
print(json.dumps(out), 'full inherited replay' if full else 'incremental replay against pinned established ancestor')
