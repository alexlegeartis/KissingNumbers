"""Exact 82D–84D Γ72–E8 stock plus a compatible D4 sidecar."""
from pathlib import Path
from hashlib import sha256
from collections import Counter
from itertools import combinations,permutations
from functools import lru_cache
import json,runpy,sys
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
up=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else here.parent/'dim80-gamma72-partial-fibres'
base=up.parent/'dim77-79-gamma72-seven-charts'
pins={
  "verify.py": "c90daad443667a237429129bbf1a7e4cec67f7a224498409a1c47066240e5e26",
  "witness.json": "76271222491054d5596ed16c245371a7ac72053bb9370513ebe2d7c0b982b3d1"
}
assert all(sha256((up/n).read_bytes()).hexdigest()==h for n,h in pins.items())
sys.argv=[str(up/'verify.py'),str(base),str(base.parent/'dim73-82-gamma72-depth6-lattice')]
R=runpy.run_path(sys.argv[0]);P=R['P0'];M=P['M4'];p0=R['p'];prior=R['bound']
assert prior==6228431712 and len(p0)==39
T=np.array(P['T'],dtype='i8');pk=P['p24']
def rank(xs):
 p={}
 for x in xs:
  while x:
   i=x.bit_length()-1
   if i not in p:p[i]=x;break
   x^=p[i]
 return p
def arr(x):
 return np.array([(x>>i)&1 for i in range(72)],dtype='i8')
def syndrome_trip(x):
 y=(arr(x)@T%2).reshape(3,24)
 A=(y@P['A']%2)@pk; B=(y@P['B']%2)@pk
 return list(zip(map(int,A),map(int,B)))
def mi(a,b):
 if a==b==0:return 0
 if b in M.get(a,()):return 4
 if not P['orth'](a,b):return 6
 return 8

Q0={}
W0=P['W']; G=P['F']; AA=P['A']; BB=P['B']; pk=P['p24']
for a,b,w in zip((W0@AA%2)@pk,(W0@BB%2)@pk,W0):Q0[int(a),int(b)]=w
def span(B):
 V=[0]
 for x in B:V+=[r^x for r in V]
 return V
UA={a:i for i,a in enumerate(span(rank(M).values()))}
UB={a:i for i,a in enumerate(span(rank(sum((int(x)&1)<<i for i,x in enumerate(r))for r in P['B']).values()))}
assert len(UA)==len(UB)==4096
def enc(a,b,t):
 return UA[a]+4096*UA[b]+4096**2*UB[t]
@lru_cache(None)
def rawshort(a,b):
 assert mi(a,b)==8
 for u,S in M.items():
  if not u^a:continue
  xy=next((q for q in S if (u^a,q^b)in Q0),None)
  if xy is not None:
   x=Q0[u,xy];y=Q0[u^a,xy^b];break
 else:raise AssertionError((a,b))
 assert x@G@y==0
 v=x+y; dd=W0@G@v
 N8=np.unique(np.vstack([v,-v,v-2*W0[dd==4]]),axis=0)
 N12=np.unique(np.vstack([v-2*W0[dd==3],-v+2*W0[dd==3]]),axis=0)
 assert len(N8)==48 and len(N12)==4096
 for u,n in [(N8,8),(N12,12)]:
  assert np.all(((u@G)*u).sum(1)==n)
  assert np.all((u@AA%2)@pk==a) and np.all((u@BB%2)@pk==b)
 return N8,N12
pers={}; in0=[]
for p in p0:
 raw=(arr(p)@T).reshape(3,24)
 s=tuple(int(x)for x in (raw%2).ravel())
 if not any(s):in0.append(p)
 else:pers.setdefault(s,[]).append(p)
assert len(pers)==3 and len(in0)==6
print('depth4 cross bad parity groups',list(map(len,pers.values())),flush=True)
def local(r,y,v):
 Z=r-y
 assert np.all(Z%2==0)
 Z=Z//2
 u=(Z@AA%2)@pk;t=(Z@BB%2)@pk
 x=np.array([UA[int(a)]for a in u],dtype='u8')
 b=np.array([UB[int(a)]for a in t],dtype='u8')
 return np.unique((x*(1 if v==0 else 4096 if v==1 else 4097)) ^ (b<<24))
def offset(q):
 z=arr(q)@T
 assert np.all(z%2==0)
 return P['syndrome'](z//2)
# For each raw parity pattern, enumerate the exact four possible short profiles.
# A sidecar phase is in the original integral-block chart, so its raw parity is zero.
badsets=[]
for ps in pers.values():
 p=ps[0]; y=(arr(p)@T).reshape(3,24)
 sg=syndrome_trip(p)
 assert all(mi(*w)==8 for w in sg)
 Z=[rawshort(*x)for x in sg]
 LB=[[local(x,y[j],j)for x in p]for j,p in enumerate(Z)]
 shifts=[enc(*offset(q^p))for q in ps]
 checks=[]
 for ip in [-1,0,1,2]:
  i=ip if ip>=0 else 2
  j,k=[u for u in range(3)if u!=i]
  pair=np.unique((LB[j][0][:,None]^LB[k][0][None,:]).ravel())
  checks.append((pair,LB[i][int(ip>=0)]))
 badsets.append((checks,shifts))
def cross(q):
 c=enc(*q)
 for checks,shifts in badsets:
  for shift in shifts:
   for pair,z in checks:
    ask=z^np.uint64(c^shift); i=np.minimum(np.searchsorted(pair,ask),len(pair)-1)
    if np.any(pair[i]==ask):return False
 return True

def xor(a,b):return tuple(u^v for u,v in zip(a,b))
def pass4(p):
 a,b,t=p
 if not(a or b):return bool(t)
 if a and b and a!=b:return True
 u=a or b
 return t not in M[u] and P['orth'](u,t)
def raw(a,b,t):
 z=[(w>>i)&1 for w in (a^t,b,0)for i in range(24)]
 assert P['syndrome'](z)==(a,b,t)
 return z
def pass6(p):
 a,b,t=p
 if not(a and b and a!=b):return False
 V=(a,b,a^b)
 if any(x^t in M[y]for s,y in combinations(V,2)for x in M[s]):return False
 return all(P['orth'](y,x^t)for x0,y in permutations(V,2)for x in M[x0])
W=json.loads((here/'witness.json').read_text())
D=[(tuple(q),tuple(z))for q,z in W['sidecar']]
E=set()
for i,j in combinations(range(4),2):
 for a in (-1,1):
  rr=[0]*4;rr[i]=1;rr[j]=a;E.add(tuple(rr))
assert len(D)==12 and set(z for q,z in D)==E
def dot(a,b):return sum(x*y for x,y in zip(a,b))
assert all(dot(z,z)==2 for q,z in D)
old0=[offset(r)for r in in0]
for q,z in D:
 # The pinned count theorem applies to each phase.
 P['prove'](raw(*q),6)
 assert P['contact'](raw(*q),6)==136704
 assert cross(q) and all(pass4(xor(q,p))for p in old0)
for (p,z),(q,w) in combinations(D,2):
 h=abs(dot(z,w));assert h in (0,1)
 assert pass6(xor(p,q)) if h else pass4(xor(p,q))
shapes={81:[(1,-1,0,0)],82:sorted(z for z in E if z[3]==0 and sum(z)==0),83:sorted(z for z in E if z[3]==0),84:sorted(E)}
assert [len(x)for x in shapes.values()]==[1,3,6,12]
expect={}
for d,n0 in [(81,6228705122),(82,6229187942),(83,6229748820),(84,6231393584)]:
 rays=len(shapes[d])
 n=prior+2*rays*(136704+1)
 expect[str(d)]={'prior':n0,'lower_bound':n,'gain':n-n0,'sidecar_rays':rays}
 print('PASS',d,n,'gain',n-n0)
assert expect==W['results']

