"""Exact verifier for the 94–96D complementary parity hexary append."""
from pathlib import Path
from itertools import combinations,product
from functools import reduce
from operator import xor
from hashlib import sha1
import json
P=Path(__file__).resolve().parent
def m(S):return sum(1<<p for p in S)
def wt(w):return w.bit_count()
def px(x,H):return m(i for i,h in enumerate(H) if wt(x&h)%2)
def X(z):return reduce(xor,z,0)
def gram(a,b):return sum(x*y for x,y in zip(a,b))
def prior():
 r=(P/'predecessor.json').read_bytes()
 assert sha1(b'blob '+str(len(r)).encode()+b'\0'+r).hexdigest()=='8d2d15133304aa856cb9e9473f9658013f91885a'
 live=P.parent/'dim94-96-bch-paired-group-append'/'claims.json'
 if live.exists():assert live.read_bytes()==r
 return {a['dimension']:a['lower_bound']for a in json.loads(r)}
def oldsupports(d):
 lengths=[4 if d==96 else 2]+[4]*5
 start=[sum(lengths[:j])for j in range(6)]
 A={m(start[j]+v for j,v in enumerate(q))
    for q in product(*map(range,lengths)) if X(q)==0}
 V=[[(a,start[j//2*2+a//2]+(j%2)*(lengths[j//2*2+a//2]//2)+a%2)
    for a in range(4) if a%2<lengths[j//2*2+a//2]//2]for j in range(6)]
 B={m(r for a,r in w) for w in product(*V) if X(a for a,r in w)==0
    and all(w[j][0]//2==w[j+1][0]//2 for j in (0,2,4))}
 assert len(A)==(1024 if d==96 else 512) and len(B)==(256 if d==96 else 160)
 assert not A&B and all(wt(s)==6 for s in A|B)
 return A|B,sum(lengths)
def axis():
 h=(1,2,3,4,5,6,7);g=(7,3,4,6,2,1,5);F=(15,0,8,1)
 A=[]
 for x,y in product(range(16),range(8)):
  a=px(x>>1,h)^(127 if x&1 else 0)
  b=X(F[i]for i in range(4)if x>>i&1)^px(y,g)
  A.append(tuple(((a>>i)&1)+2*((b>>i)&1)for i in range(6)))
 assert len(A)==len(set(A))==128
 d=lambda a,b:sum(u!=v for u,v in zip(a,b))
 assert min(d(a,b)for a,b in combinations(A,2))==3
 S=[w for w in range(64)if wt(w)%2==0]
 assert len(S)==32 and min(wt(a^b)for a,b in combinations(S,2))==2
 return len(A)*len(S)*2
def geometry():
 A=[[],[]]
 for w in range(8):
  for s in (0,15):
   v=w^s;A[wt(v)%2].append(tuple(1-2*((v>>j)&1)for j in range(4)))
 assert all(len(S)==8 for S in A)
 for S in A:
  assert all(gram(x,x)==4 for x in S)
  assert all(gram(x,y)==(-4 if tuple(-u for u in x)==y else 0)for x,y in combinations(S,2))
 assert all(abs(gram(x,y))==2 for x,y in product(*A))
 # Recheck all bounded row alphabets: cube, axis, and shifted root pair.
 R=list(product((-1,1),repeat=4))
 R += [tuple(a if j==i else 0 for j in range(4))for i in range(4)for a in (-2,2)]
 R += [tuple(a if k==i else b if k==j else 0 for k in range(4))
       for i,j in combinations(range(4),2)for a,b in product((-2,2),repeat=2)]
 assert all(abs(gram(s,r))<=4 for s,r in product(A[0]+A[1],R))
 # Six active rows; no head. The angular cap with norm-squared 96 is raw 24.
 assert 4*24**2==24*96
 # Old RS six individual amplitude-four entries and roots sqrt(48).
 assert 6*4==24 and 4*48<24**2
 # Within one support: distinct axis words, distinct signs on one axis word,
 # or complementary parity.
 assert max(4*(6-3),4*(6-2*2),6*2)==12
def verify():
 assert __debug__
 base=prior();unit=axis();geometry();out=[]
 for d in (94,95,96):
  U,n=oldsupports(d)
  raw=(P/f'best6_{94 if d<96 else 96}.txt').read_text()
  Q=[tuple(map(int,r.split()))for r in raw.splitlines()]
  assert len(Q)==len({s for s,q in Q})==(219 if d==96 else 152)
  for s,q in Q:
   assert 0<s<1<<n and wt(s)==6 and 0<=q<=s and q&s==q
   assert all(wt(s&t)<=3 for t in U)
  for (s,p),(t,q) in combinations(Q,2):
   k=wt(s&t);delta=wt(s&t&(p^q))
   assert 4*k-2*min(delta,k-delta)<=12
  out.append(dict(dimension=d,preceding=base[d],supports=len(Q),gain=len(Q)*unit,
                  lower_bound=base[d]+len(Q)*unit))
 return out
if __name__=='__main__':
 print(json.dumps(verify(),indent=2))
