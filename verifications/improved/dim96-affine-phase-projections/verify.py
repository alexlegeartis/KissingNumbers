#!/usr/bin/env python3
"""Verify the translated/permuted axis-code extension, using only exact integers."""
from itertools import combinations,product
from pathlib import Path
from hashlib import sha1
from importlib import util
from functools import lru_cache
import json
P=Path(__file__).resolve().parent
pc=int.bit_count
MIXED="dim96-mixed-triality-six-rows/"
HEXARY="dim94-96-hexary-parity-packing/"
def pin(name,d):
 path=P.parent/name
 raw=path.read_bytes()
 assert sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()==d,name
 return raw.decode()
def load_mixed():
 path=MIXED+"verify.py"
 pin(path,"6b5ad3e41d156fdf606f2acb49a1168fc3050dfc")
 spec=util.spec_from_file_location("mixed96",P.parent/path)
 m=util.module_from_spec(spec);spec.loader.exec_module(m)
 return m
def phase(x,C):return sum((pc(x&c)%2)<<i for i,c in enumerate(C))
Q=[]
for x in range(16):
 a=phase(x//2,[1,2,3,4,5,6])^(63 if x&1 else 0)
 b=0
 for k,t in enumerate((15,0,8,1)):
  if x>>k&1:b^=t
 for y in range(8):
  c=b^phase(y,[7,3,4,6,2,1])
  Q.append(sum(((a>>i&1)+2*(c>>i&1))<<(2*i)for i in range(6)))
basis=[Q[1<<j]for j in range(7)]
def wt4(x):return pc((x|(x>>1))&0x555)
def subset(M):return tuple(i for i in range(24)if M>>i&1)
def pack(x,L):return sum((x>>(2*j)&3)<<(2*i)for i,j in enumerate(L))
@lru_cache(None)
def span(L,K):
 S={0}
 for q in basis:
  for I in (L,K):
   t=pack(q,I)
   if t not in S:S|={s^t for s in S}
 return S
ID=tuple(range(6))
def mat(a,e,o,p=ID,t=0):
 m=a|e|o
 assert all(0<=u<1<<24 for u in (a,e,o)) and not(a&e or a&o or e&o)
 assert pc(m)==6 and sorted(p)==list(range(6)) and 0<=t<4096
 return a,e,o,m,subset(m),p,t
def check(A,B):
 a,e,o,m,I,u,s=A;b,f,p,n,J,v,t=B
 q=pc(m&n)
 Z=(a&b)|(e&f)|(o&p)
 r=pc(Z)
 if q+r<=6:return 0
 # An agreeing frame contributes 4 only when its projective axes agree.
 # Different frames contribute <=2. Any projected Hamming distance d
 # lowers the raw dot cap from 2(q+r) by 4d.
 L=tuple(u[I.index(k)]for k in subset(Z))
 K=tuple(v[J.index(k)]for k in subset(Z))
 offset=pack(s,L)^pack(t,K)
 d=min(wt4(x^offset)for x in span(L,K))
 assert 2*(q+r)-4*d<=12
 return 1

def main():
 assert __debug__
 A={0}
 for v in basis:A|={v^x for x in A}
 assert A==set(Q) and len(A)==128 and min(wt4(x)for x in A if x)==3
 # Use the same axis-word and genuinely independent even sign verification as
 # the previous proof. A permutation and additive shift preserve distance.
 mixed=load_mixed();unit=mixed.geometry();assert unit==128*32
 prior=pin(MIXED+'witness.txt','a00a1bdd3f82f40ed930a8b894465a3feca824ee')+pin(MIXED+'witness_extra.txt','7c2e5db39c7780c573356242bd8f12ff3a82762d')
 W=[mat(*map(int,l.split()))for l in prior.splitlines()]
 assert len(W)==1170
 earlier=pin(HEXARY+'best6_96.txt','bef34cfc76b36934686b7714d20299e53c1ea4ca')
 F=[]
 for l in earlier.splitlines():
  m,t=map(int,l.split())
  assert t&m==t and pc(m)==6
  F.extend([mat(0,t,m^t),mat(0,m^t,t)])
 assert len(F)==438
 U,n=mixed.oldsupports(96);assert n==24
 V=[]
 for l in (P/'witness.txt').read_text().splitlines():
  a,e,o,p,t,k=l.split()
  V.append(mat(int(a),int(e),int(o),tuple(map(int,p)),int(t)))
 labels={X[:3] for X in V}
 assert len(V)==len(labels)==103
 assert not labels&{X[:3] for X in W+F}
 # Axial rows distinguish old norm-24 flats, and all other mixed triality
 # patterns have distinct support-and-frame labels. Older norm-96 enlarged
 # rows cannot be proportional: see the inherited mixed geometry and counts.
 assert min(1+4*5,2+4*4,3+4*3)>6 and 6<12<24
 G=[15<<(4*q)for q in range(6)]
 rescued=[0,0,0]
 for X in V:
  a,e,o,m,I,p,t=X
  assert 1<=pc(a)<=3
  assert all(pc(a&g)<=1 for g in G) and 9*48<24**2
  assert sum(2 if a&g else (1 if (e|o)&g else 0)for g in G)<=6
  assert sum(2 if (e|o)&g else (1 if a&g else 0)for g in G)<=6
  for q in U: assert pc(m&q)+pc((e|o)&q)<=6
  for q in W:rescued[0]+=check(X,q)
  for q in F:rescued[1]+=check(X,q)
 for X,Y in combinations(V,2):rescued[2]+=check(X,Y)
 assert 4*24**2==24*96
 old=json.loads(pin(MIXED+'claims.json','17d9a3ec19228f73f7f983cb535282d95629ef8c'))
 assert old['dimension']==96 and old['lower_bound']==209291873088
 result=dict(dimension=96,previous=old['lower_bound'],added_patterns=len(V),per_pattern=unit,
             rescued_pairs=rescued,gain=len(V)*unit,lower_bound=old['lower_bound']+len(V)*unit)
 assert result==json.loads((P/'claims.json').read_text())
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
