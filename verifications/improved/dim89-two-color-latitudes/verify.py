#!/usr/bin/env python3
"""Exact replay of four-size two-color latitude extension of the 89D code."""
import json, hashlib, subprocess
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
P=Path(__file__).resolve().parent
gen=P/"z8_generator.json"
assert hashlib.sha256(gen.read_bytes()).hexdigest()=="903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73"
G=json.loads(gen.read_text())["generator"]
assert len(G)==12 and all(len(row)==24 for row in G)
assert all(G[i][i]%2 and all(v==0 for v in G[i][:i]) for i in range(12))
assert all(tuple(G[i][:2])==(0,0)for i in range(2,12))
assert len({(a*G[0][0]%8,(a*G[0][1]+b*G[1][1])%8)for a,b in product(range(8),repeat=2)})==64
# Recover alphabets from Walsh characters and their single coordinate reflection.
A=[tuple((-1)**(((q>>2)+(q&1)*a+((q>>1)&1)*b)%2)for a,b in product(range(2),repeat=2))for q in range(8)]
B=[(-x[0],*x[1:])for x in A]
C=[tuple((2 if q<4 else -2)if j==q%4 else 0 for j in range(4))for q in range(8)]
T=[A,B,C];dot=lambda a,b: sum(i*j for i,j in zip(a,b))
W=[0,2,2,2,4,2,2,2]
for S in T:
 assert len(set(S))==8
 for a,b in product(range(8),repeat=2):assert dot(S[a],S[b])==4-2*W[(a-b)%8]
for S,U in combinations(T,2):
 assert max(dot(x,y)for x,y in product(S,U))==2
profiles=[
 [(0,F(0)),(1,F(0)),(2,F(0)),(3,F(1)),(7,F(-1,4))],
 [(0,F(0)),(1,F(0)),(2,F(0)),(3,F(1,4)),(7,F(-1))],
 [(0,F(0)),(1,F(0)),(2,F(0)),(3,F(0))]
]
labels=[(c,b,h)for c,L in enumerate(profiles)for b,h in L]
n=lambda h: 88+8*h*h
assert len(labels)==len({(c,b)for c,b,h in labels})==14
# Dot bound 48 for distinct words in the unpunctured same-color inherited parent.
for c,b,h in labels:
 assert 2*(40+8*h*h)<=n(h) # same fiber, other word
 assert 4*8*h*h<=n(h) # same-height fiber vs either pure pole
for (c,b,h),(d,e,k)in combinations(labels,2):
 body=44 if c!=d else 48-dot(T[c][0],T[c][0])-dot(T[c][b],T[c][e])
 ip=body+8*h*k
 assert ip>=0 and 4*ip**2<=n(h)*n(k)
# Against any of the inherited sparse weight-11/10 vectors, each occupied
# body row contributes at most 4 and their head is zero, at square radii 88/80.
for w in (10,11):
 dirs=set()
 for i,j in combinations(range(4),2):
  for x,y in product((-2,2),repeat=2):
   v=[0]*4;v[i]=x;v[j]=y;dirs.add(tuple(v))
 assert len(dirs)==24 and all(dot(v,v)==8 for v in dirs)
 assert max(dot(v,t)for v,t in product(dirs,sum(T,[])))==4
 for c,b,h in labels: assert 4*(4*w)**2<=(8*w)*n(h)
# Replay inherited appended sparse layers and use their frozen claim breakdown.
p=P/'abbreviated_parent'
sp=subprocess.run(['node',str(p/'verify.js')],check=True,text=True,capture_output=True)
cl=next(x for x in json.loads((p/'claims.json').read_text()) if x['dimension']==89)
assert cl['incumbent']==13*8**10+2 and cl['lower_bound']==14051596290
ext=629*1024*144+11*256*72
assert ext==cl['lower_bound']-cl['incumbent']==92952576
base=14*8**10+2
total=base+ext
assert total==15125338114
print(json.dumps({'dimension':89,'new_dense_fibers':14,'dense_with_poles':base,'retained_sparse':ext,'previous':cl['lower_bound'],'lower_bound':total,'gain':total-cl['lower_bound'],'exact_cross_fiber_pairs':91,'proof':'independent geometry and sparse replay; inherits pinned homogeneous minimum 24'},indent=2))
