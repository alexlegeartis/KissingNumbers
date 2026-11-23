#!/usr/bin/env python3
"""Exact additive q4 and pure q8 tails over the Γ72/B-phase construction."""
from itertools import combinations,product
from pathlib import Path
from hashlib import sha256
import json,sys,runpy
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
parent=here.parent/'dim81-85-gamma72-apn-complete-dm'
pins={
 'verify.py':'f9888cb194f2b16624bbabbccdedb802800b13c266a06a12710aef66fc190e33',
 'witness.json':'5f89656cd4916c2534f9ed4c372ae625cf95d6fe82bfca50dd91af1b911ee334',
 'claims.json':'74cb4f0944966770b5959f61a41d73fe7c0610596d114398815b694cdf522728'
}
assert all(sha256((parent/n).read_bytes()).hexdigest()==s for n,s in pins.items())
sys.argv=[str(parent/'verify.py')]
R=runpy.run_path(sys.argv[0]); P=R['P']; H=R['H']
assert len(set(H))==128 and H[0]==0
for t in H[1:]:
 raw=[(z>>i)&1 for z in (t,0,0) for i in range(24)]
 assert P['syndrome'](raw)==(0,0,t)
 P['prove'](raw,4)
 assert P['contact'](raw,4)==144

G5=[(0,1,2,5),(0,1,3,7),(0,1,6,8),(0,2,7,8),(0,5,8,9),(0,6,7,9)]
G7=[(0,1,2,4),(0,1,10,11),(0,2,10,12),(0,3,10,13),(0,1,5,12),(0,1,7,13),(0,1,8,9),(0,2,8,13),(0,2,9,11),(0,3,8,11),(0,7,11,12),(0,3,9,12),(7,10,12,13)]
def develop(q,starts):
 S={tuple(sorted((x%q+j)%q+q*(x//q) for x in v)) for j in range(q) for v in starts}
 assert len(S)==len(starts)*q
 T=[t for s in S for t in combinations(s,3)]
 assert len(T)==len(set(T))
 return S
systems=[develop(5,G5),develop(7,G7)]
prior=json.loads((parent/'claims.json').read_text())['lower_bounds']
def dot(x,y):return sum(u*v for u,v in zip(x,y))
res={}
for m,n in zip(range(9,14),[18,30,30,45,65]):
 S=sorted([s for s in systems[int(m>10)] if s[-1]<m])
 assert len(S)==n
 N=[{j for j,t in enumerate(S) if j!=i and set(s)&set(t)}for i,s in enumerate(S)]
 rem=set(range(n)); C={}
 while rem:
  i=max(rem,key=lambda k:(len({C[t]for t in N[k]if t in C}),len(N[k]),-k))
  have={C[t]for t in N[i]if t in C}
  c=next(j for j in range(n) if j not in have)
  C[i]=c;rem.remove(i)
 assert all(C[i]!=C[j]for i in range(n)for j in N[i])
 W=[]
 for i in range(m):
  for a in (-2,2):
   v=[0]*m;v[i]=a;W.append((tuple(v),1))
 for k,s in enumerate(S):
  for signs in product((-1,1),repeat=4):
   v=[0]*m
   for i,h in zip(s,signs):v[i]=h
   W.append((tuple(v),2+2*C[k]+signs.count(-1)%2))
 assert len(W)==len({u for u,t in W})==2*m+16*n
 assert all(0<t<128 and dot(u,u)==4 and sum(sorted(map(abs,u))[-2:])==2 for u,t in W)
 assert all(dot(u,v)<=(0 if a==b else 2) for (u,a),(v,b) in combinations(W,2))
 U0=[(0,1,2,3,4,5,6,7)] if m<=11 else []
 if m==12: U0=[tuple(i for i in range(m)if i//4!=j)for j in range(3)]
 if m==13:
  ex=[{0,1,2,3,4},{5,6,7,8,9},{0,5,10,11,12}]
  U0=[tuple(i for i in range(m)if i not in a)for a in ex]
 assert all(len(s)==8 for s in U0)
 assert all(len(set(x)&set(y))<=4 for x,y in combinations(U0,2))
 # In a common pure 8-set, even sign parity gives distinct distance≥2.
 pure=[x for x in range(256)if x.bit_count()%2==0]
 assert len(pure)==128 and all((a^b).bit_count()>=2 for a,b in combinations(pure,2))
 # Pure old support size2 and amplitude2 gives dot at most4; mixed tails give at most2.
 # A q4 additional tail has L1 norm≤4, so 8-set-vs-q4 dot at most4.
 assert all(sum(map(abs,u))<=4 for u,a in W)
 d=str(m+72);low=prior[d]+144*len(W)+128*len(U0)
 res[d]=low
 print('PASS',d,low,'gain',low-prior[d],'q4 tails',len(W), 'pure 8-sets',len(U0))
assert {'status':'PASS','lower_bounds':res}==json.loads((here/'claims.json').read_text())
