from pathlib import Path
import importlib.util,json,hashlib,itertools as it,random,numpy as np
from collections import defaultdict
D=Path(__file__).resolve().parent.parent
if not (D/'dim61-66-nonary-affine-caps').is_dir():D=D/'mid_49_69_20261121'/'verifications'/'improved'
N=D/'dim61-66-nonary-affine-caps'
def pos(x):return [i for i in range(64)if x>>i&1]
def physical(m,pts):return sum(1<<j for i,j in enumerate(pts)if m>>i&1)
def local(m,pts):return sum(1<<i for i,j in enumerate(pts)if m>>j&1)
def phase():
 g=[0,1,3,2];ts=['01112','03110','21310','21132']
 A=sorted({sum(g[(z*int(t[(i+j)%5]))%4]<<(2*i)for i in range(5))>>1 for t in ts for z in (-1,1)for j in range(5)})
 assert len(A)==40 and min((x^y).bit_count()for x,y in it.combinations(A,2))==3
 return A
def setup():
 def pin(p,h):
  assert hashlib.sha256(p.read_bytes()).hexdigest()==h
  return p
 x=pin(N/'certificate.json','d2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de')
 y=pin(D/'dim63-65-affine-reconstruction/section63.json','25cae31ba6b5c8746a7d2f203746e6f0b3d4f31f391019cbf80044b5fb4d4b15')
 Z=pin(D/'dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz','cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6')
 Q=pin(D/'dim61-67-compatible-sign-families/data/global64.json','8c0c35f8127c578118f1f09e5950d712a272db32f2d50f382e6fe32f46c6dbf5')
 A=pin(D/'dim63-65-affine-reconstruction/certificate.json','266af3683caadf6d23c2aa3e7409ce28d2d3f4a8644be2784263c80bd407fc7e')
 W=list(map(int,np.load(Z,allow_pickle=False)['words']))
 W += [int(s['support'])for s in json.loads(Q.read_text())['new']]
 W += [int(s['support'])for s in json.loads(A.read_text())['new_classes']]
 assert len(W)==len(set(W))==32801 and all(w.bit_count()==16 for w in W)
 old=json.loads(x.read_text())['supports']['63']
 data=json.loads(y.read_text())
 assert old and data['total']==184415700
 H=data['octad'];P=data['four_support_permutation']
 assert H==sum(1<<j for j in P[:8]) and sorted(P)==list(range(64))
 four={sum(1<<P[i] for i in (x,y,z,x^y^z))for x,y,z in it.combinations(range(64),3)}
 assert len(four)==10416
 # Old full four-support sign blocks are restricted to the same balanced hyperplane.
 f={}
 for b in sorted(four):
  J=pos(b);n=(b&H).bit_count()
  if n%2==0:
   f[b]={physical(s,J)for s in range(16)if (physical(s,J)&H).bit_count()*2==n}
 assert sum(map(len,f.values()))==102900
 A=phase();norm=[]
 for s,M in old:
  U=pos(s&H);V=pos(s^ (s&H)); t=iter(U);w=iter(V)
  J=[next(t)if i in M else next(w) for i in range(9)]
  signs=[physical(a,J)for a in A if (physical(a,J)&H).bit_count()*2==len(U)]
  norm.append((s,signs))
 assert sum(len(a)for _,a in norm)==2084
 assert all(len(set(a))==len(a)for _,a in norm)
 return W,f,H,norm
