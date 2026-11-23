"""Replay exact APN multipartite D_m fibres of Γ72."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations, permutations
import json, sys, runpy
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'dim73-82-gamma72-depth6-lattice'
pins={
  "verify.py": "35b943ae846ba62a0270f1b3d5ff4f5280c208de4f4c2c5dd5df24ce3bc787df",
  "theta.py": "168d5d052ce0a6e62085800a1fa956e88791519349bab7309409dec85420f073",
  "certificate.json": "e5e656d728ef6ca1a7242e2863c8c41447247a9165a7b5aa862ea15d11a1584e",
  "results.json": "f4042a060aad10f6a3766aff8347e05fd5a9541a517516d3af3f00e49d1b5ccf",
  "comparisons.json": "d04816fd2d1cd4fa1efcc880245a3fcd104d0254c74dfaebe2b1dcf7cdeb7e4d",
  "catalog_gram.json": "1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab"
}
assert all(sha256((base/f).read_bytes()).hexdigest()==h for f,h in pins.items())
sys.path.insert(0,str(base))
P=runpy.run_path(str(base/'verify.py'))
W=json.loads((here/'witness.json').read_text())
M=P['M4'];Or=P['orth'];K=(23214,668618,682340)
a,b,c=K;assert a^b==c
def span(B):
 L=[0]
 for x in B:L+=[x^y for y in L]
 return L
def basis(A):
 E={}
 for u in A:
  while u:
   k=u.bit_length()-1
   if k not in E:E[k]=u;break
   u^=E[k]
 return list(E.values())
AB=span(basis(M));RA={v:i for i,v in enumerate(AB)}
Blist=span(basis(sum(int(t)<<j for j,t in enumerate(v)) for v in P['B']))
assert len(AB)==len(Blist)==4096
H=span(W['B_basis']);assert len(H)==len(set(H))==128 and set(H)<=set(Blist)
for ts in H:
 assert all(Or(u,ts)for u in K)
 assert all(not any(x^ts in M[v]for x in M[u])for u,v in combinations(K,2))
assert all(Or(v,x)for u in K for x in M[u]for v in K)
assert P['prove']([(w>>i)&1 for w in (a,b,0)for i in range(24)],6) is None
RB={b:i for i,b in enumerate(Blist)}
HE=np.array([RB[h]for h in H],dtype='u2'); HG=np.zeros(4096,dtype='bool'); HG[HE]=1
a,b=K[:2]
products=[]
MD={i:np.array([RB[x]for x in M[h]],dtype='u2')for i,h in enumerate(AB)if h}
MD0=np.zeros((4096,4096),dtype='bool')
for i,a1 in MD.items():MD0[i,a1]=True
for h in M:
 if h==a or h==b:continue
 A0,A1,A2=RA[h^a],RA[h^b],RA[h]
 M0,M1,M2=MD[A0],MD[A1],MD[A2]
 net=M0[:,None,None]^M1[None,:,None]^M2[None,None,:]
 i,j,k=np.where(HG[net])
 if len(i):
  products.append(np.column_stack([net[i,j,k],np.full(len(i),A2,dtype='u2'),M0[i],M1[j],M2[k]]).astype('u2'))
B0=np.vstack(products)
so=np.argsort(B0[:,0],kind='stable');B0=B0[so]
T={}
for r in HE:
 a0=np.searchsorted(B0[:,0],r);b0=np.searchsorted(B0[:,0],r,side='right')
 V=B0[a0:b0,1:]
 assert len(V)==13824,(r,len(V))
 A=np.column_stack([V[:,0]^RA[a],V[:,0]^RA[b],V[:,0]])
 T[int(r)]=(A,V[:,1:])

def mul16(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&16:a^=19
 return out
def decode(x):
 return x&3,(x>>2)|(mul16(mul16(x,x),x)<<2)
comparison={81:6228705122,82:6229251942,83:6230072172,84:6231712632,85:6235993850}
out={}
for record in W['rows']:
 m=record['dimension']-72; X=record['vertices']; Ds=record['rays']
 assert len(X)==len(set(X))==m and all(0<=x<16 for x in X)
 S=[]
 for i,j in combinations(range(m),2):
  c0,b0=decode(X[i]);c1,b1=decode(X[j])
  if c0!=c1:
   for e in range(2):S.append([i,j,e,c0^c1,H[b0^b1^(64 if e else 0)]])
 assert S==[v[:5]for v in Ds]
 assert len(set(tuple(s[3:])for s in S))==len(S)
 assert all(len(set(v[6]))==len(v[6]) and all(0<=t<len(Ds)for t in v[6])for v in Ds)
 for i,j in combinations(range(len(Ds)),2):
  x,y=Ds[i],Ds[j]
  is_block=(x[3]==y[3] and len(set(x[:2])&set(y[:2]))==1)
  assert is_block==(j in x[6] or i in y[6])
 assert all(i not in v[6]for i,v in enumerate(Ds))
 total=0
 for i,v in enumerate(Ds):
  _,_,_,_,t,n,bl=v
  if not bl:assert n==136704
  else:
   A,B=T[RB[t]]
   bad=np.zeros(len(A),dtype='bool')
   for u in bl:
    assert Ds[u][3]==v[3] and len(set(v[:2])&set(Ds[u][:2]))==1
    d=RB[Ds[u][4]^t]
    assert d and d in T
    bad|=MD0[A,B^d].any(axis=1)
   assert n==8*(len(A)-int(bad.sum()))
  total+=n
 d=72+m;lower=6218175600+2*total+2*m*(m-1)
 assert lower==record['lower_bound'] and d in comparison and lower>comparison[d]
 out[str(d)]={'lower_bound':lower,'prior':comparison[d],'gain':lower-comparison[d]}
 print('PASS',d,lower,'gain',lower-comparison[d],'antipodal rays',len(Ds),flush=True)
assert out==W['results']
