from itertools import combinations,product
from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
S=[tuple(i for i,v in enumerate(b)if v=='1')for b in (P/'best51.txt').read_text().split()]
assert len(S)==len(set(S))==51 and all(len(s)==4 for s in S)
assert all(len(set(s)&set(t))<=2 for s,t in combinations(S,2))
A=[]
for i in range(12):
 u=[0]*12;u[i]=2;A.append(tuple(u))
for s in S:
 for t in product((-1,1),repeat=3):
  u=[0]*12
  for i,h in zip(s,(1,*t)):u[i]=h
  A.append(tuple(u))
A=np.asarray(sorted(A),dtype='i8');G=A@A.T
n=len(A);assert n==420 and np.all(np.diag(G)==4) and np.abs(G-4*np.eye(n,dtype='i8')).max()<=2
H=(G!=0).astype('i8')-np.eye(n,dtype='i8')
deg=H.sum(axis=1)
print('full exact P12a',n,'degree range',deg.min(),deg.max(),'mean',deg.mean())
np.save(P/'P12a_gram.npy',G)
(P/'P12a_tails.json').write_text(json.dumps(A.tolist()))
I=np.ones(n,dtype=bool)
for k in range(n-384):
 a=np.flatnonzero(I)
 j=a[np.argmax(deg[a])];I[j]=0;deg-=H[j]
assert I.sum()==384
# Strictly edge-reducing active/inactive swaps.
while True:
 a=np.flatnonzero(I);b=np.flatnonzero(~I)
 diff=deg[a][:,None]-deg[b][None,:]+H[a[:,None],b[None,:]]
 i,j=np.unravel_index(np.argmax(diff),diff.shape)
 if diff[i,j]<=0:break
 x,y=a[i],b[j];I[x]=0;I[y]=1;deg+=H[y]-H[x]
nd=(H[I][:,I]).sum(axis=1)
print('thin exact P12a',I.sum(),'degree range',nd.min(),nd.max(),'mean',nd.mean())
np.save(P/'P12a_384_gram.npy',G[I][:,I])
(P/'P12a_384_tails.json').write_text(json.dumps(A[I].tolist()))
(P/'P12a_384_ids.json').write_text(json.dumps(np.flatnonzero(I).tolist()))
