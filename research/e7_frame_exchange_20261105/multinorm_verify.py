"""Exact local capacity for norm-3 and norm-18 cheap candidates, with integer-shell candidates included."""
from pathlib import Path
import numpy as np,time
p=Path(__file__).resolve().parent/'cache';z=np.load(p/'candidates14.npz');A=z['Q'].astype(np.int16);AB=z['B'];old=z['old'].astype(np.int16);sel=(AB!=65535).sum(1)==4;A=A[sel];AB=AB[sel,:4];w=np.load(p/'multinorm9.npz');Q=w['Q'].astype(np.int16);B=w['B'];norm=w['norm'];X=w['X'].astype(np.int16);cost=(B!=65535).sum(1);R=Q[norm==3];RB=B[norm==3];T=Q[(norm==18)&(cost==4)];TB=B[(norm==18)&(cost==4)];start=time.time();affected=np.unique(TB[TB!=65535]);canonical=affected[affected<966];assert len(canonical)*2==len(affected)==220
assert set(map(int,RB[RB!=65535]))<=set(map(int,affected))
assert set(map(int,AB.ravel()))<=set(map(int,canonical))

def independent(G,k):
 n=len(G)
 if n<k:return False
 D=~G.copy();np.fill_diagonal(D,False);N=[sum(1<<int(j) for j in np.flatnonzero(D[i])) for i in range(n)]
 def rec(P,r):
  if r==0:return True
  if P.bit_count()<r:return False
  while P:
   b=P&-P;u=b.bit_length()-1;P^=b
   if rec(P&N[u],r-1):return True
  return False
 return rec((1<<n)-1,k)
for v in canonical:
 ids=np.flatnonzero((AB==v).any(1));P=A[ids].copy();P[(P@old[v])<0]*=-1;ti=np.flatnonzero((TB==v).any(1));ri=np.flatnonzero((RB==v).any(1));V=np.r_[P,T[ti]];n=np.r_[np.full(len(P),8),np.full(len(ti),18)];G=V@V.T>np.floor(np.sqrt(n[:,None]*n[None,:])/2);np.fill_diagonal(G,False);assert not independent(G,5)
 for j in ri:
  allow=V@R[j]<=np.floor(np.sqrt(n*3)/2);assert not independent(G[np.ix_(allow,allow)],3)
 if len(ri)==2 and R[ri[0]]@R[ri[1]]<=1:
  allow=(V@R[ri[0]]<=np.floor(np.sqrt(n*3)/2))&(V@R[ri[1]]<=np.floor(np.sqrt(n*3)/2));assert not independent(G[np.ix_(allow,allow)],2)
assert all(np.count_nonzero((RB==v).any(1))<=2 for v in range(len(X)))
print('PASS exact weighted local capacity for',len(canonical),'canonical affected points; seconds',round(time.time()-start,1))
