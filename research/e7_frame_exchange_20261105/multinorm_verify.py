"""Exact local capacity for norm-3 and norm-18 cheap candidates, with integer-shell candidates included."""
from pathlib import Path
import numpy as np,time
from scipy.sparse.csgraph import connected_components
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
local={}
for v in canonical:
 ids=np.flatnonzero((AB==v).any(1));P=A[ids].copy();P[(P@old[v])<0]*=-1;ti=np.flatnonzero((TB==v).any(1));ri=np.flatnonzero((RB==v).any(1));V=np.r_[P,T[ti]];n=np.r_[np.full(len(P),8),np.full(len(ti),18)];G=V@V.T>np.floor(np.sqrt(n[:,None]*n[None,:])/2);np.fill_diagonal(G,False);assert not independent(G,5)
 for j in ri:
  allow=V@R[j]<=np.floor(np.sqrt(n*3)/2);assert not independent(G[np.ix_(allow,allow)],3)
 local[int(v)]=(np.r_[V,R[ri]],np.r_[n,np.full(len(ri),3)])
 if len(ri)==2 and R[ri[0]]@R[ri[1]]<=1:
  allow=(V@R[ri[0]]<=np.floor(np.sqrt(n*3)/2))&(V@R[ri[1]]<=np.floor(np.sqrt(n*3)/2));assert not independent(G[np.ix_(allow,allow)],2)
assert all(np.count_nonzero((RB==v).any(1))<=2 for v in range(len(X)))
high=Q[(norm==18)&(cost==9)];HB=B[(norm==18)&(cost==9)];free=np.unique([int(v) for row in HB for v in row if v!=65535 and int(v) not in set(affected)]);canonical_free=free[free<966];assert len(canonical_free)*2==len(free)==288
assert all(sum(int(v) in set(free) for v in row if v!=65535)==6 for row in HB)
for v in canonical_free:
 ids=np.flatnonzero((HB==v).any(axis=1));assert len(ids)==72;P=high[ids];G=P@P.T>9;np.fill_diagonal(G,False);k,lab=connected_components(G);assert k==2 and sorted(np.bincount(lab))==[24,48]
 for component in range(k):
  H=G[np.ix_(lab==component,lab==component)];cap=4 if len(H)==24 else 8;assert not independent(H,cap+1)
# Each nine-conflict candidate forces three cheap-source witnesses. At each
# such source point, compatibility with it reduces cheap weighted capacity
# from 12 to at most 9. Verify all 10,368 local claims by bitset enumeration.
for q,blockers in zip(high,HB):
 witness=[int(v) for v in blockers if v!=65535 and int(v) in set(affected)];assert len(witness)==3
 for v in witness:
  base=v%966;P,n=local[base];P=P if v<966 else -P
  allow=P@q<=np.floor(np.sqrt(n*18)/2);P=P[allow];n=n[allow];P4=P[n!=3];n4=n[n!=3];P3=P[n==3];G=P4@P4.T>np.floor(np.sqrt(n4[:,None]*n4[None,:])/2);np.fill_diagonal(G,False)
  assert not independent(G,4)
  for u in P3:
   ok=P4@u<=np.floor(np.sqrt(n4*3)/2);assert not independent(G[np.ix_(ok,ok)],2)
  for i in range(len(P3)):
   for j in range(i):
    if P3[i]@P3[j]>1:continue
    ok=(P4@P3[i]<=np.floor(np.sqrt(n4*3)/2))&(P4@P3[j]<=np.floor(np.sqrt(n4*3)/2));assert not np.any(ok)
print('PASS cheap and high local capacities and all',3*len(high),'high-candidate witness cuts; seconds',round(time.time()-start,1))
