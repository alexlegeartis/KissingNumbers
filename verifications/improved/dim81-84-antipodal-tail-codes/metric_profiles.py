"""Separate horizontal 4+8 and 6+6 profiles by actual selected intersections."""
import numpy as np
from functools import lru_cache
ALL24=(1<<24)-1

def improve(R, rows, clips, dimension=None,metric=None):
 H=R['H']; RB=R['RB']; AA6=R['A6']; layout=R['layouts']; pd=R['pack24']
 V8=R['V8']; V8F=(V8@R['F']).astype('f8'); W8=V8.astype('f8')
 V6=R['V6']; V6F=V6@R['F']; HH=R['HH']
 A6=np.array([a for a,b in AA6],dtype='u8');B6=np.array([RB[b] for a,b in AA6],dtype='u8')
 A8=np.array([a for a,b in R['A8']],dtype='u8')
 key48={};key66={}; perms48={};perms66={}
 for t in H:
  U,V=layout[t]
  key48[t]=(U[:,0].astype('u8')<<36)|(A8[U[:,2]]<<12)|U[:,1].astype('u8')
  key66[t]=(A6[V[:,0]]<<36)|(A6[V[:,1]]<<12)|B6[V[:,0]]
  for keys,perms in ((key48,perms48),(key66,perms66)):
   c=keys[t];p=np.argsort(c);assert len(np.unique(c))==len(c)
   perms[t]=(c[p],p)

 def match(keys,perms,t,u,shift):
  Q=keys[t]^shift; val,order=perms[u]
  loc=np.minimum(np.searchsorted(val,Q),len(val)-1)
  i=np.flatnonzero(val[loc]==Q)
  return i,order[loc[i]]

 masks={}
 def matrix(a,b):
  # Every coordinate and product is an integer with sum of absolute
  # summands < 2^53 (same V8-on-V8 bound as pinned parent).
  ij=(a,b)
  if ij not in masks:
   M=np.abs(V8F[a]@W8[b].T)==4
   masks[ij]=pd(M)
   masks[b,a]=pd(M.T)
  return masks[ij]

 @lru_cache(640)
 def links(t,u):
  d=RB[t^u];assert d
  Ut,Vt=layout[t]; Uu,Vu=layout[u]
  a0,b0=match(key48,perms48,t,u,d)
  assert np.all(Ut[a0,2]==Uu[b0,2]) # q8 block unchanged
  a1,b1=match(key48,perms48,t,u,0)
  # Exact changes of the q8 block; q4 block unchanged by the keyed match.
  M=np.array([matrix(i,j)for i,j in zip(Ut[a1,2],Uu[b1,2])],dtype='u4').reshape(-1,24)
  z=[]
  for k in (0,1):
   a,b=match(key66,perms66,t,u,0 if k else d)
   assert np.all(Vt[a,1-k]==Vu[b,1-k])  
   f=np.abs((V6F[Vt[a,k]]*V6[Vu[b,k]]).sum(axis=1))==4
   z.append((a[f],b[f]))
  return a0,b0,a1,b1,M,z

 occ48=[];occ66=[]
 for (v,c,t),bs in zip(rows,clips):
  U,V=layout[t]; r1,r2,r3=U.T
  m=np.full(len(U),ALL24,dtype='u4'); n=np.ones(len(V),bool)
  for delta in set(bs):
   p=HH[delta]
   m[R['MD0'][r1,r2^RB[delta]]]=0
   m&=np.uint32(ALL24)^R['BAD8'][p,r3]
   n&=~R['BAD6'][p,V].any(axis=1)
  occ48.append(m);occ66.append(n)
 def size(a,b):return sum(int(x).bit_count()for x in a)+int(np.count_nonzero(b))
 base=sum(size(a,b)for a,b in zip(occ48,occ66))
 V=np.array([q[0]for q in rows]);colors=np.array([q[1]for q in rows]);G=V@V.T if metric is None else metric
 E=(colors[:,None]==colors[None,:])&(G!=0);np.fill_diagonal(E,False)
 env=[np.flatnonzero(r)for r in E]
 order=sorted(range(len(rows)),key=lambda i:(size(occ48[i],occ66[i]),i))
 for i in order:
  t=rows[i][2]
  m=np.full(len(occ48[i]),ALL24,dtype='u4')^occ48[i]
  n=~occ66[i]
  for j in env[i]:
   if not np.any(m) and not np.any(n):break
   u=rows[j][2];a0,b0,a1,b1,M,z=links(t,u)
   m[a0]&=np.uint32(ALL24)^occ48[j][b0]
   if len(a1):
    q=(M&occ48[j][b1,None])!=0
    m[a1]&=np.uint32(ALL24)^pd(q)
   for a,b in z:n[a]&=~occ66[j][b]
  occ48[i]|=m;occ66[i]|=n

 # Full literal replay of edge compatibility, now based only on final occupancy.
 for i in range(len(rows)):
  t=rows[i][2]
  for j in env[i]:
   if j>=i:break
   u=rows[j][2];a0,b0,a1,b1,M,z=links(t,u)
   assert not np.any(occ48[i][a0]&occ48[j][b0])
   if len(a1):
    q=(M&occ48[j][b1,None])!=0
    assert not np.any(occ48[i][a1]&pd(q))
   for a,b in z:assert not np.any(occ66[i][a]&occ66[j][b])
 end=sum(size(a,b)for a,b in zip(occ48,occ66))
 import os
 if os.environ.get('EXPORT_PROFILE'):
  from pathlib import Path
  p=Path(os.environ['EXPORT_PROFILE'])
  Ainv=np.array(sorted(R['RA'],key=R['RA'].get),dtype='u4')
  Binv=np.array(sorted(RB,key=RB.get),dtype='u4')
  common=p/'profile_stock.npz'
  if not common.exists():
   U0=np.array([layout[t][0]for t in H]);U0[:,:,0]=Ainv[U0[:,:,0]];U0[:,:,1]=Binv[U0[:,:,1]]
   np.savez(common,H=H,F=R['F'],A8=R['A8'],A6=R['A6'],V8=V8,V6=V6,U=U0,V=np.array([layout[t][1]for t in H]))
  np.savez(p/f'profile_{dimension}.npz',occ48=np.array(occ48),occ66=np.array(occ66),tails=np.array([a for a,b,c in rows]),colors=np.array([b for a,b,c in rows]),phase=np.array([HH[c]for a,b,c in rows]),base=base,end=end)
 return base,end,len(masks),links.cache_info()
