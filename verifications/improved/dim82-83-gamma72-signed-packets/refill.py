"""Literal horizontal conflict replay for append-only 4+4+4, 4+8 and 6+6 packets."""
from functools import lru_cache
import numpy as np
K=24
ALL=np.uint32((1<<K)-1)
def pack24(X):
 q=np.packbits(X,axis=1,bitorder='little').astype('u4')
 return q[:,0]|(q[:,1]<<8)|(q[:,2]<<16)
def sizes(P):
 s,m,n=P
 return (int(s.sum()),sum(int(x).bit_count() for x in m),int(n.sum()))

class Model:
 def __init__(self,H,A8,A6,F,V8,V6,U,V,Q0):
  H=list(map(int,H));self.H=H
  assert len(H)==128 and len(set(H))==128
  self.U=U;self.V=V
  A8=np.asarray(A8,dtype='i8'); A6=np.asarray(A6,dtype='i8')
  U=np.asarray(U,dtype='i8'); V=np.asarray(V,dtype='i8');Q0=np.asarray(Q0,dtype='i8')
  self.N=(Q0.shape[1],U.shape[1],V.shape[1])
  assert self.N[0]==13824 and 24*self.N[1]+self.N[2]==6528
  A=np.unique(np.r_[U[:,:,0].ravel(), A8[:,0], A6[:,0], Q0[:,:,0].ravel()])
  B=np.unique(np.r_[U[:,:,1].ravel(), A8[:,1], A6[:,1], Q0[:,:,1:].ravel(),H])
  # every original B word lies in the common 4096 element binary space;
  # form its complete basis, then prove the numeric compression is linear.
  S={0}
  for h in B:
   h=int(h)
   if h not in S:S|={k^h for k in S}
  B=np.array(sorted(S),dtype='u4');assert len(B)==4096
  hidx={int(b):int(i) for i,b in enumerate(B)}
  assert all(np.array_equal(B[np.arange(4096)^hidx[h]],B^h) for h in H)
  self.B=B; self.bidx=hidx
  self.delta=np.array([[hidx[h^u] for u in H] for h in H],dtype='u8')
  au=np.searchsorted(A,U[:,:,0]).astype('u8');bu=np.searchsorted(B,U[:,:,1]).astype('u8')
  a6=np.searchsorted(A,A6[:,0]).astype('u8');b6=np.searchsorted(B,A6[:,1]).astype('u8')
  a8=np.searchsorted(A,A8[:,0]).astype('u8')
  aq=np.searchsorted(A,Q0[:,:,0]).astype('u8');bq=np.searchsorted(B,Q0[:,:,1:]).astype('u8')
  self.key=[(aq<<24)|(bq[:,:,0]<<12)|bq[:,:,1], (au<<36)|(a8[U[:,:,2]]<<12)|bu,
            (a6[V[:,:,0]]<<36)|(a6[V[:,:,1]]<<12)|b6[V[:,:,0]]]
  self.p=[np.argsort(x,axis=1) for x in self.key]
  self.z=[np.take_along_axis(x,p,axis=1) for x,p in zip(self.key,self.p)]
  assert all(np.all(x[:,1:]>x[:,:-1]) for x in self.z)
  self.V8=V8.astype('f8');self.R8=(V8@F).astype('f8')
  self.V6=V6;self.R6=V6@F
  assert 24*int(np.max(np.abs(V8)))*int(np.max(np.abs(self.R8)))<2**53
  assert 24*int(np.max(np.abs(V6)))*int(np.max(np.abs(self.R6)))<2**60
 def match(self,k,t,u,shift):
  q=self.key[k][t]^shift
  z=self.z[k][u];p=self.p[k][u]
  j=np.minimum(np.searchsorted(z,q),len(z)-1)
  ii=np.flatnonzero(z[j]==q)
  return ii,p[j[ii]]
 @lru_cache(8192)
 def mat8(self,i,j):
  return pack24(np.abs(self.R8[i]@self.V8[j].T)==4)
 @lru_cache(1024)
 def links(self,t,u):
  d=self.delta[t,u];assert d
  Ut,Vu=self.U[t],self.V[u]
  Uu,Vt=self.U[u],self.V[t]
  a0,b0=self.match(1,t,u,d)
  assert np.all(Ut[a0,2]==Uu[b0,2])
  a1,b1=self.match(1,t,u,np.uint64(0))
  M=np.array([self.mat8(int(i),int(j)) for i,j in zip(Ut[a1,2],Uu[b1,2])],dtype='u4').reshape(-1,24)
  z=[]
  for k in (0,1):
   a,b=self.match(2,t,u,np.uint64(0) if k else d)
   assert np.all(Vt[a,1-k]==Vu[b,1-k])
   f=np.abs((self.R6[Vt[a,k]]*self.V6[Vu[b,k]]).sum(axis=1))==4
   z.append((a[f],b[f]))
  return a0,b0,a1,b1,M,z
 def apron(self,t,u,m,n,P):
  _,b,f=P
  a0,b0,a1,b1,M,z=self.links(t,u)
  m[a0]&=ALL^b[b0]
  if len(a1):
   m[a1]&=ALL^pack24((M&b[b1,None])!=0)
  for a,j in z:n[a]&=~f[j]
 def add(self,i,c,t,rows,occ,G,do_triple=True):
  ns,nm,nn=self.N
  s=np.full(ns,do_triple,bool);m=np.full(nm,ALL,dtype='u4'); n=np.ones(nn,bool)
  for (j,d,u),(q,b,f) in zip(rows,occ):
   if c!=d:continue
   if t==u:
    s&=~q # identical horizontal point, irrespective of tail angle
    m&=ALL^b
    n&=~f
   elif G[i,j]!=0:
    if s.any():
     # Deliberately no cached triple link (48k bytes per pair).
     ii=np.flatnonzero(s)
     e=self.delta[t,u]
     Q=self.key[0][t][ii][None,:]^np.array([0,e,int(e)<<12],dtype='u8')[:,None]
     zz=self.z[0][u];p=self.p[0][u]
     k=np.minimum(np.searchsorted(zz,Q),len(zz)-1)
     s[ii]&=~((zz[k]==Q)&q[p[k]]).any(axis=0)
    if m.any() or n.any():self.apron(t,u,m,n,(q,b,f))
   if not (s.any() or m.any() or n.any()):break
  return s,m,n
