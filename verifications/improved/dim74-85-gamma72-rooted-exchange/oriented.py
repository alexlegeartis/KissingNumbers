"""Orientation-preserving tail replay for the 74–76 latitude variants."""
from collections import Counter
from fractions import Fraction as F
from math import ceil
import numpy as np
import ico

def edges(parent,d):
 if d==75:
  return {j:np.asarray(v,dtype='i8')for j,v in ico.edges(ico.palette()).items()}
 R=np.loadtxt(parent/('polygon13.txt' if d==74 else 'rays.txt'),dtype=object)
 R=np.asarray([[int(v)for v in x]for x in R],dtype=object);Y=R@R.T;XX=np.diag(Y)
 assert R.shape=={74:(13,2),76:(246,4)}[d] and all(v>0 for v in XX)
 return {j:np.asarray((Y>0)&(100*Y*Y>(25-4*j)**2*(XX[:,None]*XX[None,:])),dtype='i8') if 4<=j<=6 else np.full(Y.shape,j==10,dtype='i8')for j in (-10,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,10)}

def replay(parent,state,row,key,model):
 d=row['d'];K,pos,den,z,spec=state
 E=edges(parent,d)
 M=len(E[0]);assert M=={74:13,75:60,76:246}[d] and (E[4]==np.eye(M)).all()
 assert F(16,25)*10+F(8,5)==8 and F(4,5)*5==4
 assert F(16,25)*3+F(8,5)<4 and F(8,5)*8<16
 assert F(16,25)*10-F(8,5)>4
 assert all(k.split(',')[1]=='10'for k in row['free'])
 J=(4,5,6);N={j:model.P[j]//2 for j in J}
 degrees=[M-1+sum(int(E[j][i].sum())*N[j]for j in J) for i in range(M)]
 assert degrees==row['degrees']
 final=F(0)
 for i,(D,w) in enumerate(zip(degrees,row['weights'],strict=True)):
  O=[0]*len(z);c=(M-1)*D
  for j in J:
   P,E0=spec[j];bins=P[10]+E0[10]
   L={x:E[x][i]for x,y in bins};R={y:E[j][i]@E[y]for x,y in bins}
   H={(x,y):int(L[x]@R[y]) for x,y in bins}
   for x,y in P[10]:O[pos[key(j,10,x,y)]]+=N[j]*H[x,y]
   c+=N[j]*sum(H[x,y]for x,y in E0[10])
  V=tuple(sum(w*r[t]for w,r in zip(O,z)if w) for t in range(4))
  ix,weights=w;u=list(map(F,weights))
  assert len(ix)==len(u)==len(set(ix)) and all(0<=j<len(z)for j in ix) and min(u,default=0)>=0
  H=[F(V[t])+sum(y*z[j][t]for j,y in zip(ix,u)) for t in range(4)]
  D2=F(H[3]+model.N*sum(max(0,x)for x in H[:3]),den)
  assert D2>=0
  final=max(final,D2+c)
 T=ceil(final/2);S=ceil(sum(F(model.N*r,D+1)*model.low_factor(D,T)for D,r in Counter(degrees).items()))
 assert (T,S)==(row['T'],row['S'])
 return T,S,model.N8+S+{74:6,75:12,76:24}[d]
