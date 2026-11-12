"""A genuinely different actual normal subsystem: 2A4 inside the retained 2D5.

All numerical LP and head-search outputs are discovery, not certificates.
"""
from pathlib import Path
from itertools import product,permutations
from functools import lru_cache
from fractions import Fraction as F
from math import prod,gcd
import json,sys,numpy as np
from scipy.optimize import linprog
from scipy.linalg import qr

HERE=Path(__file__).parent
ROOT=HERE.resolve().parents[2]
N=6218175600

def build():
 G=np.diag([8]*4)+np.diag([-4]*3,1)+np.diag([-4]*3,-1)
 E=json.loads((ROOT/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json').read_text())
 U=np.array(E['normals'],dtype=np.int64)[[1,0,3,4]];B=np.array(E['catalogue_gram'],dtype=np.int64)
 assert np.array_equal(U@B@U.T,G)
 aset=[a for a in product(range(5),repeat=5) if min(a)==0]
 cells=np.array([tuple(a[i]-a[i+1] for i in range(4)) for a in aset],dtype=np.int64)
 def key(a):
  aa=tuple(sorted(a));bb=tuple(sorted(max(a)-v for v in a));return min(aa,bb)
 keys=sorted({key(a) for a in aset});gid={a:i for i,a in enumerate(keys)};groupid=np.array([gid[key(a)] for a in aset]);groups=[np.flatnonzero(groupid==i).tolist() for i in range(len(keys))]
 z=(np.array(aset)-np.mean(aset,axis=1)[:,None])/2
 q=np.sum(z*z,axis=1)
 roots=[]
 for start in range(4):
  for end in range(start,4):
   for sign in [-1,1]:
    a=np.array([sign if start<=i<=end else 0 for i in range(4)],dtype=np.int64);assert a@G@a==8;roots.append(a@G)
 @lru_cache(None)
 def gaussian(e):
  if not sum(e):return 1
  i=next(i for i,v in enumerate(e) if v);f=list(e);f[i]-=1;total=0
  for j,v in enumerate(f):
   if v:
    f[j]-=1;total+=v*int(G[i,j])*gaussian(tuple(f));f[j]+=1
  return total
 rows=[];rhs=[];tags=[]
 for e in product(range(11),repeat=4):
  deg=sum(e)
  if deg>10 or deg%2:continue
  vals=np.prod(cells**np.array(e),axis=1)
  row=[int(v) for v in np.bincount(groupid,weights=vals,minlength=len(groups))]
  b=F(N*8**(deg//2)*gaussian(e),prod(72+2*j for j in range(deg//2)))-sum(prod(int(x)**a for x,a in zip(r,e)) for r in roots)
  den=b.denominator;row=[v*den for v in row];bb=int(b*den);gg=gcd(*row,bb) or 1
  rows.append([v//gg for v in row]);rhs.append(bb//gg);tags.append(e)
 A=np.array(rows,dtype=float);b=np.array(rhs,dtype=float)/N
 qq=np.array([q[g[0]] for g in groups]);colscale=np.exp(-4.5*qq);colscale/=sum(len(g)*v for g,v in zip(groups,colscale))
 scales=np.maximum(np.max(np.abs(A*colscale[None,:]),axis=1),np.abs(b));scales[scales==0]=1
 As=A*colscale[None,:]/scales[:,None];bs=b/scales
 _,R,piv=qr(As.T,pivoting=True,mode='economic');rank=np.linalg.matrix_rank(As);choose=sorted(piv[:rank].tolist())
 As=As[choose];bs=bs[choose];scales=scales[choose];rows=[rows[i] for i in choose];rhs=[rhs[i] for i in choose];tags=[tags[i] for i in choose]
 base=linprog(np.zeros(len(groups)),A_eq=As,b_eq=bs,bounds=(0,None),method='highs',options={'threads':1})
 assert base.success,base.message
 print('A4 MODEL cells',len(cells),'groups',len(groups),'rank',rank,flush=True)
 return dict(gram=G,cells=cells,z=z,q=q,groups=groups,groupid=groupid,rows=rows,rhs=rhs,tags=tags,As=As,bs=bs,scales=scales,colscale=colscale,normals=U.tolist())

def solve(M,chosen):
 obj=np.bincount(M['groupid'][chosen],minlength=len(M['groups']))
 r=linprog(obj*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1});assert r.success,r.message
 return r.fun*N,r,obj

if __name__=='__main__':
 M=build();zero=np.flatnonzero(M['q']<1e-8)
 print('central',solve(M,zero)[0])
 for j,g in enumerate(M['groups']):
  a=np.zeros(len(M['groups']));a[j]=1
  lo=linprog(a*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1})
  hi=linprog(-a*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1})
  print(j,'q',M['q'][g[0]],'orbit',len(g),'range',lo.fun*N,-hi.fun*N)
