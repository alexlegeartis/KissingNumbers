"""Exact replay of twelve harmonic Γ72 latitudes."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
from math import ceil
import gzip,json,sys
import numpy as np
import model as p
from common import check_profiles,rn
import ico
HERE=Path(__file__).resolve().parent
def rational_edges(path,d):
 R=np.loadtxt(path,dtype=object)
 V=np.asarray([[int(x)for x in row]for row in R],dtype=object)
 dots=V@V.T;sq=np.diag(dots)
 assert V.shape=={74:(13,2),76:(246,4)}[d] and all(x>0 for x in sq)
 return {j:np.asarray((dots>0)&(100*dots*dots>(25-4*j)**2*(sq[:,None]*sq[None,:])),dtype='i8') if 4<=j<=6 else np.full(dots.shape,j==10,dtype='i8')for j in (-10,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,10)}
def lp_bound(C,P,D):
 _,_,K,piv,Q,G,B,U=P
 y=list(map(F,D));assert len(y)==len(B) and min(y)>=0
 base=sum(C[k]*b for k,b in zip(piv,B))
 coeff=[C[k]-sum(C[l]*row[i]for l,row in zip(piv,G))for i,k in enumerate(Q)]
 residual=[c-sum(t*row[i]for t,row in zip(y,G))for i,c in enumerate(coeff)]
 return base+sum(t*b for t,b in zip(y,B))+sum(max(F(0),d)*u for d,u in zip(residual,U))
def costvec(V,P):
 A,_,K,*_=P;loc={z:i for i,z in enumerate(K)};C=[0]*len(K)
 for x,z in A[10]:C[loc[10,*rn(x,z)]]+=V[x,z]
 return tuple(C)
def oriented(d,profiles,duals,row):
 E=ico.edges(ico.palette()) if d==75 else rational_edges(HERE/('polygon13.txt' if d==74 else 'rays.txt'),d)
 M=len(E[0]); E={j:np.asarray(A,dtype='i8')for j,A in E.items()}
 assert M=={74:13,75:60,76:246}[d] and (E[4]==np.eye(M)).all()
 assert F(16,25)*10+F(8,5)==8 and F(4,5)*5==4
 assert F(16,25)*3+F(8,5)<4 and F(8,5)*8<16
 assert F(16,25)*10-F(8,5)>4
 D=[M-1+sum((p.P[j]//2)*int(E[j][a].sum())for j in (4,5,6))for a in range(M)]
 assert D==row['degrees'];T2=[F((M-1)*v)for v in D]
 for j,P in profiles.items():
  A,S=P[:2];pairs=set(A[10]+S[10]);cache={}
  for a in range(M):
   r={x:E[x][a]for x,z in pairs}; s={z:E[j][a]@E[z]for x,z in pairs}
   V={(x,z):int(r[x]@s[z])for x,z in pairs}; C=costvec(V,P)
   if C not in cache:cache[C]=lp_bound(C,P,duals[str(j)][str(a)])
   T2[a]+=(cache[C]+sum(V[z]for z in S[10]))*(p.P[j]//2)
 T=ceil(max(T2)/2)
 s=ceil(sum(F(p.N*c,d+1)*p.low_factor(d,T)for d,c in Counter(D).items()))
 assert (T,s)==(row['triangles'],row['selected'])
def projective(d,profiles,duals,row):
 k=d-72;sh,X=p.PALETTES[d];ts=p.tails_of(k,sh);M=len(ts)
 r=np.array(list(sh)+[0]*(k-len(sh)),dtype='i8')
 norm=sum(t*t for t in sh);Q=8-10*X
 assert (ts*ts).sum(1).tolist()==[norm]*M
 assert M==len(set(map(tuple,ts))) and set(map(tuple,ts))==set(map(tuple,-ts))
 assert 0<Q and 25*X<=16 and 3*X+Q<=4
 assert 4*(sh[0]+sh[1])**2*Q/norm<=16
 H=X.numerator*norm;B=8*X.denominator-10*X.numerator;L=4*X.denominator*norm
 def gate(j,dot):return np.abs(H*j+B*dot)>L
 a=ts@r; same=int(gate(10,a).sum())-1
 D=same+sum(int(gate(j,a).sum())*(p.P[j]//2)for j in (4,5,6))
 assert D==row['degree'];T2=F(same*D)
 # Stabilizer of r: permute equal positive coordinates or permute and
 # change signs of all zero coordinates. It is transitive on each key.
 marks=[np.flatnonzero(r==x)for x in sorted(set(sh))]
 def orbit(t):return tuple(tuple(sorted(t[ix]))for ix in marks)+(tuple(sorted(np.abs(t[len(sh):]))),)
 for j,P in profiles.items():
  A,S=P[:2];pairs=set(A[10]+S[10]);reps={}
  for t in ts[gate(j,a)]:
   z=orbit(t)
   if z not in reps:reps[z]=[t,0]
   reps[z][1]+=1
  V={q:0 for q in pairs}; U={x:gate(x,a)for x,y in pairs}
  for t,m in reps.values():
   dt=ts@t;W={y:gate(y,dt)for x,y in pairs}
   for x,y in pairs: V[x,y]+=m*int(np.count_nonzero(U[x]&W[y]))
  assert all(V[x,y]==V[-x,-y] for x,y in pairs)
  C=costvec(V,P)
  T2+=(lp_bound(C,P,duals[str(j)])+sum(V[z]for z in S[10]))*F(p.P[j]//2,2)
 T=ceil(T2/2);assert T==row['triangles']
 s=2*ceil(F(p.N//2*M,D+1)*p.low_factor(D,T));assert s==row['selected']
def verify(dims=()):
 p.independently_replay_parent_counts()
 profiles={j:check_profiles(j) for j in (4,5,6)}
 with gzip.open(HERE/'duals.json.gz','rt')as f:duals=json.load(f)
 ans=[]
 for row in json.loads((HERE/'results.json').read_text()):
  d=row['dimension'];k=d-72
  if dims and d not in dims:continue
  (oriented if d<=76 else projective)(d,profiles,duals[str(d)],row)
  add=6 if d==74 else 2*k*(k-1)
  s=p.N8+row['selected']+add
  assert s==row['lower_bound'] and s-row['previous']==row['gain']>0
  ans.append(s);print(d,s,row['gain'],flush=True)
 return ans
if __name__=='__main__':
 assert __debug__
 verify(tuple(map(int,sys.argv[1:])))
