"""Independent exact replay, no numerical optimization or floating geometry."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
from math import ceil
from sympy.polys.matrices import DomainMatrix
import gzip,json
import numpy as np
import model as p
from harmonic import system,rn
HERE=Path(__file__).resolve().parent

def check_profiles(j):
 # Weighted theta: Theta_H/Delta^4 has weight deg(H)−12.  At
 # degrees k=12,14,16,18,20,22 the basis is 1,none,E4,E6,E4²,E4E6.
 # q^1 coefficients of E4,E6 are 240,−504; Delta^4=q^4−96q^5+...
 P,E,K,rows=system(j)
 A=DomainMatrix.from_list_sympy(len(rows),len(rows[0]),rows).to_dense()
 rw,d,piv=A.rref_den(method='FF')
 L=[[int(y)for y in row]for row in rw.to_Matrix().tolist() if any(row)]
 assert len(piv)=={4:114,5:106,6:94}[j] and max(piv)<len(K)
 Q=[i for i in range(len(K))if i not in piv]; d=int(d)
 assert all(L[i][k]==d for i,k in enumerate(piv))
 G=[tuple(F(r[k],d)for k in Q)for r in L]
 B=[F(r[-1],d)for r in L]
 U=[F(p.N if K[k][0]==10 else p.N8)for k in Q]
 return P,E,K,piv,Q,G,B,U

def replay():
 p.independently_replay_parent_counts()
 assert F(16,25)*10+F(8,5)==8 and F(4,5)*5==4
 assert F(16,25)*3+F(8,5)<4 and F(8,5)*8<16
 assert F(16,25)*10-F(8,5)>4
 R=np.loadtxt(HERE/'rays.txt',dtype=object);R=np.array([[int(x)for x in r]for r in R],dtype=object)
 dots=R@R.T; norms=np.diag(dots)
 M=len(R);assert R.shape==(246,4) and all(x>0 for x in norms)
 edges={j:np.asarray((dots>0)&(100*dots*dots>(25-4*j)**2*(norms[:,None]*norms[None,:])),dtype='i8') if 4<=j<=6 else np.full(dots.shape,j==10,dtype='i8')for j in (-10,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,10)}
 assert (edges[4]==np.eye(M)).all()
 D=[M-1+sum((p.P[j]//2)*int(edges[j][a].sum())for j in (4,5,6))for a in range(M)]
 T2=[F((M-1)*d)for d in D]
 with gzip.open(HERE/'duals.json.gz','rt') as f:cert=json.load(f)
 assert cert['degrees']==D
 for j in (4,5,6):
  P,E,K,piv,Q,G,B,U=check_profiles(j)
  pairs=set(P[10]+E[10]); loc={z:i for i,z in enumerate(K)}; cache={}
  for a in range(M):
   r={x:edges[x][a]for x,z in pairs}; s={z:edges[j][a]@edges[z]for x,z in pairs}
   V={(x,z):int(r[x]@s[z])for x,z in pairs}
   C=[0]*len(K)
   for x,z in P[10]:C[loc[10,*rn(x,z)]]+=V[x,z]
   C=tuple(C)
   if C not in cache:
    y=list(map(F,cert['duals'][str(j)][str(a)]));assert len(y)==len(B) and min(y)>=0
    base=sum(C[k]*b for k,b in zip(piv,B))
    coeff=[C[k]-sum(C[l]*row[i]for l,row in zip(piv,G))for i,k in enumerate(Q)]
    residual=[c-sum(t*row[i]for t,row in zip(y,G))for i,c in enumerate(coeff)]
    cache[C]=base+sum(t*b for t,b in zip(y,B))+sum(max(F(0),d)*u for d,u in zip(residual,U))
   T2[a]+=(cache[C]+sum(V[z]for z in E[10]))*(p.P[j]//2)
  print('verified joint systems and duals for j',j,flush=True)
 T=ceil(max(T2)/2);assert T==cert['triangles']
 selected=ceil(sum(F(p.N*c,d+1)*p.low_factor(d,T)for d,c in Counter(D).items()))
 assert selected==cert['selected']==459824
 result=p.N8+selected+24
 assert result==6218635448
 out={'dimension':76,'previous':6218400136,'lower_bound':result,'gain':result-6218400136,'tails':M,'selected':selected,'triangles':T}
 assert json.loads((HERE/'results.json').read_text())==out
 print(json.dumps(out));return out
if __name__=='__main__':
 assert __debug__
 replay()
