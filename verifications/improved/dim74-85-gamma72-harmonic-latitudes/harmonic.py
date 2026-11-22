import model as p
from fractions import Fraction as F
from math import prod,gcd,lcm
from itertools import product
from functools import lru_cache
from sympy.polys.matrices import DomainMatrix
from sympy import Matrix,ZZ,QQ
import numpy as np,time,json,sys
@lru_cache(None)
def smom(a,b,j,s):
 k=(a+b)//2
 return F((p.N if s==10 else p.N8)*s**k*p.wick(a,b,j),prod(72+2*i for i in range(k)))
@lru_cache(None)
def lap_exp(a,b,j):
 q={}
 if a>=2:q[(a-2,b)]=10*a*(a-1)
 if a and b:q[(a-1,b-1)]=2*j*a*b
 if b>=2:q[(a,b-2)]=10*b*(b-1)
 return q
@lru_cache(None)
def pol(a,b,j,s):
 m={(a,b):1};res={(a,b):F(1)};c=F(1)
 k=a+b
 for l in range(1,k//2+1):
  n={}
  for (u,v),h in m.items():
   for t,q in lap_exp(u,v,j).items():n[t]=n.get(t,0)+h*q
  c*=F(-s,2*l*(2*k-2*l+70));m=n
  res.update((t,c*h)for t,h in m.items())
 fac=lcm(*(v.denominator for v in res.values()))
 return res

def rep(a,b):
 return min((a,b),(-a,-b),(b,a),(-b,-a)) # never assume swap symmetry of pair itself! aggregate objective only
def rn(a,b):
 return min((a,b),(-a,-b))
def system(j,with_swap=False):
 rel={12:-96,14:0,16:144,18:-600,20:384,22:-360}
 r=rep if with_swap else rn
 p10,e10,_,_=p.profile_data(j)
 def q(s,a,b,x,y,z):return 10*(x*x+y*y)+s*z*z+2*(j*x*y+a*x*z+b*y*z)
 E={8:([(4,-4),(-4,4)]if j==6 else []),10:e10}
 P={10:p10,8:[]}
 for a,b in product(range(-5,6),repeat=2):
  if (a,b) in E[8]:continue
  if 640 - 8*j*j-10*(a*a+b*b)+2*j*a*b<0:continue
  if any(0<q(8,a,b,x,y,z)<8 for x,y,z in product(range(-2,3),repeat=3) if z and (x or y)):continue
  P[8].append((a,b))
 def harmonic(a,b,s):
  f=pol(a,b,j,s)
  ds=[sum(v*(F(x)**i)*(F(y)**k)for (i,k),v in f.items())for x,y in P[s]+E[s]]
  return ds[:len(P[s])],sum(ds[len(P[s]):])
 keys=[(s,*r(a,b)) for s in (8,10)for a,b in P[s]]
 K=sorted(set(keys)); ix={x:i for i,x in enumerate(K)}
 rows=[];rhs=[]
 def push(A,t):
  R=[F(0)]*len(K)
  for s,coeffs in A.items():
   for (a,b),c in zip(P[s],coeffs,strict=True):R[ix[s,*r(a,b)]]+=c
  fac=lcm(t.denominator,*(v.denominator for v in R))
  w=[int(x*fac)for x in R];z=int(t*fac);g=gcd(z,*w)
  if g:rows.append([v//g for v in w]);rhs.append(z//g)
 for s in (8,10):
  for k in range(0,11,2):
   for a in range(k+1):
    b=k-a; mm=lambda x,y:x**a*y**b
    push({s:[mm(x,y) for x,y in P[s]]},smom(a,b,j,s)-sum(mm(x,y)for x,y in E[s]))
 for k,ratio in rel.items():
  for a in range(k+1):
   b=k-a
   d,v=harmonic(a,b,8);e,w=harmonic(a,b,10)
   push({8:[-ratio*x for x in d],10:e},F(-w+ratio*v))
 # Translation exact across the entire same coset in projected plane.
 for x,y in product(range(-2,3),repeat=2):
  if x==0 and y==0:continue
  for s,a,b in [(s,a,b)for s in P for a,b in P[s]+E[s]]:
   t=q(s,a,b,x,y,1);c=a+10*x+j*y;d=b+j*x+10*y
   if t not in P:continue
   A={(s,a,b):1,(t,c,d):-1}
   if len(A)<2:continue
   R=[0]*len(K);v=0
   for (s,a,b),w in A.items():
    if (a,b) in E[s]:v-=w
    elif (a,b) in P[s]:R[ix[s,*r(a,b)]]+=w
   if any(R) or v:rows.append(R);rhs.append(v)
 nz=[]
 for rr,u in zip(rows,rhs):
  if any(rr) or u:nz.append((*rr,u))
 nz=list(dict.fromkeys(nz));return P,E,K,nz

def test(j):
 P,E,K,z=system(j)
 print('j,old new support',j,len(P[8]),len(P[10]),len(K),len(z),flush=True)
 A=DomainMatrix.from_list_sympy(len(z),len(z[0]),z).to_dense()
 r,d,pv=A.rref_den(method='FF')
 print('rref',j,len(pv),d,flush=True)
 mat=r.to_Matrix();assert len(K) not in pv
 exact={}
 for i,k in enumerate(pv):
  if sum(bool(v)for v in mat[i,:-1])==1 and int(mat[i,-1])%int(d)==0:
   exact[K[k]]=int(mat[i,-1])//int(d)
 print('single exact',len(exact), 'pop unique',len([u for u in exact if u[0]==10]),flush=True)
 for s in P:
  if all((s,*rn(a,b)) in exact for a,b in P[s]):
   print(s,sum(exact[s,*rn(a,b)] for a,b in P[s])+len(E[s]),flush=True)
 json.dump({'positions':P,'exceptions':E,'keys':K,'rows':z,'rref':[[int(y)for y in t]for t in mat.tolist()],'d':str(d),'pivots':pv,'exact':[[*x,c]for x,c in exact.items()]},open(f'harmonic{j}.json','w'))
if __name__=='__main__':
 for j in map(int,sys.argv[1:] or [4,5,6]):test(j)
