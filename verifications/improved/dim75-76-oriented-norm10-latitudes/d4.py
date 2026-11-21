from collections import Counter
from fractions import Fraction as F
from itertools import chain, combinations
from math import ceil
from pathlib import Path
import json,sys
import numpy as np
import importlib.util
_parent = Path(__file__).resolve().parent.parent / "dim76-85-mixed-shape-norm10" / "verify.py"
_spec = importlib.util.spec_from_file_location("norm10_parent", _parent)
p = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p)
if '--generate' in sys.argv:
 from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent

# Arithmetic in Q[phi], phi^2=phi+1. No trigonometric or rounded comparisons.
def add(p,q):return p[0]+q[0],p[1]+q[1]
def mul(p,q):return p[0]*q[0]+p[1]*q[1],p[0]*q[1]+p[1]*q[0]+p[1]*q[1]
def scale(p,q):return p[0]*q,p[1]*q
def sign(p):
 a,b=2*p[0]+p[1],p[1]
 if a>=0 and b>=0:return int(bool(a or b))
 if a<=0 and b<=0:return -int(bool(a or b))
 h=a*a-5*b*b
 assert h!=0
 return (1 if a>0 else -1)*(1 if h>0 else -1)
def dot(x,y):
 w=(0,0)
 for i,j in zip(x,y):w=add(w,mul(i,j))
 return w
def palette600():
 from itertools import permutations,product
 out=[]
 for i in range(4):
  for s in (-2,2):
   out.append(tuple((s if i==j else 0,0)for j in range(4)))
 out.extend(tuple((x,0) for x in signs)for signs in product((-1,1),repeat=4))
 for p in permutations(range(4)):
  if sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))%2:continue
  for a,b,c in product((-1,1),repeat=3):
   t=((0,0),(a,0),(0,b),(-c,c))
   out.append(tuple(t[j]for j in p))
 assert len(out)==len(set(out))==120 and {dot(v,v)for v in out}=={(4,0)}
 return out

def paletteD4():
 from itertools import product
 out=[tuple((t,0)for t in v)for v in product(range(-3,4),repeat=4)if sum(t*t for t in v)==10]
 assert len(out)==144
 return out

def paletteD4q8():
 from itertools import product
 out=[tuple((t,0)for t in v)for v in product(range(-2,3),repeat=4)if sum(t*t for t in v)==8]
 assert len(out)==24
 return out

def edges(pts):
 M=len(pts); g={}
 for j in [0,1,2,3,4,5,6,10]:
  mat=np.zeros((M,M),bool)
  c=F(25-4*j,10)
  for i,x in enumerate(pts):
   for k,y in enumerate(pts):
    if j<4: continue
    if j==10:mat[i,k]=True;continue
    u=dot(x,y)
    # t = dot(x,y) / sqrt(norm(x)^2 norm(y)^2); c is positive.
    assert c>0
    v=add(mul(u,u),scale(mul(dot(x,x),dot(y,y)),-c*c))
    mat[i,k]=sign(u)>0 and sign(v)>0
  g[j]=mat
 for j in [4,5,6,10]:assert (g[j]==g[j].T).all()
 # Oriented—not antipodally identified—shell vertices. Negative head products
 # cannot conflict, since even head product 3 cannot contribute an edge.
 for j in [1,2,3,4,5,6,10]:g[-j]=g[0]
 assert (g[4]==np.eye(M,dtype=bool)).all()
 return g
def main(generate=False):
 p.independently_replay_parent_counts()
 assert F(16,25)*25==16 and F(8,5)*8<16
 G=edges(paletteD4());M=len(G[0])
 D=[M-1+sum((p.P[j]//2)*int(G[j][r].sum())for j in (4,5,6))for r in range(M)]
 T2=[F((M-1)*x)for x in D]
 cert=HERE/'d4_76.json'
 prev={} if generate else json.loads(cert.read_text())['duals']
 new={}
 for j in (4,5,6):
  pts,ex,rhs,features=p.profile_data(j)
  sc=[F(p.N)/abs(v) if v else F(1) for v in rhs]
  if generate:
   af=np.array([[float(t*s)for t,s in zip(row,sc)]for row in features]).T
   rf=[float(t*s/p.N)for t,s in zip(rhs,sc)]
  cache={}
  for a in range(M):
   for b in np.flatnonzero(G[j][a]):
    def obj(x,y):return int(np.sum(G[x][a]&G[y][b]))
    o=tuple(obj(x,y) for x,y in pts)
    if o not in cache:
     label=f'{j}:{len(cache)}'
     costs=[F(v)for v in o]
     if generate:
      sol=linprog([-float(x)for x in costs],A_eq=af,b_eq=rf,method='highs')
      assert sol.success,sol.message
      y=[F(float(-u)).limit_denominator(10**12)*s for u,s in zip(sol.eqlin.marginals,sc)]
      y[0]+=max(F(0),max(t-sum(a*u for a,u in zip(row,y))for row,t in zip(features,costs)))
      new[label]=[str(v)for v in y]
     else:y=list(map(F,prev[label]))
     assert len(y)==len(rhs)
     assert all(sum(a*u for a,u in zip(row,y))>=t for row,t in zip(features,costs))
     cache[o]=sum(y_i*x for y_i,x in zip(y,rhs))
    C=cache[o]+sum(F(obj(x,y))for x,y in ex)
    T2[a]+=C*(p.P[j]//2)
 T=ceil(max(T2)/2)
 alpha=sum(F(p.N,d+1)*p.low_factor(d,T) for d in D)
 cnt=ceil(alpha);out=dict(dimension=76,degree_counts=dict(Counter(D)),triangles=T,selected=cnt,kissing_lower_bound=p.N8+cnt+24)
 if generate:cert.write_text(json.dumps(dict(out,duals=new),indent=2)+'\n')
 else:
  s=json.loads(cert.read_text())
  assert all(out[k]==s[k] for k in out if k!='degree_counts')
  assert {int(k):v for k,v in s['degree_counts'].items()}==out['degree_counts']
 print(out)
if __name__=='__main__':
 assert __debug__
 main('--generate' in sys.argv)
