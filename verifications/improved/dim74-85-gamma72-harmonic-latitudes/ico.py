from collections import Counter
from fractions import Fraction as F
from itertools import chain, combinations
from math import ceil
from pathlib import Path
import json,sys
import numpy as np
import importlib.util
_parent = Path(__file__).resolve().parent / "model.py"
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
def palette():
 ico=[]
 for a,b in [(a,b)for a in (-1,1)for b in (-1,1)]:
  r=((0,0),(a,0),(0,b))
  ico += [r[i:]+r[:i] for i in range(3)]
 assert len(set(ico))==12 and {dot(x,x) for x in ico}=={(2,1)}
 def cross(x,y):
  return tuple(add(mul(x[(j+1)%3],y[(j+2)%3]),scale(mul(x[(j+2)%3],y[(j+1)%3]),-1))for j in range(3))
 out=[]
 for u in ico:
  for v in ico:
   if dot(u,v)==(0,1):
    out.append(tuple(add(add(scale(a,32),scale(b,12)),scale(c,7))for a,b,c in zip(u,v,cross(u,v))))
 assert len(set(out))==60 and len(set(dot(v,v)for v in out))==1
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
 G=edges(palette());M=60
 D=[M-1+sum((p.P[j]//2)*int(G[j][r].sum())for j in (4,5,6))for r in range(M)]
 T2=[F((M-1)*x)for x in D]
 cert=HERE/'ico60_75.json'
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
 cnt=ceil(alpha);out=dict(dimension=75,degree_counts=dict(Counter(D)),triangles=T,selected=cnt,kissing_lower_bound=p.N8+cnt+12)
 if generate:cert.write_text(json.dumps(dict(out,duals=new),indent=2)+'\n')
 else:
  s=json.loads(cert.read_text())
  assert all(out[k]==s[k] for k in out if k!='degree_counts')
  assert {int(k):v for k,v in s['degree_counts'].items()}==out['degree_counts']
 print(out)
if __name__=='__main__':
 assert __debug__
 main('--generate' in sys.argv)

