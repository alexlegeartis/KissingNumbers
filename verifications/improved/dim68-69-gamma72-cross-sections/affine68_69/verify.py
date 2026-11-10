"""Standard-library verification, independent of both LP and reduction implementations."""
from fractions import Fraction as Q
from itertools import product
from functools import lru_cache
from math import prod,gcd,lcm
from pathlib import Path
import sys,json,time,hashlib
p=Path(__file__).parent;K=int(sys.argv[1]) if len(sys.argv)>1 else 3;C=json.load(open(sys.argv[2] if len(sys.argv)>2 else p/f'certificate{K}.json'));st=time.time();N=6218175600;assert C['population']==N;G=C['gram'];assert len(G)==K and C['dimension']==72-K
embedding_bytes=(p/'embedding.json').read_bytes();assert hashlib.sha256(embedding_bytes).hexdigest()=='15f3643961b517312c136faeb359353f5f202f08158191d8db09cf47a2851e41';E=json.loads(embedding_bytes);Q72=E['catalogue_gram'];U=E['normals'][:K];assert len(Q72)==72 and all(len(x)==72 for x in Q72);assert all(Q72[i][j]==Q72[j][i] for i in range(72) for j in range(72));assert all(isinstance(x,int) for row in U for x in row)
UG=[[sum(U[i][k]*Q72[k][j] for k in range(72)) for j in range(72)] for i in range(K)];assert [[sum(UG[i][k]*U[j][k] for k in range(72)) for j in range(K)] for i in range(K)]==G
# Exact inverse and positivity of the normal Gram.
aug=[[Q(x) for x in row]+[Q(int(i==j)) for j in range(K)] for i,row in enumerate(G)]
for i in range(K):
 assert aug[i][i]>0;v=aug[i][i];aug[i]=[x/v for x in aug[i]]
 for j in range(K):
  if j!=i:
   v=aug[j][i];aug[j]=[a-v*b for a,b in zip(aug[j],aug[i])]
GI=[row[K:] for row in aug]
def proj(a,b):return sum(Q(a[i])*GI[i][j]*b[j] for i in range(K) for j in range(K))
M1={0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800}
for j in range(6):
 lhs=(M1[0] if j==0 else 0)+sum(2*M1[x]*x**(2*j) for x in range(1,5))+2*8**(2*j);rhs=Q(N*64**j*prod(range(1,2*j,2)),prod(72+2*t for t in range(j)));assert lhs==rhs
# Lattice-range necessary constraints and all minima in the generated sublattice within the proof's coefficient box.
W=[];span=[];known={}
for a in product(range(-3,4),repeat=K):
 nu=sum(a[i]*G[i][j]*a[j] for i in range(K) for j in range(K))
 if 0<nu<=48:W.append((a,nu))
 if nu==8 and max(abs(x) for x in a)<=2:
  profile=tuple(sum(a[i]*G[i][j] for i in range(K)) for j in range(K));span.append((a,profile));known[profile]=1
cells=list(map(tuple,C['cells']));assert len(set(cells))==len(cells);idx={c:i for i,c in enumerate(cells)};possible=set();singular=set()
for c in product([0,1,-1,2,-2,3,-3,4,-4,8,-8],repeat=K):
 nn=proj(c,c)
 if nn>8:continue
 if any(2*abs(sum(x*y for x,y in zip(a,c)))>nu and not(nu==8 and abs(sum(x*y for x,y in zip(a,c)))==8) for a,nu in W):continue
 if c in known:continue
 possible.add(c)
 if nn==8:singular.add(c)
assert set(cells)==possible
# Independently eliminate true cell identifications, not approximate row dependencies.
parent=list(range(len(cells)));fixed=[]
def root(x):
 while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
 return x
def unite(a,b):
 a=root(a);b=root(b)
 if a!=b:parent[max(a,b)]=min(a,b)
for j,c in enumerate(cells):unite(j,idx[tuple(-x for x in c)])
for a,profile in span:
 for j,c in enumerate(cells):
  if sum(x*y for x,y in zip(a,c))!=4:continue
  cp=tuple(x-y for x,y in zip(profile,c))
  if cp in idx:unite(j,idx[cp])
  else:fixed.append((j,Q(known.get(cp,0))))
fix={}
for j,value in fixed:
 r=root(j)
 if r in fix:assert fix[r]==value
 fix[r]=value
roots=sorted(set(root(j) for j in range(len(cells)))-set(fix));groups=[[j for j in range(len(cells)) if root(j)==r] for r in roots];fixedcells={str(j):str(fix[root(j)]) for j in range(len(cells)) if root(j) in fix};assert groups==C['groups'] and fixedcells==C['fixed_cells'];upper=[1 if any(cells[j] in singular for j in group) else None for group in groups];assert upper==C['upper']
@lru_cache(None)
def gaussian(ps):
 if sum(ps)==0:return Q(1)
 if sum(ps)%2:return Q(0)
 i=next(i for i,x in enumerate(ps) if x);q=list(ps);q[i]-=1;total=Q(0)
 for j,count in enumerate(q):
  if count:
   r=q.copy();r[j]-=1;total+=count*G[i][j]*gaussian(tuple(r))
 return total
def primitive(row,rhs):
 den=lcm(*[x.denominator for x in row],rhs.denominator);v=[int(x*den) for x in row];b=int(rhs*den);d=0
 for x in v+[b]:d=gcd(d,abs(x))
 if d:v=[x//d for x in v];b//=d
 return [(i,x) for i,x in enumerate(v) if x],b
D=int(C['denominator']);coef=list(map(int,C['multipliers']));offsets=list(map(tuple,C['offsets']));assert len(offsets)==len(set(offsets));objective=[sum(cells[j] in offsets for j in group) for group in groups];constant=sum(Q(value) for j,value in fixedcells.items() if cells[int(j)] in offsets);assert objective==C['objective'] and str(constant)==C['objective_constant'];slack=[D*x for x in objective];numerator=constant*D
assert len(coef)==len(C['tags'])
for nr,(a,tag) in enumerate(zip(coef,C['tags'])):
 kind=tag[0]
 if kind=='mom':
  ps=tag[1];assert len(ps)==K and min(ps)>=0 and sum(ps)<=10 and sum(ps)%2==0;degree=sum(ps)
  def val(c):return prod(Q(x,8)**power for x,power in zip(c,ps))
  full=[val(c) for c in cells];rhs=Q(N)*gaussian(tuple(ps))/Q(8**(degree//2)*prod(72+2*j for j in range(degree//2)))-sum(val(c) for c in known)
 elif kind=='marg':
  i,x=tag[1:];assert 0<=i<K and abs(x)<=4;full=[Q(int(c[i]==x)) for c in cells];rhs=Q(M1[abs(x)]-sum(c[i]==x for c in known))
 else:raise AssertionError(kind)
 rhs-=sum(full[int(j)]*Q(value) for j,value in fixedcells.items());row=[sum(full[j] for j in group) for group in groups];ir,ib=primitive(row,rhs)
 if 'rows' in C:assert [list(x) for x in ir]==C['rows'][nr] and ib==C['rhs'][nr]
 numerator+=a*ib
 for j,v in ir:slack[j]-=a*v
assert all(x>=0 for x,u in zip(slack,upper) if u is None);numerator-=sum(u*max(0,-s) for s,u in zip(slack,upper) if u is not None);assert numerator==int(C['numerator']);bound=-((-numerator)//D);assert bound==C['bound'];geometry=[]
for i,x in enumerate(offsets):
 for j,y in enumerate(offsets[:i+1]):
  nx=8-proj(x,x);ny=8-proj(y,y);P=4-proj(x,y);assert nx>0 and ny>0;assert P<=0 or 4*P*P<=nx*ny;geometry.append(str(nx*ny-4*P*P))
baseline={3:627822180,4:375222168}[K];report=dict(status='PASS',dimension=72-K,bound=int(bound),baseline=baseline,gain=int(bound)-baseline,cells=len(cells),groups=len(groups),constraints=len(coef),slice_count=len(offsets),geometric_pairs=len(geometry),minimum_squared_margin=str(min(map(Q,geometry))),embedding_checked=True,assumptions=['Published Gamma72 catalogue Gram, minimum 8, population 6,218,175,600 and 11-design property.'],scope='Exact affine-slice existence lower bound; no full shell enumeration or worldwide-priority claim.')
if len(sys.argv)<3:(p/f'independent{K}.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
