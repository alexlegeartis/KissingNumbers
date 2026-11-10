"""Independent exact verification. No optimizer and no imports from source LP code."""
from fractions import Fraction as Q
from itertools import product,combinations
from math import gcd,lcm,prod
from functools import lru_cache
from pathlib import Path
import json,sys
p=Path(__file__).parent;C=json.load(open(sys.argv[1] if len(sys.argv)>1 else p/'certificate.json'));N=C['population'];G=C['gram'];assert N==6218175600 and G==[[8,4],[4,8]];cells=[tuple(x) for x in C['cells']];idx={c:i for i,c in enumerate(cells)};M1={int(k):v for k,v in C['onepoint'].items()};mu=8;n=72
# Independent one-point design equations. First five form a nonsingular Vandermonde system in0,1,4,9,16.
for k in range(6):
 lhs=(M1[0] if k==0 else 0)+sum(2*M1[c]*c**(2*k) for c in range(1,5))+2*8**(2*k)
 rhs=Q(N*64**k*prod(range(1,2*k,2)),prod(n+2*j for j in range(k)))
 assert lhs==rhs
assert M1[4]>0 # Guarantees existence of two minimal vectors with inner product4.
# Enumerate all potentially occupied joint cells, using only PSD and lattice minimum.
known={};span=[];W=[]
for a in product(range(-3,4),repeat=2):
 nu=sum(a[i]*G[i][j]*a[j] for i in range(2) for j in range(2))
 if 0<nu<=48:W.append((a,nu))
 if nu==8:
  c=tuple(sum(a[i]*G[i][j] for i in range(2)) for j in range(2));known[c]=1;span.append((a,c))
def projected(c,d):return Q(2*c[0]*d[0]-c[0]*d[1]-c[1]*d[0]+2*c[1]*d[1],12)
possible=[];uppers={}
for c in product([0,1,-1,2,-2,3,-3,4,-4,8,-8],repeat=2):
 nn=projected(c,c)
 if nn>8:continue
 admissible=True
 for a,nu in W:
  ip=sum(x*y for x,y in zip(a,c))
  if 2*abs(ip)>nu and not(nu==8 and abs(ip)==8):admissible=False;break
 if not admissible:continue
 if c in known:continue
 possible.append(c);uppers[c]=1 if nn==8 else None
assert set(possible)==set(cells) and len(possible)==len(cells)
assert [uppers[c] for c in cells]==C['upper']
@lru_cache(None)
def gaussian(a,b):
 if a+b==0:return Q(1)
 if a+b<0 or a<0 or b<0 or (a+b)%2:return Q(0)
 if a:return (a-1)*G[0][0]*gaussian(a-2,b)+b*G[0][1]*gaussian(a-1,b-1)
 return (b-1)*G[1][1]*gaussian(0,b-2)
def primitive(row,rhs):
 den=lcm(*[x.denominator for x in row],rhs.denominator);v=[int(x*den) for x in row];b=int(rhs*den);g=0
 for x in v+[b]:g=gcd(g,abs(x))
 if g:v=[x//g for x in v];b//=g
 return [[i,x] for i,x in enumerate(v) if x],b
for no,tag in enumerate(C['tags']):
 kind=tag[0];row=[Q(0)]*len(cells);rhs=Q(0)
 if kind=='mom':
  a,b=tag[1];deg=a+b;assert deg%2==0 and deg<=10
  def value(c):return Q(c[0],8)**a*Q(c[1],8)**b
  row=[value(c) for c in cells];sphere=gaussian(a,b)/Q(8**(deg//2)*prod(n+2*j for j in range(deg//2)));rhs=N*sphere-sum(value(c) for c in known)
 elif kind=='marg':
  i,x=tag[1:];row=[Q(int(c[i]==x)) for c in cells];rhs=Q(M1[abs(x)]-sum(c[i]==x for c in known))
 elif kind=='sym':
  c=tuple(tag[1]);row[idx[c]]=1;row[idx[tuple(-x for x in c)]]=-1
 elif kind=='tr':
  a,c=tuple(tag[1]),tuple(tag[2]);assert sum(a[i]*G[i][j]*a[j] for i in range(2) for j in range(2))==8;assert sum(x*y for x,y in zip(a,c))==4
  profile=tuple(sum(a[i]*G[i][j] for i in range(2)) for j in range(2));cp=tuple(x-y for x,y in zip(profile,c));row[idx[c]]=1
  if cp in idx:row[idx[cp]]=-1
  elif cp in known:rhs=Q(1)
 else:raise AssertionError(kind)
 row=[Q(x) for x in row];rr,bb=primitive(row,rhs);assert rr==C['rows'][no] and bb==C['rhs'][no],(no,tag)
# Replay the supplied exact dual solely in integer arithmetic.
D=int(C['denominator']);a=list(map(int,C['multipliers']));objective=[int(tuple(c) in map(tuple,C['offsets'])) for c in cells];assert objective==C['objective'];slack=[D*x for x in objective]
for aa,row in zip(a,C['rows']):
 for j,v in row:slack[j]-=aa*v
assert all(s>=0 for s,u in zip(slack,C['upper']) if u is None)
num=sum(x*y for x,y in zip(a,C['rhs']))-sum(u*max(0,-s) for s,u in zip(slack,C['upper']) if u is not None);assert num==int(C['numerator']);bound=-((-num)//D);assert bound==C['bound']==1360986112
# Exact projected geometric compatibility, including within each slice.
offsets=C['offsets'];checks=[]
for i in range(len(offsets)):
 for j in range(i+1):
  x,y=offsets[i],offsets[j];xx=projected(x,x);yy=projected(y,y);xy=projected(x,y);assert xx<8 and yy<8;P=4-xy;assert P<=0 or 4*P*P<=(8-xx)*(8-yy);checks.append(dict(i=i,j=j,squared_margin=str((8-xx)*(8-yy)-4*P*P)))
report=dict(status='PASS',dimension=70,bound=bound,baseline=1249778250,gain=bound-1249778250,cells=len(cells),constraints=len(C['tags']),geometry_checks=checks,assumptions=['Gamma72 is an even unimodular 72D lattice of minimum 8 with 6,218,175,600 minimal vectors.','Its minimal shell is an 11-design.'],scope='Exact existence lower bound from three affine slices; no full shell enumeration and no global priority claim.')

if len(sys.argv)==1:(p/'independent_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
