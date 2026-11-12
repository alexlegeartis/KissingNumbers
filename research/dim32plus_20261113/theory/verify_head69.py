"""Exact 69D head-compression proof using only Python's standard library.

Reconstructs the profile domain, all orbit identifications, the Gram embedding,
Gaussian moment recurrences, the rational dual, and every geometric pair type.
No optimizer or numerical array package is imported.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from math import prod,gcd,lcm
import json,hashlib,sys

p=Path(__file__).parent;root=p.resolve().parents[2]
C=json.loads(Path(sys.argv[1] if len(sys.argv)>1 else p/'head69_certificate.json').read_text())
G=C['gram'];k=4;N=6218175600
assert C['dimension']==69 and C['population']==N and C['norm_squared']==8
assert G==[[8,-4,-4,-4],[-4,8,0,0],[-4,0,8,0],[-4,0,0,8]]
def transpose(A):return [list(x) for x in zip(*A)]
def matmul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def bilinear(a,B,b):return sum(a[i]*B[i][j]*b[j] for i in range(len(a)) for j in range(len(b)))
def inverse(A):
 n=len(A);B=[[F(A[i][j]) for j in range(n)]+[F(i==j) for j in range(n)] for i in range(n)]
 for i in range(n):
  pivot=next(j for j in range(i,n) if B[j][i]);B[i],B[pivot]=B[pivot],B[i]
  z=B[i][i];B[i]=[v/z for v in B[i]]
  for j in range(n):
   if i!=j:
    z=B[j][i];B[j]=[v-z*w for v,w in zip(B[j],B[i])]
 return [row[n:] for row in B]
H=inverse(G)
embedding=root/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json'
raw=embedding.read_bytes();assert hashlib.sha256(raw).hexdigest()=='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842'
E=json.loads(raw);U=E['normals'][:k]
assert matmul(matmul(U,E['catalogue_gram']),transpose(U))==G
cells=list(map(tuple,C['cells']));idx={c:i for i,c in enumerate(cells)}
assert len(cells)==len(idx)
rootcoeff=[a for a in product(range(-2,3),repeat=k) if bilinear(a,G,a)==8]
known=[tuple(sum(a[i]*G[i][j] for i in range(k)) for j in range(k)) for a in rootcoeff]
assert len(known)==24 and len(set(known))==24
possible=[c for c in product(range(-4,5),repeat=k) if c not in known and all(abs(sum(a[i]*c[i] for i in range(k)))<=4 for a in rootcoeff)]
assert set(possible)==set(cells)
generators=[]
for j in range(k):
 T=[[int(i==l) for l in range(k)] for i in range(k)]
 for i in range(k):
  assert G[i][j]%4==0;T[i][j]-=G[i][j]//4
 generators.append(T)
for a,b in ((1,2),(2,3)):
 T=[[int(i==l) for l in range(k)] for i in range(k)];T[a],T[b]=T[b],T[a];generators.append(T)
act=lambda T,c:tuple(sum(row[j]*c[j] for j in range(k)) for row in T)
for T in generators:
 assert matmul(matmul(T,G),transpose(T))==G
 assert {act(T,c) for c in known}==set(known)
 assert {act(T,c) for c in cells}==set(cells)
parent=list(range(len(cells)))
def find(a):
 while a!=parent[a]:parent[a]=parent[parent[a]];a=parent[a]
 return a
def union(a,b):
 a=find(a);b=find(b)
 if a!=b:parent[max(a,b)]=min(a,b)
for i,c in enumerate(cells):
 union(i,idx[tuple(-v for v in c)])
 for a,v in zip(rootcoeff,known):
  if sum(x*y for x,y in zip(a,c))==4:
   d=tuple(x-y for x,y in zip(v,c));assert d in idx;union(i,idx[d])
 for T in generators:union(i,idx[act(T,c)])
parts={}
for i in range(len(cells)):parts.setdefault(find(i),[]).append(i)
groups=sorted(parts.values(),key=lambda x:x[0]);assert groups==C['groups']
@lru_cache(None)
def moment(exponents):
 if not sum(exponents):return 1
 i=next(i for i,a in enumerate(exponents) if a)
 b=list(exponents);b[i]-=1;answer=0
 for j,m in enumerate(b):
  if m:
   b[j]-=1;answer+=m*G[i][j]*moment(tuple(b));b[j]+=1
 return answer
def primitive(row,b):
 D=lcm(b.denominator,*[F(v).denominator for v in row]);v=[int(F(x)*D) for x in row];b=int(b*D)
 g=gcd(*v,b) or 1
 return [x//g for x in v],b//g
for sparse,b,tag in zip(C['rows'],C['rhs'],C['tags']):
 assert tag[0]=='mom';powers=tuple(tag[1]);d=sum(powers);assert d%2==0 and d<=10
 val=lambda c:prod(x**a for x,a in zip(c,powers))
 rhs=F(N*8**(d//2)*moment(powers),prod(72+2*j for j in range(d//2)))-sum(val(c) for c in known)
 row=[sum(val(cells[j]) for j in group) for group in groups]
 row,rhs=primitive(row,rhs)
 assert sparse==[[j,v] for j,v in enumerate(row) if v] and b==rhs
offsets=list(map(tuple,C['offsets']));assert len(offsets)==len(set(offsets))==10 and all(c in idx for c in offsets)
objective=[sum(cells[j] in offsets for j in group) for group in groups];assert objective==C['objective']
y=list(map(F,C['multipliers']));assert len(y)==len(C['rows'])==len(C['rhs'])==len(C['tags'])
slack=list(map(F,objective))
for yi,row in zip(y,C['rows']):
 for j,v in row:slack[j]-=yi*v
assert min(slack)>=0 and list(map(str,slack))==C['slack']
value=sum(a*b for a,b in zip(y,C['rhs']));assert str(value)==C['exact_lower_bound']
bound=-((-value.numerator)//value.denominator);assert bound==C['bound']==743522304
heights=list(map(F,C['heights']));alpha2=list(map(F,C['alpha_squared']));assert len(heights)==len(alpha2)==10
for c,h,a in zip(offsets,heights,alpha2):
 radius=8-bilinear(c,H,c);assert radius>0 and a>0 and a*radius+h*h==8
gaps=[]
for i,c in enumerate(offsets):
 for j,d in enumerate(offsets[:i+1]):
  left=4-bilinear(c,H,d);right=4-heights[i]*heights[j]
  gap=right*right-alpha2[i]*alpha2[j]*left*left
  assert right>=0 and (left<=0 or gap>=0),(i,j,gap)
  gaps.append({'i':i,'j':j,'squared_gap':str(gap)})
assert gaps==C['geometry_gaps']
report={'status':'PASS','dimension':69,'bound':bound,'previous':689452922,'gain':bound-689452922,'geometry_pair_types':len(gaps),'profile_cells':len(cells),'averaged_variables':len(groups),'exact_moment_rows':len(y),'dual_value':str(value),'dual_slack_min':str(min(slack)),'inherited_inputs':['Nebe Gamma72 existence, integral even lattice, minimum8, minimal-shell population6218175600, spherical11-design, and identity of the supplied catalogue Gram matrix.'],'scope':'Existence of one favorable orientation of the D4 normals, with explicit algebraic slice/head maps; not an enumeration of all743522304points.'}
print(json.dumps(report,indent=2))
