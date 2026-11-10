"""Second independent audit: polynomial moments, graph components, integer dual."""
import itertools,json,math,sys,re,hashlib,subprocess
from fractions import Fraction as F
from pathlib import Path
import sympy as s
p=Path(__file__).parent;K=int(sys.argv[1]);C=json.loads(Path(sys.argv[2] if len(sys.argv)>2 else p/f'certificate{K}.json').read_text());N=6218175600
assert __debug__;assert C['population']==N and C['dimension']==72-K
G=s.Matrix(C['gram']);H=G.inv();cells=list(map(tuple,C['cells']));idx={c:j for j,c in enumerate(cells)}
assert len(cells)==len(idx);assert all(G[:i,:i].det()>0 for i in range(1,K+1))
embedding_bytes=(p/'embedding.json').read_bytes();assert hashlib.sha256(embedding_bytes).hexdigest()=='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842';E=json.loads(embedding_bytes);B=s.Matrix(E['catalogue_gram']);U=s.Matrix(E['normals'][:K]);assert U*B*U.T==G
# Independent polynomial coefficient moments; no Wick/matching routine.
x=s.symbols('x:'+str(K));quad=sum(G[i,j]*x[i]*x[j] for i in range(K) for j in range(K))
polys=[s.Poly(quad**j,*x) for j in range(6)]
roots=[];known=[]
for a in itertools.product(range(-2,3),repeat=K):
 if sum(a[i]*G[i,j]*a[j] for i in range(K) for j in range(K))==8:
  roots.append(a);known.append(tuple(G*s.Matrix(a)))
# A simpler sufficient domain: constraints from all known norm-8 roots.
possible=[]
for c in itertools.product(range(-4,5),repeat=K):
 if c in known:continue
 if all(abs(sum(a[i]*c[i] for i in range(K)))<=4 for a in roots):possible.append(c)
assert set(possible)==set(cells),('Root-only domain differs',len(possible),len(cells))
# Every profile with coordinate +/-8 is +/- one chosen normal, hence known.
# Averaging generators are integral isometries of the normal Gram.
# Each transformed normal is an actual integer combination of the same lattice normals.
# Positive-definite integral isometries form a finite group. Average the actual
# cell populations over this group; at least one orientation attains the average.
generators=[]
for k in range(K):
 T=s.eye(K)
 for i in range(K):T[i,k]-=G[i,k]/4
 assert all(v.q==1 for v in T) and T*G*T.T==G
 generators.append(T)
if K==4 and G[1,2]==G[1,3]==G[2,3]==0:
 for a,b in [(1,2),(2,3)]:
  T=s.eye(K);T.row_swap(a,b);assert T*G*T.T==G;generators.append(T)
for T in generators:
 assert {tuple(T*s.Matrix(c)) for c in known}==set(known)
# Reconstruct identifications as an undirected graph and find its components.
adj=[set() for c in cells]
for j,c in enumerate(cells):
 adj[j].add(idx[tuple(-v for v in c)])
 for a,v in zip(roots,known):
  if sum(ai*ci for ai,ci in zip(a,c))==4:
   other=tuple(vi-ci for vi,ci in zip(v,c));assert other in idx
   adj[j].add(idx[other])
for T in generators:
 TT=[[int(T[i,j]) for j in range(K)] for i in range(K)]
 for j,c in enumerate(cells):
  other=tuple(sum(TT[i][k]*c[k] for k in range(K)) for i in range(K));assert other in idx;adj[j].add(idx[other])
remaining=set(range(len(cells)));components=[]
while remaining:
 todo=[min(remaining)];comp=set()
 while todo:
  j=todo.pop()
  if j in comp:continue
  comp.add(j);todo.extend(adj[j]-comp)
 remaining-=comp;components.append(sorted(comp))
assert components==C['groups'];assert C.get('fixed_cells',{})=={};assert C.get('objective_constant','0')=='0'
assert all(v is None for v in C.get('upper',[None]*len(components)))
# Known one-point counts independently satisfy all six even moments.
M={0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800}
for k in range(6):
 assert (M[0] if k==0 else 0)+2*sum(M[i]*i**(2*k) for i in range(1,5))+2*8**(2*k)==F(N*64**k*math.prod(range(1,2*k,2)),math.prod(72+2*j for j in range(k)))
assert len(C['rows'])==len(C['rhs'])==len(C['tags'])==len(C['multipliers'])
for row,b,tag in zip(C['rows'],C['rhs'],C['tags']):
 if tag[0]=='mom':
  powers=tag[1];d=sum(powers);k=d//2;assert d<=10 and d%2==0
  val=lambda c:math.prod(v**a for v,a in zip(c,powers))
  coeff=int(polys[k].coeff_monomial(math.prod(z**a for z,a in zip(x,powers))))
  multinomial=math.factorial(d)//math.prod(math.factorial(a) for a in powers)
  rhs=F(N*8**k*math.prod(range(1,2*k,2))*coeff,math.prod(72+2*j for j in range(k))*multinomial)-sum(val(c) for c in known)
 elif tag[0]=='marg':
  i,v=tag[1:];val=lambda c:int(c[i]==v);rhs=F(M[abs(v)]-sum(val(c) for c in known))
 else:raise AssertionError(tag)
 wanted=[sum(val(cells[j]) for j in group) for group in components];actual=dict(row)
 j=next(j for j,v in enumerate(wanted) if v);f=F(actual[j],wanted[j]);assert f
 assert all(actual.get(j,0)==f*v for j,v in enumerate(wanted));assert b==f*rhs
D=int(C['denominator']);a=list(map(int,C['multipliers']));offs=list(map(tuple,C['offsets']));assert len(set(offs))==len(offs) and all(c in idx for c in offs)
obj=[sum(cells[j] in offs for j in comp) for comp in components];assert obj==C['objective'];sl=[D*v for v in obj]
for aa,row in zip(a,C['rows']):
 for j,v in row:sl[j]-=aa*v
assert min(sl)>=0;num=sum(aa*b for aa,b in zip(a,C['rhs']));assert num==int(C['numerator']);bound=-((-num)//D);assert bound==C['bound']
projection=lambda c,d:(s.Matrix(c).T*H*s.Matrix(d))[0]
radii=[8-projection(c,c) for c in offs];assert min(radii)>0
gaps=[]
for i,c in enumerate(offs):
 for j,d in enumerate(offs[:i+1]):
  inner=4-projection(c,d);gap=radii[i]*radii[j]-4*inner**2
  assert inner<=0 or gap>=0;gaps.append(gap)
R=dict(dimension=72-K,bound=bound,exact_value=str(F(num,D)),cells=len(cells),groups=len(components),rows=len(C['rows']),known_minima=len(known),slack_min=min(sl),slack_max=max(sl),radii=list(map(str,radii)),minimum_geometry_gap=str(min(gaps)),gain_over_supplied_benchmark=bound-{3:689442146,4:471366656,5:355090822}[K],status='PASS')
subprocess.run([sys.executable, '-B', str(p/'verify_orientations.py'), str(K), str(Path(sys.argv[2]) if len(sys.argv)>2 else p/f'certificate{K}.json')],check=True,capture_output=True,text=True)
(p/f'verification{K}.json').write_text(json.dumps(R,indent=2));print(json.dumps(R,indent=2))
