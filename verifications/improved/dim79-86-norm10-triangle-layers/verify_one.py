from pathlib import Path
import json,math,itertools,hashlib,sys
sys.set_int_max_str_digits(0)
from fractions import Fraction as F
import numpy as np
p=Path(__file__).resolve().parent;dim=int(sys.argv[1]) if len(sys.argv)>1 else 82;suffix=sys.argv[2] if len(sys.argv)>2 else '';raw=json.loads((p/f'exact{dim}{suffix}.json').read_text());blocked=json.loads((p/f'blocked{dim}{suffix}.json').read_text())
N=15281788354560;k=dim-72;M=16*math.comb(k,4);assert raw['dimension']==dim and raw['k']==k and raw['tail_count']==M and F(raw['scale_squared'])==F(8,13);shell={4:24994840000,5:500966912,6:2112250}
tails=[]
for supp in itertools.combinations(range(k),4):
 for signs in itertools.product((-1,1),repeat=4):
  v=[0]*k
  for i,c in zip(supp,signs):v[i]=c
  tails.append(v)
tails=np.array(tails,dtype=np.int64);assert len(tails)==M
r=np.array([1]*4+[0]*(k-4),dtype=np.int64);dotr=tails@r
D=M-1;triangle_twice=F((M-1)*raw['degree']);seen=[]
def oddfac(n):return math.prod(range(n,0,-2))
def wick(r,s,j):
 return sum(math.comb(r,l)*math.comb(s,l)*math.factorial(l)*oddfac(r-l-1)*oddfac(s-l-1)*10**((r+s-2*l)//2)*j**l for l in range(min(r,s)+1) if (r-l)%2==0 and(s-l)%2==0)
for j in (4,5,6):
 points=[]
 for a,b in itertools.product(range(-6,7),repeat=2):
  determinant=1000+2*j*a*b-10*(j*j+a*a+b*b)
  if determinant<0:continue
  ok=True
  for x,y,z in itertools.product(range(-2,3),repeat=3):
   if z==0 or(x==0 and y==0):continue
   q=10*(x*x+y*y+z*z)+2*(j*x*y+a*x*z+b*y*z)
   if 0<q<8:ok=False;break
  if ok:points.append((a,b))
 exceptions=[(10,j),(-10,-j),(j,10),(-j,-10)]
 exps=[(a,b)for a in range(11)for b in range(11-a)if(a+b)%2==0]
 rhs=[]
 for a,b in exps:
  m=(a+b)//2
  total=F(N*10**m*wick(a,b,j),math.prod(72+2*i for i in range(m)))
  rhs.append((total-sum(x**a*y**b for x,y in exceptions))/6**(a+b))
 classes={}
 for t in tails:
  key=(int(np.count_nonzero(t[:4])),int(t[:4].sum()));classes[key]=classes.get(key,0)+1
 for (overlap,dot),multiplicity in classes.items():
  if 8*j+6*dot<=52:continue
  witness=next(x for x in raw['proofs']if(x['j'],x['overlap'],x['dot'])==(j,overlap,dot));assert witness['multiplicity']==multiplicity
  du=list(map(F,witness['dual']))
  marginal=len(du)>len(exps)
  def features(a,b):
   row=[F(a**rr*b**ss,6**(rr+ss))for rr,ss in exps]
   if marginal:row.extend(F(int((a,b)[axis]==val))for axis in range(2)for val in range(-6,7))
   return row
  rhs_here=list(rhs)
  if marginal:
   counts={0:5118989848896,1:7238164864000//2,2:2503639594500//2,3:395496128000//2,4:24994840000//2,5:500966912//2,6:2112250//2}
   rhs_here.extend(F(counts[abs(val)]-sum(pt[axis]==val for pt in exceptions))for axis in range(2)for val in range(-6,7))
  assert len(du)==len(rhs_here)
  if 'row_normalizations' in witness:
   scales=list(map(F,witness['row_normalizations']));assert scales==[F(N)/abs(t)if t else F(1)for t in rhs_here]
   du=[a*b for a,b in zip(du,scales)]
  t=next(t for t in tails if np.count_nonzero(t[:4])==overlap and t[:4].sum()==dot);dott=tails@t
  def cost(a,b):return F(int(np.count_nonzero((np.abs(8*a+6*dotr)>52)&(np.abs(8*b+6*dott)>52))),2)
  margin=min(sum(c*t for c,t in zip(du,features(a,b)))-cost(a,b)for a,b in points);assert margin>=0
  common=sum(c*t for c,t in zip(du,rhs_here))+sum(cost(a,b)for a,b in exceptions)
  assert common==F(witness['common_upper'])
  D+=shell[j]//2*multiplicity;triangle_twice+=shell[j]//2*multiplicity*common
  seen.append((j,overlap,dot));print('dual PASS',seen[-1],len(points),flush=True)
assert len(seen)==len(raw['proofs']);assert D==raw['degree']
T=math.ceil(triangle_twice/2);assert T==F(raw['triangle_upper'])==F(blocked['triangles'])
g=F(0);last=None
for block in blocked['blocks']:
 a,b=block['first_degree'],block['last_degree'];L=b-a+1
 assert last is None or a==last+1
 den=a*a+1-2*T;assert den>a+1 and L>0
 B=(1-F(a+1,den))/b;C=F(2*(a+T),a*den)
 assert B==F(block['B'])and C==F(block['C'])and B>=0 and 0<C<1
 g=max(F(0),1-L*C)*g+L*B/(1+L*C);last=b
assert last==D
vertices=N//2*M;expected=F(vertices,D+1)*(1+g);pairs=math.ceil(expected)
assert expected==F(blocked['expected_pairs']);assert pairs==blocked['pairs']
assert blocked['bound']==6218175600+2*pairs+2*k*(k-1)
# Geometric compatibility: equator |<e,v>|<=5; head-scale squared8/13.
assert F(25*8,13)<16
# Pure-tail D10 points2(+-ei+-ej): squared norm8; cross dot squared<=16*6/13.
assert F(16*6,13)<16
assert F(10*8+4*6,13)==8
out={'status':'EXACT_AUDIT_PASS_WITH_PUBLISHED_SOURCE_PREMISES','dimension':dim,'kissing_bound':blocked['bound'],'pairs':pairs,'degree':D,'triangles':T,'graph_vertices':vertices,'dual_checks':len(seen),'input_sha256':{f:hashlib.sha256((p/f).read_bytes()).hexdigest()for f in(f'exact{dim}{suffix}.json',f'blocked{dim}{suffix}.json')}}
(p/f'audit_result{dim}{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
