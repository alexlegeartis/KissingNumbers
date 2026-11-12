"""Adversarial third review of the Gamma72 69D slice/head theorem.

Rebuilds finite orientation orbits by graph traversal, computes moments from
factorials, compares both coordinate descriptions, and attacks the installed
verifier with mathematical corruptions. No optimization result is trusted.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from math import factorial,prod,gcd,lcm
from collections import deque
import copy,hashlib,importlib.util,json,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PACKAGE=ROOT/'verifications/improved/dim69-gamma72-head-compression'
PATH=PACKAGE/'orthonormal_certificate.json'
C=json.loads(PATH.read_text())
N=6218175600
def require(ok,what):
 if not ok: raise ValueError(what)

# An independent traversal of the generators on individual integer profiles.
# Signed permutations and H/2 each permute the 24 roots and are involutions.
# Root preservation proves preservation of their Z-span (2D4); spanning makes
# the generated action faithful on this finite root set, so the group is finite.
def transform(h,j):
 h=list(h)
 if j<4:h[j]=-h[j];return tuple(h)
 if j<7:
  j-=4;h[j],h[j+1]=h[j+1],h[j];return tuple(h)
 a,b,c,d=h
 z=[a+b+c+d,a+b-c-d,a-b+c-d,a-b-c+d]
 require(all(v%2==0 for v in z),'Hadamard integrality')
 return tuple(v//2 for v in z)

domain=set()
for h in product(range(-8,9),repeat=4):
 if len({v%2 for v in h})!=1:continue
 if all(abs(h[i])+abs(h[j])<=8 for i,j in combinations(range(4),2)):domain.add(h)
normal_roots=set()
for i,j in combinations(range(4),2):
 for si,sj in product((-1,1),repeat=2):
  h=[0]*4;h[i]=8*si;h[j]=8*sj;normal_roots.add(tuple(h))
require(len(domain)==1681 and len(normal_roots)==24,'Complete domains')
for j in range(8):
 require({transform(h,j) for h in domain}==domain,'Profile domain invariant')
 require({transform(h,j) for h in normal_roots}==normal_roots,'Root lattice invariant')
 require(all(transform(transform(h,j),j)==h for h in normal_roots),'Involution')
unseen=set(domain);orbits=[]
while unseen:
 start=min(unseen);todo=deque([start]);orb={start};unseen.remove(start)
 while todo:
  h=todo.popleft()
  for j in range(8):
   t=transform(h,j)
   if t not in orb:orb.add(t);unseen.remove(t);todo.append(t)
 orbits.append(frozenset(orb))
given=[frozenset(map(tuple,g)) for g in C['groups']]
require(set(orbits)==set(given) and len(orbits)==16,'Exact finite averaging orbits')

# Recompute all primitive integer rows from the sphere moment identity.
# We use (2a)!/(2^a a!) rather than the producer's double-factorial product.
rows=[];rhs=[]
for pp in C['moment_powers']:
 require(len(pp)==4 and all(isinstance(a,int) and a>=0 and a%2==0 for a in pp),'Powers')
 degree=sum(pp)//2;require(degree<=5,'Eleven-design applicability')
 gaussian=prod(factorial(a)//(2**(a//2)*factorial(a//2)) for a in pp)
 moment=F(N*128**degree*gaussian,prod(range(72,72+2*degree,2)))
 value=lambda h:prod(v**a for v,a in zip(h,pp))
 b=moment-sum(value(h) for h in normal_roots)
 raw=[sum(value(h) for h in orb) for orb in given]
 scaled=[a*b.denominator for a in raw];g=gcd(*scaled,b.numerator)
 rows.append([a//g for a in scaled]);rhs.append(b.numerator//g)
require(rows==C['rows'] and rhs==C['rhs'],'Independent moment system')
profiles=list(map(tuple,C['profiles']));objective=[len(orb.intersection(profiles)) for orb in given]
require(objective==C['objective'],'Selected-profile population objective')
y=[F(a,C['dual_denominator']) for a in C['dual_multipliers']]
slack=[F(objective[j])-sum(y[i]*rows[i][j] for i in range(len(rows))) for j in range(16)]
require(all(a>=0 for a in slack),'Weak duality')
value=sum(a*b for a,b in zip(y,rhs));bound=-((-value.numerator)//value.denominator)
require(bound==743522304,'Certified bound')

# Independently compare the root-basis and orthonormal descriptions.
other=json.loads((PACKAGE/'certificate.json').read_text())
B=[[1,-1,0,0],[0,1,-1,0],[0,1,1,0],[-1,0,0,1]]
offsets=[]
for h in profiles:
 a=[sum(v*w for v,w in zip(row,h)) for row in B]
 require(all(t%2==0 for t in a),'Integral simple-root profile')
 offsets.append([t//2 for t in a])
require(offsets==other['offsets'] and C['heights']==other['heights'],'Two coordinate descriptions agree')
heights=list(map(F,C['heights']))
inner=lambda a,b:F(sum(x*y for x,y in zip(a,b)),16)
q=[inner(h,h) for h in profiles];alpha=[(8-t*t)/(8-r) for r,t in zip(q,heights)]
require(list(map(str,alpha))==other['alpha_squared'],'Two exact scale lists agree')
def pair_gaps(hh):
 aa=[(8-t*t)/(8-r) for r,t in zip(q,hh)]
 gaps=[]
 for i in range(10):
  for j in range(i+1):
   left=4-inner(profiles[i],profiles[j]);right=4-hh[i]*hh[j]
   require(left>=0 and right>=0,'Signs before squaring')
   gaps.append(right*right-aa[i]*aa[j]*left*left)
 return gaps
gaps=pair_gaps(heights);require(all(g>=0 for g in gaps),'Every distinct-source pair type')
require(all(a*(8-r)+t*t==8 for a,r,t in zip(alpha,q,heights)),'Image norms')

# Exercise the installed verifier. Inputs are copied; published files stay intact.
spec=importlib.util.spec_from_file_location('reviewed_orthonormal69',PACKAGE/'verify_orthonormal.py')
verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
baseline=verifier.check(PATH);require(baseline['bound']==bound,'Installed verifier agrees')
tests=[]
with tempfile.TemporaryDirectory(prefix='gamma72-adversarial-') as td:
 td=Path(td)
 def reject(name,mutate):
  c=copy.deepcopy(C);mutate(c);p=td/(name+'.json');p.write_text(json.dumps(c))
  try:verifier.check(p)
  except (AssertionError,ValueError,RuntimeError,IndexError,KeyError,ZeroDivisionError) as e:
   tests.append({'name':name,'rejected':True,'exception':type(e).__name__});return
  raise ValueError('Corruption incorrectly accepted: '+name)
 reject('false_count',lambda c:c.__setitem__('bound',c['bound']+1))
 reject('wrong_shell_population',lambda c:c.__setitem__('population',N+1))
 reject('duplicate_slice',lambda c:c['profiles'].__setitem__(1,c['profiles'][0]))
 reject('wrong_dual_parity',lambda c:c['profiles'].__setitem__(0,[0,0,0,1]))
 reject('missing_profile_cell',lambda c:c['groups'][0].pop())
 reject('false_moment_rhs',lambda c:c['rhs'].__setitem__(0,c['rhs'][0]+1))
 reject('false_moment_coefficient',lambda c:c['rows'][0].__setitem__(0,c['rows'][0][0]+1))
 reject('false_objective',lambda c:c['objective'].__setitem__(0,1))
 reject('negative_dual_denominator',lambda c:c.__setitem__('dual_denominator',-c['dual_denominator']))
 reject('false_dual_multiplier',lambda c:c['dual_multipliers'].__setitem__(1,c['dual_multipliers'][1]+1))
 reject('excessive_same_slice_height',lambda c:c['heights'].__setitem__(0,'2'))
 # Find a corruption that preserves every within-slice radius requirement
 # but violates a cross-slice pair. This targets the actual gluing inequality.
 cross=None
 for i in range(10):
  for t in [F(a,8) for a in range(-8,9)]:
   if t*t>q[i] or t==heights[i]:continue
   hh=heights.copy();hh[i]=t
   if min(pair_gaps(hh))<0:cross=(i,t);break
  if cross:break
 require(cross is not None,'Nontrivial geometric corruption exists')
 i,t=cross;reject('cross_slice_collision_risk',lambda c:c['heights'].__setitem__(i,str(t)))
 # Refresh the embedding hash in this isolated imported module so rejection
 # reaches the integer Gram calculation rather than stopping at a checksum.
 embedding=json.loads((PACKAGE/'embedding.json').read_text())
 embedding['normals'][0][0]+=1;raw=json.dumps(embedding).encode()
 (td/'embedding.json').write_bytes(raw)
 old_here,old_hash=verifier.HERE,verifier.EMBEDDING_SHA256
 verifier.HERE=td;verifier.EMBEDDING_SHA256=hashlib.sha256(raw).hexdigest()
 try:reject('wrong_normal_embedding_with_refreshed_hash',lambda c:None)
 finally:verifier.HERE,verifier.EMBEDDING_SHA256=old_here,old_hash

report={
 'status':'PASS','dimension':69,'certified_lower_bound':bound,
 'previous_registered_bound':689452922,'gain':bound-689452922,
 'exact_dual_value':str(value),'profile_cells':len(domain),'normal_roots':len(normal_roots),
 'averaging_orbits':len(orbits),'moment_rows':len(rows),'pair_types':len(gaps),
 'minimum_positive_squared_gap':str(min(g for g in gaps if g)),
 'certificate_sha256':hashlib.sha256(PATH.read_bytes()).hexdigest(),
 'root_basis_certificate_sha256':hashlib.sha256((PACKAGE/'certificate.json').read_bytes()).hexdigest(),
 'corruption_tests':tests,
 'mathematical_findings':[
  'Finite normal-coordinate averaging requires no automorphism of the ambient Gamma72 lattice.',
  'Integral products with the embedded 2D4 roots imply exactly the common-parity profile lattice.',
  'The domain covers every non-normal-root minimal vector; all 24 omitted normal roots are subtracted in the moments.',
  'Within a profile the projection map is injective. Between profiles the verified norm and pair bounds rule out image collisions.',
  'All squared geometric inequalities have explicitly checked nonnegative sides.',
  'The integer lower bound follows from a favorable finite orientation and exact rational weak duality, not a numerical LP status.'
 ],
 'inherited_premises':[
  'The supplied catalogue Gram matrix describes the cited integral even 72-dimensional Gamma72 lattice.',
  'Its minimum is 8, its minimal shell has population 6218175600, and that shell is a spherical 11-design.'
 ],
 'scope':'Third exact arithmetic and mathematical review. It does not reprove the published lattice existence/minimum/design theorems or enumerate all source points.'
}
(HERE/'review69_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
