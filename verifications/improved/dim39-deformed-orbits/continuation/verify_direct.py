"""Independent certificate checker; direct orbit products, no Walsh code.
Binary64 BLAS below is exact integer arithmetic: each input/product is an
integer, absolute sum of all 39 terms is <2**53, so every multiply/add/FMA
and every reassociation is exact. No tolerances or normalization floats.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,hashlib,time,sys,itertools
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent
S=Path(sys.argv[1]) if len(sys.argv)>1 else H
if not __debug__:raise RuntimeError('Run without -O.')
t=time.monotonic(); cert=json.loads((S/'orbit_certificate.json').read_text()); raw=(S/'orbit_certificate.json').read_bytes()
V=np.array(cert['representatives'],dtype=np.int64);G=cert['group_words'];lo=cert['norm_lower'];hi=cert['norm_upper']
assert V.shape==(2563,39) and len(G)==len(set(G))==128 and 0 in G
assert all(type(x)is int and 0<=x<2**39 for x in G)
gset=set(G);assert all(a^b in gset for a in G for b in G)
assert all(type(x)is int for row in cert['representatives'] for x in row)
assert abs(V).max()<=10**7 and 39*10**14<2**53
n2=(V*V).sum(1);assert np.all(n2>=lo*lo) and np.all(n2<=hi*hi)
a=np.sort(abs(V),axis=1)[:,::-1]
slacks={'middle':int((2*lo*lo-a[:,:8].sum(1)**2).min()),'root':int((lo*lo-2*a[:,:2].sum(1)**2).min()),'heavy':int((2*lo*lo-(2*a[:,0]+a[:,1:5].sum(1))**2).min())}
assert min(slacks.values())>=0,slacks
# Independently reconstruct sparse layer. Norm^2=32, inner products<=16.
b=json.loads((S/'base_certificate.json').read_text()); supports=[tuple(i for i,c in enumerate(row) if c=='1') for row in (S/'source/supports_3424.txt').read_text().splitlines()]
for r in b['support_replacements']:
 assert supports[r['remove_index']]==tuple(r['remove']);supports[r['remove_index']]=tuple(r['add'])
assert len(supports)==len(set(supports))==3424
seen=set()
for support in supports:
 assert len(support)==len(set(support))==8 and all(0<=i<39 for i in support)
 for subset in itertools.combinations(sorted(support),5):
  assert subset not in seen;seen.add(subset)
del seen
removed=set();stars=[]
for star in b['heavy_stars']:
 c=star['center']; leaves=star['leaves']; assert len(set([c]+leaves))==5 and all(0<=x<39 for x in [c]+leaves)
 removed.update(tuple(sorted((c,i))) for i in leaves)
 for signs in itertools.product([-1,1],repeat=5):
  if np.prod(signs[1:])!=1:continue
  v=[0]*39;v[c]=4*signs[0]
  for j,sgn in zip(leaves,signs[1:]):v[j]=2*sgn
  stars.append(tuple(v))
assert len(stars)==len(set(stars))==64
roots=[p for p in itertools.combinations(range(39),2) if p not in removed]
assert len(roots)*4==2908
for i,v in enumerate(stars):
 assert sum(x*x for x in v)==32
 assert all(sum(x*y for x,y in zip(v,w))<=16 for w in stars[:i])
 assert max(4*(abs(v[a])+abs(v[b])) for a,b in roots)<=16
 for support in supports:
  z=[v[j] for j in support];score=2*sum(abs(x) for x in z)
  if sum(x<0 for x in z)%2 and all(z):score-=4*min(abs(x) for x in z)
  assert score<=16
# Middle-middle: different supports intersect <=4, identical support even
# parity implies >=2 sign changes. Root-middle and root-root <=16 by shape.
sparse=128*len(supports)+4*len(roots)+len(stars);assert sparse==441244
# Direct relative products exhaust all dense-dense interactions because G
# is a group. Only the literal identical representative/group pair is omitted.
signs=np.array([[1-2*((g>>j)&1) for g in G] for j in range(39)],dtype=np.float64)
identity=G.index(0);maximum=-10**18;checked=0;witness=None
for i in range(len(V)):
 D=(V[i:]*V[i]).astype(np.float64)@signs
 D[0,identity]=-10**18
 m=int(D.max())
 if m>maximum:
  maximum=m;j,k=np.unravel_index(D.argmax(),D.shape);witness=[i,i+int(j),G[int(k)]]
 assert 2*m<=lo*lo,(i,m)
 checked+=D.size-1
 if i%500==0:print('rows',i,'maximum',maximum,flush=True)
# Pairwise <=1/2 proves no duplicates, both within and between orbits; sparse
# shapes distinguish sparse families, and dense-sparse <=1/2 excludes overlap.
result={'passed':True,'dimension':39,'total':sparse+len(V)*len(G),'sparse':sparse,'dense':len(V)*len(G),'representatives':len(V),'group_order':len(G),'integer_norm_squared_min':int(n2.min()),'integer_norm_squared_max':int(n2.max()),'sparse_slacks':slacks,'maximum_integer_dot':maximum,'allowed_integer_dot':lo*lo//2,'witness':witness,'complete_relative_products':checked,'certificate_sha256':hashlib.sha256(raw).hexdigest(),'seconds':time.monotonic()-t,'method':'Independent sparse reconstruction and complete direct matrix orbit enumeration. Binary64 operations proven exact integers; no Walsh transform, no tolerances.'}
assert result['total']==cert['claimed_size']==769308
(H/'audit_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
