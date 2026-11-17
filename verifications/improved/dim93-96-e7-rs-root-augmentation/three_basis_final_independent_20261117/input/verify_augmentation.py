from pathlib import Path
import numpy as np,itertools,json
P=Path(__file__).resolve().parent
b=json.loads((P/'base/certificate.json').read_text());e=json.loads((P/'e7/candidate.json').read_text())
for r in b['compressed_results']:assert np.max(np.abs(r['heads']))<=2
assert np.max(np.abs(b['basis_alphabets']))<=2
EH=np.array([r['head'] for r in e['selected']]);assert np.max(np.abs(EH))<=2 and np.all(EH[:,0]+EH[:,4]==0)
# Orthonormalhead coordinates: sqrt(2)*x0, x1,x2,x3,x5,x6,x7; special coordinate has squared max8<9.
assert 8<9
base={93:42949672960,94:64424509440,95:107374182400,96:206158430208}
def mul(a,b):
 r=0
 while b:
  if b&1:r^=a
  a<<=1
  if a&16:a^=19
  b>>=1
 return r
M=np.array([[mul(a,b) for b in range(16)] for a in range(16)],dtype=np.uint8)
inv=[0]+[next(b for b in range(1,16) if M[a,b]==1) for a in range(1,16)]
for a,b0,c in itertools.product(range(16),repeat=3):
 assert M[M[a,b0],c]==M[a,M[b0,c]] and M[a,b0^c]==(M[a,b0]^M[a,c])
V=np.ones((4,6),dtype=np.uint8)
for j in range(1,4):
 for x in range(6):V[j,x]=M[V[j-1,x],x]
def rank(A):
 A=A.copy();r=0
 for c in range(A.shape[1]):
  z=np.flatnonzero(A[r:,c])
  if not len(z):continue
  k=r+int(z[0]);A[[r,k]]=A[[k,r]];A[r]=M[inv[int(A[r,c])],A[r]]
  for i in range(A.shape[0]):
   if i!=r:A[i]^=M[A[i,c],A[r]]
  r+=1
  if r==A.shape[0]:break
 return r
assert all(rank(V[:,cols])==4 for cols in itertools.combinations(range(6),4))
Signs=np.array([s for s in itertools.product((0,1),repeat=6) if sum(s)%2==0],dtype=np.uint8);DS=np.count_nonzero(Signs[:,None,:]!=Signs[None,:,:],axis=2);assert len(Signs)==32 and DS[DS>0].min()==2
out={}
for n in range(93,97):
 excluded=[88] if n==95 else [];eligible=[i for i in range(n) if i not in excluded];r=len(eligible)-80;assert r in [13,14,16]
 groups=[eligible[:r]]+[eligible[r+16*j:r+16*(j+1)] for j in range(5)]
 assert list(itertools.chain.from_iterable(groups))==eligible and len(set(eligible))==len(eligible)
 co=np.array(list(itertools.product(range(r),range(16),range(16),range(16))),dtype=np.uint8)
 W=np.zeros((len(co),6),dtype=np.uint8)
 for j in range(4):W^=M[co[:,j,None],V[j]]
 assert np.array_equal(W[:,0],co[:,0])
 supports=np.column_stack([np.array(groups[j])[W[:,j]] for j in range(6)])
 assert len(np.unique(supports,axis=0))==r*16**3 and np.all(np.diff(np.sort(supports,axis=1),axis=1)>0)
 edges=[pair for g in groups for pair in itertools.combinations(g,2)]
 if excluded:edges += [(excluded[0],i) for i in eligible]
 assert len(set(tuple(sorted(p)) for p in edges))==len(edges)
 # Everyedge is within one RSgroup, or uses the excludedcoordinate; RS meets it at most once.
 group_of={i:j for j,g in enumerate(groups) for i in g}
 assert all((x in excluded or y in excluded or group_of[x]==group_of[y]) for x,y in edges)
 rs=len(supports)*32;roots=4*len(edges)
 assert 6*16==96 and 3*16==48 and 4*16*2>=96
 assert 6*4*2==48 and 16*48<48**2 and 25*48<48**2
 out[n]=dict(base=base[n],eligible_coordinates=eligible,excluded_coordinates=excluded,groups=groups,rs_count=rs,root_count=roots,total=base[n]+rs+roots)
result=dict(status='PASS finiteRS16 and exactroot/cross geometry; dense proofs inherited from includedbase/e7',field_modulus=19,evaluation_generator=V.tolist(),results=out)
(P/'augmentation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({n:{k:v for k,v in row.items() if k in ['rs_count','root_count','total']} for n,row in out.items()},indent=2))
