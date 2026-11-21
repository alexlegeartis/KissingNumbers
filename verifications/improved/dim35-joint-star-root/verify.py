"""Independent exact replay of a 415,360-point kissing code in R^35."""
import hashlib,itertools as it,json
from pathlib import Path
import numpy as np
if not __debug__:
    raise RuntimeError("Run without -O")
P=Path(__file__).resolve().parent
shared=P.parent/'dim35-shared-root-heads'
switched=P.parent/'dim35-octad-switches'
src={'code_data.json':shared/'code_data.json','witness.json':switched/'witness.json',
     'shared_positive.npz':shared/'shared_positive.npz',
     'replacement.npz':P/'replacement.npz'}
N=35
for f,h in {
'code_data.json':'0f44ca2ff564941d58f92cbb0c47ada139612196fd3d69a547d5318f78009d05',
'witness.json':'3e8a212b8a8ee4c6d81a45134dae67b99fbfbc4a4f990fc174c42af987e01e16',
'shared_positive.npz':'9d4b35ac5d18f40011282ec2e83753bdffee947351584a77a910baef875768cf',
'replacement.npz':'ffb73e3082a0ed5a1f94ce8e5c32febb0683e151bcd1aee3838eb52a6e76320c'}.items():
 assert hashlib.sha256(src[f].read_bytes()).hexdigest()==h
C=json.loads(src['code_data.json'].read_text()); W=json.loads(src['witness.json'].read_text())
F=np.load(src['shared_positive.npz'],allow_pickle=False)
Z=np.load(P/'replacement.npz',allow_pickle=False)
S,R,H=[Z[k] for k in ('stars','roots','heads')]; D=F['dense_positions']
assert np.array_equal(H,F['heads'])
assert S.shape==(262,N) and R.shape==(2158,N) and H.shape==(12,N)
assert all(z.dtype.kind in 'iu' for z in (S,R,H))
for z,shape in zip((S,R,H),([0]*30+[2]*4+[4],[0]*33+[4,4],[1]*33+[3,3])):
 assert np.array_equal(np.sort(abs(z),axis=1),np.tile(shape,(len(z),1)))
 assert len(z)==len(set(map(tuple,z)))
assert len(D)==32 and len(set(D.tolist()))==32 and min(D)>=0 and max(D)<35
v=[int(m) for m,_ in C['full_parity_supports']]
def mask(s):return sum(1<<j for j in s)
seen=set()
for i,b,a in W['octads']:
 assert 0<=i<len(v) and i not in seen and v[i]==mask(b)
 seen.add(i);v[i]=mask(a)
assert len(v)==len(set(v))==2202
assert all(x.bit_count()==8 and 0<=x<1<<N for x in v)
assert all((x&y).bit_count()<=4 for x,y in it.combinations(v,2))
G=np.r_[S,R,H].astype(np.int64); norm=(G*G).sum(1)
dots=G@G.T; np.fill_diagonal(dots,0)
assert np.all(4*np.maximum(dots,0)**2<=norm[:,None]*norm)
parity_counts=[0,0]
for m in v:
 ix=[i for i in range(N) if m>>i&1]
 # A wrong parity costs the smallest sign change, zero on a missing site.
 X=G[:,ix]; Q=2*abs(X).sum(1); pp=(X<0).sum(1)%2
 cand=[]
 for p in (0,1):
  M=Q-4*abs(X).min(1)*(pp!=p)
  if np.all(4*M*M<=32*norm):cand.append(p)
 assert cand, tuple(ix)
 parity_counts[cand[0]]+=1
# Common-parity signings on one octad differ at >=2 signs; different octads
# intersect in <=4 sites. Any dense/ordinary-octad product is at most 16.
words=[0]
for a in map(int,C['dense_generator']):
 assert 0<a<1<<32
 words.extend([w^a for w in words])
assert len(words)==len(set(words))==131072
assert min(w.bit_count() for w in words[1:])>=8
shift=int(C['dense_coset_shift']); assert 0<=shift<1<<32
assert (abs(S[:,D]).sum(1)<=16).all() and (abs(R[:,D]).sum(1)<=16).all()
for block in np.array_split(np.asarray(words,dtype=np.int64)^shift,32):
 dense=1-2*((block[:,None]>>np.arange(32))&1)
 q=H[:,D]@dense.T
 assert np.all(4*np.maximum(q,0)**2<=32*(H*H).sum(1)[:,None])
# Support cardinalities 32,8,5,2,35 preclude duplicate rays across layers.
total=len(words)+128*len(v)+len(S)+len(R)+len(H)
assert total==415360
print(json.dumps(dict(status='PASS',dimension=N,total=total,octad_parity_counts=parity_counts,stars=len(S),roots=len(R),heads=len(H))))
