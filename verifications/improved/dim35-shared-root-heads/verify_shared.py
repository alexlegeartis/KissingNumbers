#!/usr/bin/env python3
"""Exact independent verifier, no importer from the construction search."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import itertools as I,json
from pathlib import Path
from hashlib import sha256
import numpy as np
if not __debug__:
    raise SystemExit("Run without -O: verification assertions must be enabled.")
P=Path(__file__).resolve().parent
c=json.loads((P/'code_data.json').read_text())
S=[int(a) for a,_ in c['full_parity_supports']]
z=np.load(P/'shared_positive.npz',allow_pickle=False)
def ints(k):
 x=z[k]
 assert x.dtype.kind in 'iu'
 return x.astype(np.int64)
V,H,D,R=map(ints,['new_stars','heads','dense_positions','removed_roots'])
assert len(D)==len(set(map(int,D)))==32 and 0<=D.min()<=D.max()<35
assert V.ndim==H.ndim==R.ndim==2 and V.shape[1]==H.shape[1]==R.shape[1]==35
assert (np.abs(V)<=4).all() and (np.abs(H)<=3).all() and np.all(np.isin(R,[-4,0,4]))
assert (np.sort(abs(V),axis=1)==np.array([0]*30+[2]*4+[4])).all()
assert (np.sort(abs(H),axis=1)==np.array([1]*33+[3]*2)).all()
assert np.all(np.abs(H[:,D])==1)
assert len({tuple(x) for x in V})==len(V)
assert len({tuple(x) for x in H})==len(H)
assert len({tuple(x) for x in R})==len(R)
A=[]
for a,b in I.combinations(range(35),2):
 for u,v in I.product([-4,4],repeat=2):
  x=[0]*35;x[a]=u;x[b]=v;A.append(x)
A=np.array(A,dtype=np.int64);mV=(V@A.T>16).any(0);mH=(H@A.T>20).any(0)
assert {tuple(x)for x in A[mV|mH]}=={tuple(x)for x in R}
U=A[~(mV|mH)]
def maxima(X,Y,nX,nY,same=False):
 w=X@Y.T
 if same:np.fill_diagonal(w,0)
 m=int(w.max());assert 4*max(m,0)**2<=nX*nY
 return m
pair={'star_star':maxima(V,V,32,32,True),'head_head':maxima(H,H,51,51,True),'star_head':maxima(V,H,32,51),'star_root':maxima(V,U,32,32),'head_root':maxima(H,U,51,32)}
# Dense words encoded by parity-check independent generation; min distance is a full exhaustive nonzero difference check.
G=[int(x)for x in c['dense_generator']]
B=[0]
for g in G:B+=[x^g for x in B]
assert len(B)==len(set(B))==131072 and min(x.bit_count()for x in B[1:])>=8
B=np.array(B,dtype=np.int64)^int(c['dense_coset_shift']);B=1-2*((B[:,None]>>np.arange(32))&1)
pair['head_dense']=maxima(H[:,D],B,51,32)
assert abs(V[:,D]).sum(1).max()<=16
# All other root/dense and octad/dense dot bounds follow from magnitude; root/root is exact.
assert len(S)==len(set(S))==2202
assert all(len(set([j for j in range(35)if s>>j&1]))==8 and 0<s<1<<35 for s in S)
assert max((a&b).bit_count()for a,b in I.combinations(S,2))<=4
pd=[]
for s in S:
 ix=[i for i in range(35)if s>>i&1]
 assert 2*abs(V[:,ix]).sum(1).max()<=16
 h=H[:,ix];dots=2*abs(h).sum(1);mins=4*abs(h).min(1);neg=(h<0).sum(1)&1
 p=next((p for p in range(2)if (dots-mins*(neg!=p)).max()<=20),None)
 assert p is not None
 pd.append([s,p])
# Distinct shape categories, because dense all have support 32, heads support 35,
# middle support 8, stars 5, roots 2. Exact whole count.
total=len(B)+128*len(S)+len(U)+len(V)+len(H)
assert total==415336,total
report={'status':'PASS','dimension':35,'points':total,'gain':total-415334,'root_star_count':len(V),'head_count':len(H),'delta_deleted_roots':len(R),'roots_retained':len(U),'roots_paid_by_stars':int(mV.sum()),'roots_paid_by_heads':int(mH.sum()),'root_overlap':int((mV&mH).sum()),'extra_coordinates':sorted(set(range(35))-set(map(int,D))),'pair_maxima_checked':pair,'certificate_sha256':sha256((P/'shared_positive.npz').read_bytes()).hexdigest()}
(P/'verified.json').write_text(json.dumps(report,indent=2)+'\n')
(P/'full_parities.json').write_text(json.dumps(pd,separators=(',',':'))+'\n')
print(json.dumps(report,indent=2))
