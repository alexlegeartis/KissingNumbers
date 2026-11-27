#!/usr/bin/env python3
"""Construct and independently check the fixed-stock root-star matching exchange."""
import hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

here=Path(__file__).resolve().parent
source=Path(sys.argv[1])
assert hashlib.sha256(source.read_bytes()).hexdigest()=="90d130ba1bb96172fae7b7ae10d1169747b76f2e720946bc2350c99a08d427c9"
with np.load(source,allow_pickle=False) as f:
 data=dict(f)
s=data['sparse'].astype(np.int32)
w=s[np.count_nonzero(s,axis=1)==5]
o=s[np.count_nonzero(s,axis=1)==8]
assert w.shape==(172,33) and o.shape==(2216,33)
r=[]
for i,j in itertools.combinations(range(33),2):
 for x,y in itertools.product((-4,4),repeat=2):
  a=np.zeros(33,np.int32);a[i]=x;a[j]=y;r.append(a)
r=np.array(r)
er=w@r.T>16
assert (er.sum(axis=1)==4).all()
g=csr_matrix(er)
m=maximum_bipartite_matching(g,perm_type='column')
m1=m>=0;assert np.count_nonzero(m1)==152 and len(np.unique(m[m1]))==152
assert er[np.flatnonzero(m1),m[m1]].all()
back=np.full(len(r),-1)
for i,j in enumerate(m):
 if j>=0:back[j]=i
S=set(np.flatnonzero(~m1));R=set();stack=list(S)
while stack:
 i=stack.pop()
 for j in g.indices[g.indptr[i]:g.indptr[i+1]]:
  if j in R:continue
  R.add(int(j));k=int(back[j])
  if k>=0 and k not in S:S.add(k);stack.append(k)
assert len(S)==132 and len(R)==112
present=np.ones(len(r),bool);present[sorted(R)]=False
new=np.vstack((o,w[sorted(S)],r[present])).astype(np.int16)
G=new.astype(np.int32)@new.astype(np.int32).T
np.fill_diagonal(G,-32)
assert G.max()<=16
data['sparse']=new
np.savez_compressed(here/'witness.npz',**data)
h=hashlib.sha256((here/'witness.npz').read_bytes()).hexdigest()
out=dict(dimension=33,count=370556,parent_count=370528,witness_sha256=h)
(here/'manifest.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(dict(matching_size=152,selected_stars=len(S),blocked_roots=len(R),
                      signed_octads=len(o),retained_roots=int(present.sum()),witness_sha256=h),indent=2))
