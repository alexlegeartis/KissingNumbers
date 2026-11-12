"""Collect the bounded six-term D7 extensions and test root differences for D8."""
from pathlib import Path
from collections import defaultdict
import itertools
import json
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
raw=json.loads((REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json').read_text())
B=np.asarray(raw['catalogue_gram'],dtype=np.int64)
U=json.loads((HERE/'d7_embedding.json').read_text())['normals'][:6]
M=np.asarray(U,dtype=np.int64)@B
cols=[tuple(map(int,M[:,i])) for i in range(72)]
target=(0,0,0,0,0,-4)
triples=defaultdict(list)
for i,j,k in itertools.combinations(range(72),3):
    a,b,c=cols[i],cols[j],cols[k]
    base=((i*72+j)*72+k)*8
    for signs in range(8):
        si,sj,sk=[1 if signs&(1<<x) else -1 for x in range(3)]
        triples[tuple(si*x+sj*y+sk*z for x,y,z in zip(a,b,c))].append(base+signs)
def unpack(code):
    code,s=divmod(code,8);code,k=divmod(code,72);i,j=divmod(code,72)
    return [(i,1 if s&1 else -1),(j,1 if s&2 else -1),(k,1 if s&4 else -1)]
pool=set();checked=0
for key,aa in triples.items():
    other=tuple(t-x for t,x in zip(target,key))
    if other<key or other not in triples:continue
    for a,b in itertools.product(aa,triples[other]):
        sparse={}
        for i,s in unpack(a)+unpack(b):sparse[i]=sparse.get(i,0)+s
        sparse={i:s for i,s in sparse.items() if s}
        norm=sum(si*sj*int(B[i,j]) for i,si in sparse.items() for j,sj in sparse.items())
        checked+=1
        if norm==8:pool.add(tuple(sparse.get(i,0) for i in range(72)))
pool=sorted(pool)
print('D7 extensions',len(pool),'tested matching pairs',checked,flush=True)
(HERE/'d7_normal_pool.json').write_text(json.dumps(dict(normals=pool,finite_scope='Sum of two signed, distinct-coordinate triples.'),indent=2)+'\n')
vectors=np.asarray(pool,dtype=np.int64)
products=vectors@B@vectors.T
found=False
for i in range(len(pool)):
    for j in range(i):
        if products[i,j]==4:
            v=vectors[i];w=vectors[j]-v
            ext=U+[v.tolist(),w.tolist()]
            gram=np.asarray(ext)@B@np.asarray(ext).T
            assert (gram.diagonal()==8).all()
            assert gram[-1,:-1].tolist()==[0,0,0,0,0,0,-4]
            (HERE/'d8_embedding.json').write_text(json.dumps(dict(normals=ext,gram=gram.tolist()),indent=2)+'\n')
            print('D8 from pair',i,j,'supports',sum(v!=0),sum(w!=0),flush=True)
            found=True;break
    if found:break
print('D8 difference found',found,flush=True)
