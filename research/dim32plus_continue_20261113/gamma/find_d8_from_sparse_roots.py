"""Try sums/differences of the finite signed four-coordinate root pool."""
from pathlib import Path
from collections import defaultdict
import json
import subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
raw=json.loads((REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json').read_text())
B=np.asarray(raw['catalogue_gram'],dtype=np.int64)
U=np.asarray(json.loads((HERE/'d7_embedding.json').read_text())['normals'],dtype=np.int64)
M=U@B
gramfile=HERE/'root_pool_gram.txt';poolfile=HERE/'sparse_root_pool.txt'
gramfile.write_text('\n'.join(' '.join(map(str,r)) for r in B)+'\n')
subprocess.run([str(HERE/'sparse_root_pool'),str(gramfile),str(poolfile)],check=True)
roots=[]
for line in poolfile.read_text().splitlines():
    sparse=tuple((abs(x)-1,1 if x>0 else -1) for x in map(int,line.split()) if x)
    roots.append(sparse)
    roots.append(tuple((i,-s) for i,s in sparse))
profiles=defaultdict(list)
for sparse in roots:
    profile=tuple(sum(s*int(M[j,i]) for i,s in sparse) for j in range(7))
    profiles[profile].append(sparse)
print('Root pool',len(roots),'profile buckets',len(profiles),flush=True)
target=(0,0,0,0,0,0,-4)
found=None;checked=0
for key,aa in profiles.items():
    other=tuple(t-x for t,x in zip(target,key))
    if other<key or other not in profiles:continue
    for a in aa:
        for b in profiles[other]:
            inner=sum(sa*sb*int(B[i,j]) for i,sa in a for j,sb in b)
            checked+=1
            if inner==-4:
                v=np.zeros(72,dtype=np.int64)
                for i,s in a+b:v[i]+=s
                assert v@B@v==8 and tuple(M@v)==target
                found=v.tolist();break
        if found is not None:break
    if found is not None:break
print('D8 compatible pairs tested',checked,'found',None if found is None else [(i,s) for i,s in enumerate(found) if s],flush=True)
if found is not None:
    ext=U.tolist()+[found]
    gram=np.asarray(ext)@B@np.asarray(ext).T
    (HERE/'d8_embedding.json').write_text(json.dumps(dict(normals=ext,gram=gram.tolist()),indent=2)+'\n')
