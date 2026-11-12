"""Bounded meet-in-the-middle search for exact D_k normal extensions."""
from pathlib import Path
from collections import defaultdict
import itertools
import json
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
source = REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json'
d = json.loads(source.read_text())
B = np.asarray(d['catalogue_gram'], dtype=np.int64)
U = list(d['normals'])+[json.loads((HERE/'d6_normal.json').read_text())['normal']]

def unpack(code):
    code, signs = divmod(code, 8)
    code, k = divmod(code, 72)
    i, j = divmod(code, 72)
    return [(i, 1 if signs&1 else -1), (j, 1 if signs&2 else -1), (k, 1 if signs&4 else -1)]

for rank in range(7, 10):
    M = np.asarray(U, dtype=np.int64)@B
    columns = [tuple(int(x) for x in M[:,i]) for i in range(72)]
    target = (0,)*(rank-2)+(-4,)
    triples = defaultdict(list)
    for i,j,k in itertools.combinations(range(72),3):
        a,b,c = columns[i],columns[j],columns[k]
        base = ((i*72+j)*72+k)*8
        for signs in range(8):
            si,sj,sk = [1 if signs&(1<<x) else -1 for x in range(3)]
            key = tuple(si*x+sj*y+sk*z for x,y,z in zip(a,b,c))
            triples[key].append(base+signs)
    print('D',rank,'triple signatures',len(triples),flush=True)
    found = None
    checked = 0
    for key, aa in triples.items():
        other = tuple(t-x for t,x in zip(target,key))
        if other < key or other not in triples:
            continue
        for a,b in itertools.product(aa,triples[other]):
            sparse = {}
            for i,s in unpack(a)+unpack(b):
                sparse[i] = sparse.get(i,0)+s
            sparse = {i:s for i,s in sparse.items() if s}
            norm = sum(si*sj*int(B[i,j]) for i,si in sparse.items() for j,sj in sparse.items())
            checked += 1
            if norm == 8:
                found = [sparse.get(i,0) for i in range(72)]
                assert np.array_equal(M@found,target)
                break
        if found is not None:
            break
    print('D',rank,'checked',checked,'normal',None if found is None else [(i,s) for i,s in enumerate(found) if s],flush=True)
    if found is None:
        break
    U.append(found)
    G = np.asarray(U,dtype=np.int64)@B@np.asarray(U,dtype=np.int64).T
    (HERE/f'd{rank}_embedding.json').write_text(json.dumps(dict(normals=U,gram=G.tolist()),indent=2)+'\n')
