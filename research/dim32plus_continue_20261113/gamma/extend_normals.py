"""Find sparse, exact D_k extensions in the supplied Gamma72 Gram basis."""
from pathlib import Path
import itertools
import json
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json'
d = json.loads(SOURCE.read_text())
B = np.asarray(d['catalogue_gram'], dtype=np.int64)
U = [list(x) for x in d['normals']]
U.append(json.loads((HERE/'d6_normal.json').read_text())['normal'])
for k in range(7, 13):
    profiles = np.asarray(U, dtype=np.int64) @ B
    target = np.zeros(k-1, dtype=np.int64)
    target[-1] = -4
    pairs = {}
    for i in range(72):
        for j in range(i):
            for si, sj in itertools.product((-1, 1), repeat=2):
                signature = tuple(si*profiles[:, i]+sj*profiles[:, j])
                pairs.setdefault(signature, []).append(((i, si), (j, sj)))
    found = None
    checked = 0
    for signature, aa in pairs.items():
        bb = pairs.get(tuple(target - np.asarray(signature)))
        if bb is None:
            continue
        for a, b in itertools.product(aa, bb):
            sparse = {}
            for i, s in a+b:
                sparse[i] = sparse.get(i, 0) + s
            sparse = {i:s for i,s in sparse.items() if s}
            norm = sum(si*sj*int(B[i,j]) for i,si in sparse.items() for j,sj in sparse.items())
            checked += 1
            if norm == 8:
                v = [sparse.get(i,0) for i in range(72)]
                assert np.array_equal(profiles@v, target)
                found = v
                break
        if found is not None:
            break
    print('D', k, 'checked', checked, 'found', None if found is None else [(i,x) for i,x in enumerate(found) if x], flush=True)
    if found is None:
        break
    U.append(found)
    G = np.asarray(U, dtype=np.int64)@B@np.asarray(U, dtype=np.int64).T
    (HERE/f'd{k}_embedding.json').write_text(json.dumps(dict(normals=U, gram=G.tolist()), indent=2)+'\n')
