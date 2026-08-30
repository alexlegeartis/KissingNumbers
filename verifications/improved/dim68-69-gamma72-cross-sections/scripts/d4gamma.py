#!/usr/bin/env python3
"""The D_4-patterned four-point Gram of Gamma_72 (diag 8, three edges -4, det 1024).

On the Leech the D_4 cross-section is Lambda_20 = 17400 and the A_4 one only 15540, and the
four-point LP reproduces both exactly, so D_4 is the pattern to aim for.  Scaled to Gamma_72
its determinant is 1024, well below the 1700 reachable by extending the (4,-2,0) triple."""
import numpy as np, os, itertools
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.join(HERE,'..','..','dim73-95-gamma72-caps','data')
if not os.path.exists(os.path.join(D, 'gamma72_pool.npy')):
    raise SystemExit(
        "This is the SEARCH that found the D_4 four-tuple, not a verifier.  It needs a "
        "pool of Gamma_72 minimal vectors (488 MB), which this repository does not ship. "
        "What it found IS shipped, and scripts/verify68.py checks it from scratch using "
        "nothing but the Gram.")
G=np.load(os.path.join(D,'gamma72_gram.npy')).astype(np.int64)
P=np.load(os.path.join(D,'gamma72_pool.npy'),mmap_mode='r')
X=np.asarray(P[:900000],dtype=np.int64); XG=X@G
print("pool sample %d"%len(X),flush=True)
rng=np.random.default_rng(0)
found=None
for t in range(4000):
    ci=int(rng.integers(len(X)))
    ipc=XG@X[ci]
    nb=np.nonzero(ipc==-4)[0]
    if len(nb)<3: continue
    S=X[nb]; M=S@G@S.T
    ok=np.argwhere(np.triu(M==0,1))
    for a,b in ok:
        for c3 in range(len(nb)):
            if c3 in (a,b): continue
            if M[a,c3]==0 and M[b,c3]==0:
                found=np.vstack([X[ci],S[a],S[b],S[c3]]); break
        if found is not None: break
    if found is not None: break
if found is None:
    print("no D_4 pattern found in this sample"); raise SystemExit
Mt=found@G@found.T
print("D_4 pattern FOUND, exact Gram:")
for r in Mt.tolist(): print("   ",r)
print("determinant %d  (target 1024)"%int(round(np.linalg.det(Mt.astype(float)))))
np.save(os.path.join(HERE,'gamma72_d4_k4.npy'),found)
# also try A_4: a chain u1-u2-u3-u4 with -4 on the chain and 0 elsewhere
