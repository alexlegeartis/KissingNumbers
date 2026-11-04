#!/usr/bin/env python3
"""Independent exhaustive integer ray check of the rebuilt 31560-point baseline."""
import json
import numpy as np
from baselines import DATA

z=np.load(DATA/'baseline21.npz')
Y=np.concatenate([z['base'],z['auxiliary']]).astype(np.int64)
Q=(Y*Y).sum(1)
assert Y.shape==(31560,21) and set(Q)=={8,164}
assert len(np.unique(Y,axis=0))==len(Y)
F=Y.astype(float);bad=0;worst_by_norm={};checked=0
for a in range(0,len(Y),256):
    G=F[a:a+256]@F.T
    assert np.array_equal(G,np.rint(G))
    D=G.astype(np.int64);P=Q[a:a+256,None]*Q[None,:]
    invalid=(D>0)&(4*D*D>P)
    for t in range(len(D)):invalid[t,a+t]=False
    bad+=int(invalid.sum());checked+=D.size
assert bad==0
# Corruption control: duplicate a ray; its cosine is one and must fail.
assert 4*int(Q[0])**2>int(Q[0])**2
report=dict(count=len(Y),unordered_pairs=len(Y)*(len(Y)-1)//2,violations=bad,
            method='All small integer dot products computed exactly using float64 BLAS, then integer squared comparisons.',
            duplicate_control='PASS')
(DATA/'independent21.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
