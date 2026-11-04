"""Export a nonlinear recovered 31560-point configuration and check every pair."""
import json,time
import numpy as np
from model import B,C,N,DATA,labels_from_holes,verify_labels

start=time.monotonic()
holes=np.load(DATA/'nonlinear_quotient_holes.npy')[0]
labels=labels_from_holes(holes)
global_labels=np.full(4096,-2,dtype=int)
cosets=[]
for m in range(4096):
    if global_labels[m]!=-2:continue
    ids=C^m
    assert np.all(global_labels[ids]==-2)
    global_labels[ids]=labels
    cosets.append(m)
assert len(cosets)==4 and np.sum(global_labels>=0)==3840
ids=np.flatnonzero(global_labels>=0)
extra=B['profiles'][global_labels[ids]]*(1-2*B['code'][ids,:21])
Y=np.concatenate([B['base'],extra]).astype(np.int64)
assert Y.shape==(31560,21) and len(np.unique(Y,axis=0))==len(Y)
Q=(Y*Y).sum(1)
F=Y.astype(np.float64)
violations=0
for a in range(0,len(Y),256):
    dots=F[a:a+256]@F.T
    assert np.array_equal(dots,np.rint(dots))
    D=dots.astype(np.int64)
    bad=(D>0)&(4*D*D>Q[a:a+256,None]*Q[None,:])
    for k in range(len(D)):bad[k,a+k]=False
    violations+=int(bad.sum())
assert violations==0
# Real negative control: replacing a selected word's label by the direction
# to another selected word must be rejected by the discrete verifier.
corrupt=labels.copy()
for u in np.flatnonzero(labels>=0):
    candidates=np.flatnonzero(labels[N[u]]>=0)
    if len(candidates):corrupt[u]=candidates[0];break
try: verify_labels(corrupt)
except AssertionError: negative_control='PASS'
else: raise AssertionError('invalid labeling was accepted')
# The geometric predicate must also detect a duplicated ray.
assert 4*int(Y[0]@Y[0])**2>int(Q[0])**2
np.savez_compressed(DATA/'nonlinear21_verified.npz',rays=Y.astype(np.int16),
                    labels=global_labels,component_holes=holes)
report=dict(status='PASS',dimension=21,count=len(Y),base=27720,auxiliary=3840,
            unordered_pairs=len(Y)*(len(Y)-1)//2,violations=violations,
            discrete_negative_control=negative_control,geometric_duplicate_control='PASS',
            arithmetic='exact small integer dot products via float64 BLAS; integer squared comparisons',
            seconds=time.monotonic()-start,coset_representatives=cosets,
            novelty_claim=False)
(DATA/'geometry_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
