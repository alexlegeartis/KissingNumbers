"""Verify every saved low-blocker vector of the complete six-support scan."""
from pathlib import Path
import numpy as np
p=Path(__file__).parent;z=np.load(p/'coordinates.npz');X=z['switched_twice'].astype(np.int64);assert np.all(X%18==0);X//=18;z=np.load(p/'tier_low.npz');Q=z['Q'].astype(np.int64);B=z['blockers'];assert Q.shape==(32768,17) and B.shape==(32768,12) and len(np.unique(Q,axis=0))==len(Q)
assert np.all((Q[:,:16]**2).sum(1)+2*Q[:,16]**2==288)
hist=np.bincount((B!=65535).sum(1),minlength=13);assert hist.tolist()==[0,1024,0,0,0,0,0,0,9216,4096,11264,4096,3072]
for st in range(0,len(Q),128):
 W=Q[st:st+128];D=W[:,:16]@X[:,:16].T+2*np.outer(W[:,16],X[:,16]);bad=D>144
 for j,row in enumerate(B[st:st+128]):
  ids=np.sort(row[row!=65535]);assert np.array_equal(ids,np.flatnonzero(bad[j]))
  if len(ids)==1:assert np.array_equal(W[j],X[ids[0]])
print('all 32768 low outputs verified, including 1024 source duplicates; no free extension among them')
