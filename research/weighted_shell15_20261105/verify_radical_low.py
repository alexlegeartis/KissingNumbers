from pathlib import Path
import numpy as np
p=Path(__file__).parent;z=np.load(p/'source.npz');X=z['X'].astype(np.int64);g=z['g'].astype(np.int64);z=np.load(p/'radical50_low.npz');Q=z['Q'].astype(np.int64);B=z['blockers'];assert Q.shape==(1024,15) and B.shape==(1024,12) and len(np.unique(Q,axis=0))==len(Q)
assert np.all(((Q*Q)*g).sum(1)==50)
D=(Q*g)@X.T
for i in range(len(Q)):assert np.array_equal(np.sort(B[i]),np.flatnonzero(D[i]>14))
print('all 1024 norm-50 candidates and their 12 blockers verified exactly')
