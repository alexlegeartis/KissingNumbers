from pathlib import Path
import numpy as np
p=Path(__file__).parent;X=np.fromfile(p/'search16_switched_twice_baseline.bin',np.int8).reshape(-1,16);Q=np.fromfile(p/'search16_switched_eight_Q.bin',np.int8).reshape(-1,16);B=np.fromfile(p/'search16_switched_eight_B.bin',np.uint16).reshape(-1,8);assert len(X)==4320 and len(Q)==len(B)==124416 and np.all((Q.astype(np.int32)**2).sum(1)==8);k=(B!=65535).sum(1);assert (k==1).sum()==3712 and (k==8).sum()==120704;ids=np.flatnonzero(k==1);assert np.array_equal(Q[ids],X[B[ids,0]])
rng=np.random.default_rng(160752)
for i in rng.choice(len(Q),1000,False):
 d=Q[i].astype(np.int32)@X.astype(np.int32).T;assert np.array_equal(np.sort(B[i][B[i]!=65535]),np.flatnonzero(d>4))
print('all retained norms and source duplicates, plus 1000 sampled blocker lists, verified')
