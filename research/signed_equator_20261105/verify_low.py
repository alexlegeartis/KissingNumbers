from pathlib import Path
import numpy as np,collections
root=Path(__file__).parent;z=np.load(root/'source.npz');X=np.c_[z['rational'],z['radical']].astype(np.int64);z=np.load(root/'low_points.npz');Q=z['Q'].astype(np.int64);B=z['blockers'];part=z['part'];assert len(Q)==2762 and np.array_equal(np.bincount(part,minlength=2),[1274,1488]);assert np.all((Q[:,:13]**2+3*Q[:,13:]**2).sum(1)==16) and np.all((Q[:,:13]*Q[:,13:]).sum(1)==0)
for st in range(0,len(Q),128):
 p=Q[st:st+128,:13]@X[:,:13].T+3*Q[st:st+128,13:]@X[:,13:].T-8;q=Q[st:st+128,:13]@X[:,13:].T+Q[st:st+128,13:]@X[:,:13].T;bad=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q))
 for i,row in enumerate(bad):assert set(map(int,np.flatnonzero(row)))==set(map(int,B[st+i][B[st+i]!=65535]))
S={tuple(v) for v in X};old=np.array([tuple(v) in S for v in Q]);assert old.sum()==1154;assert collections.Counter((B[~old]!=65535).sum(1))=={2:264,4:1344};a=Q[:,:13];b=Q[:,13:];z=(b==0).all(1);heads=z&(np.abs(a[:,:12])==1).all(1)&(np.abs(a[:,12])==2);eq=z&(a[:,12]==0)&(np.count_nonzero(a[:,:12],axis=1)==4)&np.isin(a[:,:12],[-2,0,2]).all(1);axes=z&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==4).any(1);assert not np.any(~old&~heads&~eq&~axes)
print('Verified all 2,762 saved points, exact source blockers, and low-cost family coverage.')
