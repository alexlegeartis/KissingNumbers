from pathlib import Path
import numpy as np,collections
p=Path(__file__).parent;z=np.load(p/'source.npz');X=z['X'].astype(np.int64);g=z['g'].astype(np.int64);assert X.shape==(2564,15) and np.array_equal(g,[1]*14+[2]) and np.all((X*X*g).sum(1)==16)
for st in range(0,len(X),128):
 D=(X[st:st+128]*g)@X.T;ii=np.arange(st,st+len(D))[:,None];jj=np.arange(len(X))[None,:];assert not np.any((D>8)&(ii<jj))
a=X[:,:14];t=X[:,14];E=(t==0)&(np.count_nonzero(a,axis=1)==4)&np.isin(a,[-2,0,2]).all(1);H=(np.abs(a)==1).all(1)&(np.abs(t)==1);L=(np.count_nonzero(a,axis=1)==2)&np.isin(a,[-2,0,2]).all(1)&(np.abs(t)==2);T=(t==0)&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==4).any(1);assert (E.sum(),(H&(t==1)).sum(),(H&(t==-1)).sum(),(L&(t==2)).sum(),(L&(t==-2)).sum(),T.sum())==(1456,512,512,28,28,28) and np.all(E.astype(int)+H+L+T==1)
print('Verified all 2,564 source vectors and pairs.')
z=np.load(p/'low_points.npz');Q=z['Q'].astype(np.int64);B=z['blockers'];part=z['part'];assert len(Q)==17572 and np.array_equal(np.bincount(part,minlength=2),[2212,15360]) and np.all((Q*Q*g).sum(1)==16)
for st in range(0,len(Q),128):
 D=(Q[st:st+128]*g)@X.T
 for i,row in enumerate(D):assert set(map(int,np.flatnonzero(row>8)))==set(map(int,B[st+i][B[st+i]!=65535]))
old={tuple(v) for v in X};isold=np.array([tuple(v) in old for v in Q]);assert isold.sum()==2564 and collections.Counter((B[~isold]!=65535).sum(1))=={4:672,6:14336};print('Verified 17,572 saved shell points and exact blocker lists.')
z=np.load(p/'pool.npz');v=np.load(p/'sparse_low.npz');Q=np.r_[z['Q'],v['Q']].astype(np.int64);B=np.r_[z['B'],v['blockers']];_,idx=np.unique(Q,axis=0,return_index=True);Q=Q[idx];B=B[idx];assert len(Q)==52373 and np.all((Q*Q*g).sum(1)==16*144) and len({tuple(v) for v in Q})==len(Q)
for st in range(0,len(Q),128):
 D=(Q[st:st+128]*g)@(12*X).T
 for i,row in enumerate(D):assert set(map(int,np.flatnonzero(row>8*144)))==set(map(int,B[st+i][B[st+i]!=65535]))
a=Q[:,:14];t=Q[:,14];E=(t==0)&(np.count_nonzero(a,axis=1)==4)&np.isin(a,[-24,0,24]).all(1);H=(np.abs(a)==12).all(1)&(np.abs(t)==12);L=(np.abs(t)==24);T=(t==0)&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==48).any(1);assert np.count_nonzero(~(E|H|L|T))==5120
print('Verified all 52,373 search candidates and their exact source blocker lists.')

z=np.load(p/'sparse_low.npz');Q=z['Q'].astype(np.int64);B=z['blockers'];part=z['part'];assert len(Q)==5120 and np.array_equal(np.bincount(part,minlength=2),[1536,3584]) and np.all((Q*Q*g).sum(1)==16*144)
for st in range(0,len(Q),128):
 D=(Q[st:st+128]*g)@(12*X).T
 for i,row in enumerate(D):assert set(map(int,np.flatnonzero(row>8*144)))==set(map(int,B[st+i][B[st+i]!=65535]))
assert collections.Counter((B!=65535).sum(1))=={8:768,9:1280,10:3072}
print('Verified all 5,120 outside-family sparse-head candidates and exact blockers.')
