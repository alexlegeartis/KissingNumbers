"""Compare sampled orbit conflicts with direct exact 17D Gram checks."""
import numpy as np
p='research/mixed_tier17_equator_exchange_20261105/cache/';z=np.load(p+'pool.npz');X=z['X'].astype(np.int32);Q=z['points'].reshape(-1,4,17).astype(np.int32);cand=np.load(p+'full_candidate_ids.npy');ng=len(z['weights']);point_group=z['point_group'];E=np.fromfile(p+'full_edges.bin',np.int32).reshape(-1,2);N=ng+len(cand);a=np.minimum(E[:,0],E[:,1]);b=np.maximum(E[:,0],E[:,1]);codes=np.sort(a.astype(np.int64)*N+b);rng=np.random.default_rng(171014);g=np.r_[np.ones(16,dtype=np.int32),2]
def expand(i):return X[point_group==i] if i<ng else Q[cand[i-ng]]
# Half graph edges and half unrestricted pairs, to exercise both decisions.
for t in range(2000):
 if t<1000:i,j=map(int,E[rng.integers(len(E))])
 else:i,j=map(int,rng.choice(N,2,replace=False))
 A=expand(i);B=expand(j);bad=bool(np.max((A*g)@B.T)>144);key=min(i,j)*N+max(i,j);at=np.searchsorted(codes,key);found=at<len(codes) and codes[at]==key
 assert bad==found,(i,j,bad,found)
print('PASS 2000 exact weighted-orbit Gram comparisons')
