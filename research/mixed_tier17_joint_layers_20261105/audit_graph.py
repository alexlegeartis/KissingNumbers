import numpy as np
z=np.load('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz');nodes=np.load('research/mixed_tier17_joint_layers_20261105/cache/joint_nodes.npz');X=z['X'].astype(np.int32);Q=z['points'].reshape(-1,4,17).astype(np.int32);eq=nodes['eqids'];cap=nodes['caprep'].astype(np.int32);ng=int(nodes['ng']);group=z['point_group'];H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];sg=np.array([[1-2*(h>>j&1) for j in range(16)]+[1] for h in K],np.int32);g=np.r_[np.ones(16,dtype=np.int32),2];E=np.fromfile('research/mixed_tier17_joint_layers_20261105/cache/joint_edges.bin',np.int32).reshape(-1,2);N=ng+len(eq)+len(cap);a=np.minimum(E[:,0],E[:,1]);b=np.maximum(E[:,0],E[:,1]);codes=np.sort(a.astype(np.int64)*N+b);rng=np.random.default_rng(170611)
def expand(i):
 if i<ng:return X[group==i]
 i-=ng
 if i<len(eq):return Q[eq[i]]
 return cap[i-len(eq)]*sg
for t in range(2000):
 if t<1000:i,j=map(int,E[rng.integers(len(E))])
 else:i,j=map(int,rng.choice(N,2,replace=False))
 A=expand(i);B=expand(j);bad=bool(np.max((A*g)@B.T)>144);key=min(i,j)*N+max(i,j);at=np.searchsorted(codes,key);found=at<len(codes) and codes[at]==key;assert bad==found,(i,j,bad,found)
print('PASS 2000 source/equator/cap subgroup Gram comparisons')
