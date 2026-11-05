import numpy as np
p='research/mixed_tier17_order8_20261105/cache/';q='research/mixed_tier17_eval_subgroup_20261105/cache/';N=7936
E=np.fromfile(p+'graph_edges.bin',np.int32).reshape(-1,2);X=np.fromfile(p+'cross_edges.bin',np.int32).reshape(-1,2)
with open(p+'anneal_edges.bin','wb') as f:
 E.tofile(f);(E+N).tofile(f);(X+np.array([0,N],np.int32)).tofile(f)
pos=np.fromfile(q+'baseline_seed.bin',np.int32);neg=np.loadtxt(p+'negative_seed.txt',dtype=np.int32)
base=[2*int(v)+j for v in pos for j in (0,1)]+[N+2*int(v)+j for v in neg for j in (0,1)];assert len(base)==176
origin=np.zeros(2*N,np.uint8);origin[base]=1;origin.tofile(p+'anneal_origin.bin');print('seed',len(base),'graph edges',2*len(E)+len(X))
