from pathlib import Path
Path('research/mixed_tier17_joint_layers_20261105/cache/').mkdir(parents=True,exist_ok=True)
import numpy as np
N=15872;P=2*N;pos=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/baseline_seed.bin',np.int32);neg=np.loadtxt('research/mixed_tier17_order8_20261105/cache/negative_seed.txt',dtype=np.int32);base=np.array([4*int(v)+j for v in pos for j in (0,1,2,3)]+[N+4*int(v)+j for v in neg for j in (0,1,2,3)],np.int32);flag=np.zeros(P,bool);flag[base]=1;cost=np.zeros(P,np.int32);E=np.fromfile('research/mixed_tier17_order4_20261105/cache/graph_edges.bin',np.int32).reshape(-1,2);X=np.fromfile('research/mixed_tier17_order4_20261105/cache/cross_edges.bin',np.int32).reshape(-1,2)
def accumulate(a,b):
 i=flag[a]&~flag[b];np.add.at(cost,b[i],1);j=flag[b]&~flag[a];np.add.at(cost,a[j],1)
for off in (0,N):accumulate(E[:,0]+off,E[:,1]+off)
accumulate(X[:,0],X[:,1]+N)
cost[flag]=999;vals,cts=np.unique(cost,return_counts=True);print('blocker histogram',list(zip(vals[:18],cts[:18])),flush=True);np.savez_compressed('research/mixed_tier17_joint_layers_20261105/cache/cap_pool.npz',base=base,cost=cost)
